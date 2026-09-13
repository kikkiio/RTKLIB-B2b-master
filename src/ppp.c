/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS tide model
*           -DOUTSTAT_AMB output ambiguity parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to rtkcmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
*           2015/05/10 1.8  add function to detect slip by MW-LC jump
*                           fix ppp solution problem with large clock variance
*           2015/06/08 1.9  add precise satellite yaw-models
*                           cope with day-boundary problem of satellite clock
*           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
*                           pppoutsolsat() -> pppoutstat()
*           2015/11/13 1.11 add L5-receiver-dcb estimation
*                           merge post-residual validation by rnx2rtkp_test
*                           support support option opt->pppopt=-GAP_RESION=nnnn
*           2016/01/22 1.12 delete support for yaw-model bug
*                           add support for ura of ephemeris
*           2018/10/10 1.13 support api change of satexclude()
*           2020/11/30 1.14 use sat2freq() to get carrier frequency
*                           use E1-E5b for Galileo iono-free LC
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */

/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?((opt)->if_model==1?2:1):(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)

#ifdef BDS2BDS3
#define NC(opt)     (NSYS+1)
#else
#define NC(opt)     (NSYS)
#endif

#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->if_model==1?(NSATGPS+1+(((opt)->navsys&SYS_GAL)?NSATGAL:0)):((opt)->if_model?0:((opt)->nf>=3?1:0)))
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IK(s,opt)   (ID(opt)+(s)-1)  /* GPS link IFCB, GPS sat numbers 1..NSATGPS */
#define IFB(opt)    (ID(opt)+NSATGPS) /* BDS receiver IF13 minus IF12 code bias */
#define IG(s,opt)   (IFB(opt)+1+(s)-NSATGPS-NSATGLO-1) /* GAL satellite-link CODE IF13 bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

static int no1_flag = 1;
static int isapplypppar = 0;
static int preweek = 0;
static double presecond = 0.;

typedef struct {
    int sat,frq,code;
    double elevation,y,model,residual,variance;
    char signal[24];
} ppp_eq_trace_t;

/* Return a stable, human-readable name for each PPP filter state. */
static void ppp_state_name(int index, const prcopt_t *opt, char *name,
                           size_t size)
{
    static const char *clock_name[]={
        "CLK_GPS","ISB_GLO","ISB_GAL","ISB_BDS","ISB_IRN","ISB_BDS3"
    };
    char id[16];
    int i,sat,frq;

    if (index<3) {
        snprintf(name,size,"POS_%c","XYZ"[index]);
        return;
    }
    if (opt->dynamics&&index<NP(opt)) {
        if (index<6) snprintf(name,size,"VEL_%c","XYZ"[index-3]);
        else         snprintf(name,size,"ACC_%c","XYZ"[index-6]);
        return;
    }
    if (IC(0,opt)<=index&&index<IC(0,opt)+NC(opt)) {
        i=index-IC(0,opt);
        snprintf(name,size,"%s",i<(int)(sizeof(clock_name)/sizeof(*clock_name))?
                 clock_name[i]:"ISB_UNKNOWN");
        return;
    }
    if (IT(opt)<=index&&index<IT(opt)+NT(opt)) {
        i=index-IT(opt);
        snprintf(name,size,"%s",i==0?"TROP_ZTD":(i==1?"TROP_GRAD_N":"TROP_GRAD_E"));
        return;
    }
    if (NI(opt)&&II(1,opt)<=index&&index<II(1,opt)+NI(opt)) {
        sat=index-II(1,opt)+1; satno2id(sat,id);
        snprintf(name,size,"ION_%s",id);
        return;
    }
    if (ND(opt)&&index>=ID(opt)&&index<ID(opt)+ND(opt)) {
        if (opt->if_model==1) {
            if (index==IFB(opt)) snprintf(name,size,"IFB_BDS_13_12");
            else if (index>IFB(opt)) snprintf(name,size,"IFB_E%02d",index-IFB(opt));
            else snprintf(name,size,"IFCB_G%02d",index-ID(opt)+1);
        }
        else snprintf(name,size,"DCB_L5");
        return;
    }
    if (index>=NR(opt)&&index<NX(opt)) {
        i=index-NR(opt); frq=i/MAXSAT; sat=i%MAXSAT+1; satno2id(sat,id);
        if (opt->ionoopt==IONOOPT_IFLC) {
            snprintf(name,size,"AMB_%s_IF1%d",id,frq+2);
        }
        else snprintf(name,size,"AMB_%s_F%d",id,frq+1);
        return;
    }
    snprintf(name,size,"STATE_%d",index);
}

/* Output the complete, accepted pre-fit PPP observation equations. */
static void trace_ppp_equations(gtime_t time, const ppp_eq_trace_t *eq,
                                int nv, const double *H, int nx)
{
    char line[8192],id[16],str[32];
    int i,j,n,first;

    if (gettracelevel()<4) return;
    time2str(time,str,3);
    trace(4,"========== PPP FILTER EQUATIONS time=%s nv=%d nx=%d ==========\n",
          str,nv,nx);
    for (i=0;i<nv;i++) {
        satno2id(eq[i].sat,id);
        n=snprintf(line,sizeof(line),
            "PPP_EQ row=%3d sat=%-3s el=%7.3f deg | %-9s %-5s "
            "v=%+11.4f var=%11.3e sigma=%9.4f y=%14.4f model=%14.4f | H={",
            i,id,eq[i].elevation,eq[i].signal,
            eq[i].code?"CODE":"PHASE",eq[i].residual,eq[i].variance,
            sqrt(eq[i].variance),eq[i].y,eq[i].model);
        first=1;
        for (j=0;j<nx&&n<(int)sizeof(line)-64;j++) {
            if (H[j+nx*i]==0.0) continue;
            n+=snprintf(line+n,sizeof(line)-n,"%s%d:%+.9g",first?"":",",
                        j,H[j+nx*i]);
            first=0;
        }
        snprintf(line+n,sizeof(line)-n,"}\n");
        trace(4,"%s",line);
    }
    trace(4,"========== END PPP FILTER EQUATIONS ==========\n");
}

/* Output the state vector and standard deviations before/after filtering. */
static void trace_ppp_filter_states(gtime_t time, const prcopt_t *opt,
                                    const double *xb, const double *sb,
                                    const double *xa, const double *Pa, int nx)
{
    char name[64],str[32];
    double sa;
    int i,active;

    if (gettracelevel()<4) return;
    time2str(time,str,3);
    trace(4,"PPP FILTER STATES: time=%s nx=%d\n",str,nx);
    trace(4,"INDEX NAME                              VALUE BEFORE -> VALUE AFTER"
            "       STD BEFORE -> STD AFTER  ACTIVE\n");
    for (i=0;i<nx;i++) {
        ppp_state_name(i,opt,name,sizeof(name));
        sa=Pa[i+i*nx]>0.0?sqrt(Pa[i+i*nx]):0.0;
        active=xa[i]!=0.0&&Pa[i+i*nx]!=0.0;
        trace(4,"%5d %-31s %+15.6f -> %+15.6f  %12.7f -> %12.7f  %s\n",
              i,name,xb[i],xa[i],sb[i],sa,active?"YES":"NO");
    }
    trace(4,"========== END PPP FILTER STATES ==========\n");
}

/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->nx]);
    return SQRT(rtk->P[i+i*rtk->nx]);
}
/* write solution status for PPP ---------------------------------------------*/
/* write solution status for PPP ---------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff)
{
    ssat_t *ssat;
    double pos[3], vel[3], acc[3], *x;
    int i, j;
    char id[32], *p = buff;

    if (!rtk->sol.stat) return 0;

    trace(3, "pppoutstat:\n");

    /* 生成时间字符串 */
    char sol_time_str[32];
    time2str(rtk->sol.time, sol_time_str, 3);  // 保留3位小数秒

    x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;

    /* receiver position */
    p += sprintf(p, "$POS,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                 sol_time_str, rtk->sol.stat, x[0], x[1], x[2],
                 STD(rtk, 0), STD(rtk, 1), STD(rtk, 2));

    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->x + 3, vel);
        ecef2enu(pos, rtk->x + 6, acc);
        p += sprintf(p, "$VELACC,%s,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
                     "%.4f,%.5f,%.5f,%.5f\n",
                     sol_time_str, rtk->sol.stat, vel[0], vel[1], vel[2],
                     acc[0], acc[1], acc[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    /* receiver clocks */
    i = IC(0, &rtk->opt);
    p += sprintf(p, "$CLK,%s,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
                 sol_time_str, rtk->sol.stat, 1, x[i] * 1E9 / CLIGHT,
                 x[i + 1] * 1E9 / CLIGHT, x[i + 2] * 1E9 / CLIGHT,
                 x[i + 3] * 1E9 / CLIGHT, STD(rtk, i) * 1E9 / CLIGHT,
                 STD(rtk, i + 1) * 1E9 / CLIGHT, STD(rtk, i + 2) * 1E9 / CLIGHT,
                 STD(rtk, i + 3) * 1E9 / CLIGHT);

    /* tropospheric parameters */
    if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TROP,%s,%d,%d,%.4f,%.4f\n",
                     sol_time_str, rtk->sol.stat, 1, x[i], STD(rtk, i));
    }
    if (rtk->opt.tropopt == TROPOPT_ESTG) {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TRPG,%s,%d,%d,%.5f,%.5f,%.5f,%.5f\n",
                     sol_time_str, rtk->sol.stat, 1, x[i + 1], x[i + 2],
                     STD(rtk, i + 1), STD(rtk, i + 2));
    }

    /* ionosphere parameters */
    if (rtk->opt.ionoopt == IONOOPT_EST) {
        for (i = 0; i < MAXSAT; i++) {
            ssat = rtk->ssat + i;
            if (!ssat->vs) continue;
            j = II(i + 1, &rtk->opt);
            if (rtk->x[j] == 0.0) continue;
            satno2id(i + 1, id);
            p += sprintf(p, "$ION,%s,%d,%s,%.1f,%.1f,%.4f,%.4f\n",
                         sol_time_str, rtk->sol.stat, id,
                         rtk->ssat[i].azel[0] * R2D, rtk->ssat[i].azel[1] * R2D,
                         x[j], STD(rtk, j));
        }
    }

    /* Compact proof of actual BDS participation in the new IF solution. */
    if (rtk->opt.bds_if==1) {
        int used=0,pilot=0,data=0,mixed=0;
        for (i=0;i<MAXSAT;i++) {
            ssat=rtk->ssat+i;
            if (satsys(i+1,NULL)!=SYS_CMP||!ssat->vsat[0]) continue;
            used++;
            if (ssat->obs_code[0]==CODE_L1P&&ssat->obs_code[1]==CODE_L5P) pilot++;
            else if (ssat->obs_code[0]==CODE_L1D&&ssat->obs_code[1]==CODE_L5D) data++;
            else mixed++;
        }
        p+=sprintf(p,"$BDSIF,%s,%d,%d,%d,%d,%d\n",sol_time_str,rtk->sol.stat,used,pilot,data,mixed);
    }
    if (rtk->opt.if_model) {
        int counts[6]={0},sys;
        for (i=0;i<MAXSAT;i++) {
            ssat=rtk->ssat+i;sys=satsys(i+1,NULL);satno2id(i+1,id);
            for (j=0;j<NF(&rtk->opt);j++) {
                if (!ssat->vsat[j]) continue;
                counts[(sys==SYS_GPS?0:sys==SYS_GAL?4:2)+j]++;
                p+=sprintf(p,"$IFPAIR,%s,%d,%s,1%d,%s,%s,%.5f,%.5f\n",
                    sol_time_str,rtk->sol.stat,id,j+2,code2obs(ssat->obs_code[0]),
                    code2obs(ssat->obs_code[j+1]),ssat->resc[j],ssat->resp[j]);
            }
            if (rtk->opt.if_model==1&&sys==SYS_GPS&&ssat->vsat[1])
                p+=sprintf(p,"$IFCB,%s,%s,%.6f,%.6f\n",sol_time_str,id,
                    x[IK(i+1,&rtk->opt)],STD(rtk,IK(i+1,&rtk->opt)));
            if (rtk->opt.if_model==1&&sys==SYS_GAL&&ssat->vsat[1])
                p+=sprintf(p,"$GALIFB,%s,%s,%.6f,%.6f\n",sol_time_str,id,x[IG(i+1,&rtk->opt)],STD(rtk,IG(i+1,&rtk->opt)));
        }
        p+=sprintf(p,"$IFSTAT,%s,%d,%d,%d,%d,%d\n",sol_time_str,rtk->sol.stat,
                   counts[0],counts[1],counts[2],counts[3]);
        if (rtk->opt.navsys&SYS_GAL) p+=sprintf(p,"$GALSTAT,%s,%d,%d,%d,BROADCAST_FNAV\n",
            sol_time_str,rtk->sol.stat,counts[4],counts[5]);
        if (rtk->opt.if_model==1) p+=sprintf(p,"$IFB,%s,BDS,%.6f,%.6f\n",sol_time_str,
            x[IFB(&rtk->opt)],STD(rtk,IFB(&rtk->opt)));
    }
    /* ambiguity parameters */
    int k;
    for (i = 0; i < MAXSAT; i++) for (j = 0; j < NF(&rtk->opt); j++) {
        k = IB(i + 1, j, &rtk->opt);
        if (rtk->x[k] == 0.0) continue;
        satno2id(i + 1, id);
        p += sprintf(p, "$AMB,%s,%d,%s,%d,%.4f,%.4f\n",
                     sol_time_str, rtk->sol.stat, id, j + 1, x[k], STD(rtk, k));
    }

    return (int)(p - buff);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;

    trace(3,"testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);

    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;

        if ((r=norm(rs+i*6,3))<=0.0) continue;

        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;

        /* sun-earth-satellite angle */
        cosa=dot3(rs+i*6,esun)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);

        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;

        trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
              obs[i].sat);

        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    *yaw=yaw_nominal(beta,mu);
    return 1;
}

/*****************************************************************************
* Name        : gravitationalDelayCorrection
* Description : Obtains the gravitational delay correction for the effect of
*               general relativity (red shift) to the GPS signal
* Parameters  :
* Name                           |Da|Unit|Description
* double  *receiverPosition       I  m    Position of the receiver
* double  *satellitePosition      I  m    Position of the satellite
* Returned value (double)         O  m    Gravitational delay correction
*****************************************************************************/
extern double gravitationalDelayCorrection (const int sys, const double *receiverPosition,
	                                        const double *satellitePosition)
{
	double	receiverModule;
	double	satelliteModule;
	double	distance;
	double  MU=MU_GPS;

	receiverModule=sqrt(receiverPosition[0]*receiverPosition[0]+receiverPosition[1]*receiverPosition[1]+
		receiverPosition[2]*receiverPosition[2]);
	satelliteModule=sqrt(satellitePosition[0]*satellitePosition[0]+satellitePosition[1]*satellitePosition[1]+
		satellitePosition[2]*satellitePosition[2]);
	distance=sqrt((satellitePosition[0]-receiverPosition[0])*(satellitePosition[0]-receiverPosition[0])+
		(satellitePosition[1]-receiverPosition[1])*(satellitePosition[1]-receiverPosition[1])+
		(satellitePosition[2]-receiverPosition[2])*(satellitePosition[2]-receiverPosition[2]));

	switch (sys) {
	case SYS_GPS:
		MU=MU_GPS;
		break;
	case SYS_GLO:
		MU=MU_GLO;
		break;
	case SYS_GAL:
		MU=MU_GAL;
		break;
	case SYS_CMP:
		MU=MU_CMP;
		break;
	default:
		MU=MU_GPS;
		break;
	}

	return 2.0*MU/(CLIGHT*CLIGHT)*log((satelliteModule+receiverModule+distance)/(satelliteModule+receiverModule-distance));
}

/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
                   const double *rs, double *exs, double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;

    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);

    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot3(esun,en));
    E=acos(dot3(es,ep));
    mu=PI/2.0+(dot3(es,esun)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;

    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
    int i;

    if (opt<=0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;

    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!normv3(r,ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek,eys,eks);
    cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot3(ek,exs)-eks[i];
        dr[i]=exr[i]-ek[i]*dot3(ek,exr)+ekr[i];
    }
    cosp=dot3(ds,dr)/norm(ds,3)/norm(dr,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot3(ek,drs)<0.0) ph=-ph;

    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, double snr_rover,
                     int f, const prcopt_t *opt, const obsd_t *obs)
{
    double a,b,e;
    double snr_max=opt->err[5];
    double fact=1.0;
    double sinel=sin(el),var;
    int frq,code;
    int curprn = 0;
    satsys(sat,&curprn);

    if(sys == SYS_CMP && !opt->if_model && (curprn==38 || curprn==40)){
        fact *= 3; 
    }

    if(sys == SYS_CMP && !opt->if_model && curprn==39){
        /* During the test, C39 showed the worst orbital accuracy, so the C39 weight was reduced by another two times. */
        fact *= 6; 
    }
		
    frq=f/2;code=f%2; /* 0=phase, 1=code */
    /* increase variance for pseudoranges */
    if (code) fact=opt->eratio[frq];
    if (fact<=0.0) fact=opt->eratio[0];
    /* adjust variances for constellation */
    switch (sys) {
        case SYS_GPS: fact*=EFACT_GPS;break;
        case SYS_GLO: fact*=EFACT_GLO;break;
        case SYS_GAL: fact*=EFACT_GAL;break;
        case SYS_SBS: fact*=EFACT_SBS;break;
        case SYS_QZS: fact*=EFACT_QZS;break;
        case SYS_CMP: fact*=EFACT_CMP;break;
        case SYS_IRN: fact*=EFACT_IRN;break;
        default:      fact*=EFACT_GPS;break;
    }
    if (sys==SYS_GPS||sys==SYS_QZS) {
        if (frq==2) fact*=EFACT_GPS_L5; /* GPS/QZS L5 error factor */
    }
    
    /* adjust variance for config parameters */
    a=fact*opt->err[1];  /* base term */
    b=fact*opt->err[2];  /* el term */

    // if(opt->ionoopt==IONOOPT_IFLC){
    //     if(opt->nf==1) fact*=SQR(0.5);
    //     else if(opt->nf==2) fact*=SQR(3.0);
    // }

    /* calculate variance */
    var=(a*a+b*b/sinel/sinel);
    // var=(a*a+b*b/sinel/sinel);
    if (opt->err[6]>0) {  /* add SNR term */
        e=fact*opt->err[6];
        var+=e*e*(pow(10,0.1*MAX(snr_max-snr_rover,0)));
    }
    if (opt->err[7]>0.0) {   /* add rcvr stdevs term */
        if (code) var+=SQR(opt->err[7]*0.01*(1<<(obs->Pstd[frq]+5))); /* 0.01*2^(n+5) */
        else var+=SQR(opt->err[7]*obs->Lstd[frq]*0.004*0.2); /* 0.004 cycles -> m) */
    }
    /* FIXME: the scaling factor is not 3 for other signals/constellations than GPS L1/L2 */
    var*=(opt->ionoopt==IONOOPT_IFLC&&!opt->if_model)?SQR(3.0):1.0;
    return var;
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) rtk->P[i+j*rtk->nx]=0.0;
    for (j=0;j<rtk->nx;j++) rtk->P[j+i*rtk->nx]=0.0;
    rtk->P[i+i*rtk->nx]=var;
}
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav)
{
    double freq1,freq2;

    freq1=sat2freq(obs->sat,obs->code[0],nav);
    freq2=sat2freq(obs->sat,obs->code[1],nav);
    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;
    return (obs->L[0]/freq1-obs->L[1]/freq2)*CLIGHT;
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav)
{
    double freq1,freq2;

    freq1=sat2freq(obs->sat,obs->code[0],nav);
    freq2=sat2freq(obs->sat,obs->code[1],nav);

    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0||
        obs->P[0]==0.0||obs->P[1]==0.0) return 0.0;
    trace(3,"mwmeas: %12.1f %12.1f %15.3f %15.3f %15.3f %15.3f %d %d\n",freq1,freq2,obs->L[0],obs->L[1],obs->P[0],obs->P[1],obs->code[0],obs->code[1]);
    return (obs->L[0]-obs->L[1])*CLIGHT/(freq1-freq2)-
           (freq1*obs->P[0]+freq2*obs->P[1])/(freq1+freq2);
}
/* Unified corrected observations: raw slots and up to two shared-anchor IFs. */
typedef struct {
    double raw_L[NFREQ],raw_P[NFREQ];
    double L[2],P[2],a[2][3],cov[2][2][2];
    int valid[NFREQ],pair[2][2],npair;
} meas_t;

/* Product conventions are independent of the number of IF combinations.
 * Keep historical SSR/file conventions for old configurations. Configured
 * GPS IF retains the documented baseline-code approximation; its IF13 link
 * IFCB is estimated in the filter, not looked up as an L2W correction. */
static int measurement_code_bias(const obsd_t *obs, const nav_t *nav,
                                 const prcopt_t *opt, int i, double *delta)
{
    int sys=satsys(obs->sat,NULL),ix=0,frq,bias_ix;
    double bias=0.0;
    *delta=0.0;
    if (opt->if_model) {
        if (sys==SYS_GAL&&!gal_code_bias(nav,obs->sat,obs->time,obs->code[i],&bias)) return 0;
        if (sys==SYS_CMP&&!bds_code_bias(obs->time,nav->B2bssr+obs->sat,obs->code[i],&bias))
            return 0;
        *delta=-bias;return 1;
    }
    if (opt->sateph==EPHOPT_SSRAPC||opt->sateph==EPHOPT_SSRCOM)
        *delta=nav->ssr[obs->sat-1].cbias[obs->code[i]-1];
    if (opt->sateph==EPHOPT_B2b) {
        if      (sys==SYS_GPS) ix=i==0?CODE_L1W:CODE_L2W;
        else if (sys==SYS_GLO) ix=i==0?CODE_L1P:CODE_L2P;
        else if (sys==SYS_GAL) ix=i==0?CODE_L1X:CODE_L7X;
        else if (sys==SYS_CMP) ix=obs->code[i];
        *delta-=nav->B2bssr[obs->sat].cbias[ix];
    }
    else {
        frq=sys==SYS_GAL&&(i==1||i==2)?3-i:i;
        if (frq<MAX_CODE_BIAS_FREQS&&(bias_ix=code2bias_ix(sys,obs->code[i]))>0)
            *delta+=nav->cbias[obs->sat-1][frq][bias_ix-1];
    }
    return 1;
}

/* One path for raw signal corrections and one path for IF construction.
 * rs==NULL is the ambiguity-initialization call: omit geometric antenna
 * terms, but retain code validity, frequencies and product bias handling. */
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *rs, const double *rr,
                      const double *e, double elapsed, double phw, meas_t *out)
{
    double freq[NFREQ]={0},rv[2][3]={{0}},delta,cbias;
    double dantr[NFREQ]={0},dants[NFREQ]={0},off[NFREQ][3]={{0}},rot[3];
    double dr,ds,ru[3],rz[3],eu[3],ez[3],nadir=0.0,cosa;
    uint8_t requested[3]={CODE_L1C,CODE_L2W,CODE_L5Q};
    int i,j,k,c,wi,sys=satsys(obs->sat,NULL),configured=opt->if_model!=0;
    int nraw=configured?(opt->if_model==1?3:2):MIN(opt->nf,NFREQ);
    prcopt_t noiseopt=*opt;
    memset(out,0,sizeof(*out));
    out->npair=opt->if_model==1?2:1;
    if (configured) {
        if (sys!=SYS_GPS&&sys!=SYS_CMP&&sys!=SYS_GAL) return;
        if (sys==SYS_GAL&&!if_parse_pairs(SYS_GAL,opt->gal_if_pairs,requested)) return;
        if (sys==SYS_CMP&&!bds_parse_freqs(opt->bds_freqs,requested)) return;
        if (rs) {
            if (!rr||!e) return;
            for (j=0;j<3;j++) {ru[j]=rr[j]-rs[j];rz[j]=-rs[j];}
            if (!normv3(ru,eu)||!normv3(rz,ez)) return;
            cosa=dot(eu,ez,3);nadir=acos(MAX(-1.0,MIN(1.0,cosa)));
        }
    }
    else {
        /* Legacy calibration/product policy, shared by IF and uncombined PPP. */
        if (sys==SYS_CMP&&opt->bds_if) for (i=0;i<2;i++) {
            if (bds_ant_index(nav->pcvs+obs->sat-1,obs->code[i],0)<0||
                (opt->sateph==EPHOPT_B2b&&!bds_code_bias(obs->time,
                 nav->B2bssr+obs->sat,obs->code[i],&cbias))) return;
        }
        if (rs) {
            if (!e) return;
            if (opt->posopt[0]&&!satantpcv(obs->sat,rs,rr,nav->pcvs+obs->sat-1,
                    sys==SYS_CMP&&opt->bds_if?obs->code:NULL,dants)) return;
            if (sys==SYS_CMP&&opt->bds_if) {
                if (!antmodel_bds(opt->pcvr,obs->code,opt->antdel[0],azel,
                                 opt->posopt[1],opt->bds_ant_fallback,dantr)) return;
            }
            else antmodel(obs->sat,opt->pcvr,opt->antdel[0],azel,opt->posopt[1],dantr);
            if (opt->sateph==EPHOPT_PREC||opt->sateph==EPHOPT_SSRCOM||opt->sateph==EPHOPT_B2b)
                satantoff1(obs->time,rs,obs->sat,nav,obs->code,off[0],off[1]);
        }
    }
    for (i=0;i<nraw;i++) {
        const char *s=code2obs(obs->code[i]);
        wi=0; /* legacy SNR mask selected its first group */
        if (configured) {
            if (s[0]!=code2obs(requested[i])[0]) continue;
            if (sys==SYS_CMP&&(s[0]=='2'||s[0]=='6'?s[1]!='I':s[1]!='P'&&s[1]!='D')) continue;
            if ((wi=signal_noise_index(sys,obs->code[i]))<0) continue;
            if (signal_ant_index(nav->pcvs+obs->sat-1,sys,obs->code[i],0)<0||
                ((opt->posopt[1]||opt->anttype[0][0]||opt->pcvr[0].type[0])&&
                 signal_ant_index(opt->pcvr,sys,obs->code[i],opt->bds_ant_fallback)<0)) continue;
        }
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);
        if (freq[i]<=0||!obs->P[i]||!obs->L[i]||
            testsnr(0,wi,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;
        if (!measurement_code_bias(obs,nav,opt,i,&delta)) continue;
        dr=dantr[i];ds=dants[i];
        if (rs&&configured) {
            if (!signal_antmodel(opt->pcvr,sys,obs->code[i],opt->antdel[0],azel,0,
                                 opt->posopt[1],opt->bds_ant_fallback,&dr)||
                !signal_antmodel(nav->pcvs+obs->sat-1,sys,obs->code[i],NULL,NULL,nadir,
                                 opt->posopt[0],0,&ds)||
                !signal_satantoff(obs->time,rs,obs->sat,nav,obs->code[i],off[i])) continue;
        }
        /* Shared carrier conversion, phase windup and geometric correction. */
        out->raw_L[i]=obs->L[i]*CLIGHT/freq[i]-phw*CLIGHT/freq[i];
        out->raw_P[i]=obs->P[i]+delta;
        if (rs) {
            rot[0]=off[i][0]+off[i][1]*OMGE*elapsed;
            rot[1]=off[i][1]-off[i][0]*OMGE*elapsed;rot[2]=off[i][2];
            out->raw_L[i]-=dr+ds+dot(rot,e,3);
            out->raw_P[i]-=dr+ds+dot(rot,e,3);
        }
        out->valid[i]=1;
        if (configured) {
            noiseopt.eratio[i]=opt->eratio[wi];
            for (c=0;c<2;c++) rv[c][i]=varerr(obs->sat,sys,azel[1],
                obs->SNR[i]*SNR_UNIT,i*2+c,&noiseopt,obs);
        }
    }
    for (k=0;k<out->npair;k++) {
        j=k+1;
        /* Only old configurations retain their historical missing-L2 fallback. */
        if (!configured&&!(sys==SYS_CMP&&opt->bds_if)&&!out->raw_L[1]) j=2;
        out->pair[k][0]=0;out->pair[k][1]=j;
        if (configured&&(!out->valid[0]||!out->valid[j])) continue;
        if (!if_coefficients(freq[0],freq[j],out->a[k],out->a[k]+j)) continue;
        if (out->raw_L[0]&&out->raw_L[j])
            out->L[k]=out->a[k][0]*out->raw_L[0]+out->a[k][j]*out->raw_L[j];
        if (out->raw_P[0]&&out->raw_P[j])
            out->P[k]=out->a[k][0]*out->raw_P[0]+out->a[k][j]*out->raw_P[j];
    }
    if (configured) for (c=0;c<2;c++) for (i=0;i<out->npair;i++) for (j=0;j<out->npair;j++)
        out->cov[c][i][j]=if_covariance(out->a[i],out->a[j],rv[c]);
}

/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i,j,nf=rtk->opt.nf;

    trace(3,"detslp_ll: n=%d\n",n);

    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<rtk->opt.nf;j++) {
        ssat_t *ssat=rtk->ssat+obs[i].sat-1;
        if ((rtk->opt.if_model||(rtk->opt.bds_if&&satsys(obs[i].sat,NULL)==SYS_CMP))&&obs[i].code[j]&&obs[i].L[j]) {
            if (ssat->obs_code[j]&&ssat->obs_code[j]!=obs[i].code[j]) {
                ssat->slip[j]|=1;
                ssat->gf[0]=ssat->mw[0]=0.0;
                if (rtk->opt.if_model) {
                    ssat->gf[1]=ssat->mw[1]=0.0;
                    if (satsys(obs[i].sat,NULL)==SYS_GPS&&rtk->opt.if_model==1)
                        initx(rtk,0.0,0.0,IK(obs[i].sat,&rtk->opt));
                    if (satsys(obs[i].sat,NULL)==SYS_GAL&&rtk->opt.if_model==1)
                        initx(rtk,0.0,0.0,IG(obs[i].sat,&rtk->opt));
                }
            }
            if (ssat->obs_code[j]!=obs[i].code[j]) {
                trace(1,"selected signal: sat=%d slot=%d code=%d (%s)\n",
                      obs[i].sat,j,obs[i].code[j],code2obs(obs[i].code[j]));
            }
            ssat->obs_code[j]=obs[i].code[j];
        }
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;

        trace(3,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);

        rtk->ssat[obs[i].sat-1].slip[j<nf?j:nf]=1;
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double g0,g1;
    int i,j;

    trace(4,"detslp_gf: n=%d\n",n);

    for (i=0;i<n&&i<MAXOBS;i++) {

        if ((g1=gfmeas(obs+i,nav))==0.0) continue;

        g0=rtk->ssat[obs[i].sat-1].gf[0];
        rtk->ssat[obs[i].sat-1].gf[0]=g1;

        trace(4,"detslip_gf: sat=%2d gf0=%8.3f gf1=%8.3f\n",obs[i].sat,g0,g1);

        if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
            trace(3,"detslip_gf: slip detected sat=%2d gf=%8.3f->%8.3f\n",
                  obs[i].sat,g0,g1);

            for (j=0;j<rtk->opt.nf;j++) rtk->ssat[obs[i].sat-1].slip[j]|=1;
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double w0,w1;
    int i,j;

    trace(4,"detslp_mw: n=%d\n",n);

    for (i=0;i<n&&i<MAXOBS;i++) {
        if ((w1=mwmeas(obs+i,nav))==0.0) continue;

        w0=rtk->ssat[obs[i].sat-1].mw[0];
        rtk->ssat[obs[i].sat-1].mw[0]=w1;

        trace(4,"detslip_mw: sat=%2d mw0=%8.3f mw1=%8.3f\n",obs[i].sat,w0,w1);

        if (w0!=0.0&&fabs(w1-w0)>THRES_MW_JUMP) {
            trace(3,"detslip_mw: slip detected sat=%2d mw=%8.3f->%8.3f\n",
                  obs[i].sat,w0,w1);

			for (j=0;j<rtk->opt.nf;j++) rtk->ssat[obs[i].sat-1].slip[j]|=1;
        }
	}
}
/* Pair-local GF/MW tests do not falsely identify which raw carrier slipped. */
static void detslp_if(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    int i,k,j;
    for (i=0;i<n&&i<MAXOBS;i++) for (k=0;k<NF(&rtk->opt);k++) {
        obsd_t pair=obs[i];
        ssat_t *s=rtk->ssat+obs[i].sat-1;
        double gf,mw;
        j=k+1;
        pair.L[1]=obs[i].L[j];pair.P[1]=obs[i].P[j];pair.code[1]=obs[i].code[j];
        gf=gfmeas(&pair,nav);mw=mwmeas(&pair,nav);
        if (!(s->slip[0]||s->slip[j])) {
            if ((gf&&s->gf[k]&&fabs(gf-s->gf[k])>rtk->opt.thresslip)||
                (mw&&s->mw[k]&&fabs(mw-s->mw[k])>THRES_MW_JUMP)) {
                s->if_slip[k]=1;
                trace(3,"IF1213 pair slip sat=%d IF1%d\n",obs[i].sat,k+2);
            }
        }
        s->gf[k]=gf;s->mw[k]=mw;
    }
}

/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
    double *F,*P,*FP,*x,*xp,pos[3],Q[9]={0},Qv[9],var=0.0;
    int i,j,*ix,nx;

    trace(3,"udpos_ppp:\n");

    /* fixed mode */
    if (rtk->opt.mode==PMODE_PPP_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
    }
    /* static ppp mode */
    if (rtk->opt.mode==PMODE_PPP_STATIC) {
        for (i=0;i<3;i++) {
            rtk->P[i*(1+rtk->nx)]+=SQR(rtk->opt.prn[5])*fabs(rtk->tt);
        }
        return;
    }
    /* kinematic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i=0;i<3;i++) {
            initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        }
        return;
    }
    /* check variance of estimated position */
    for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx];
    var/=3.0;

    if (var>VAR_POS) {
        /* reset position with large variance */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    /* generate valid state index */
    ix=imat(rtk->nx,1);
    for (i=nx=0;i<rtk->nx;i++) {
        if  (i<9||(rtk->x[i]!=0.0&&rtk->P[i+i*rtk->nx]>0.0)) ix[nx++]=i;
    }
    /* state transition of position/velocity/acceleration */
    F=eye(nx); P=mat(nx,nx); FP=mat(nx,nx); x=mat(nx,1); xp=mat(nx,1);

    for (i=0;i<6;i++) {
        F[i+(i+3)*nx]=rtk->tt;
    }
    /* include accel terms if filter is converged */
    if (var<rtk->opt.thresar[1]) {
        for (i=0;i<3;i++) {
            F[i+(i+6)*nx]=SQR(rtk->tt)/2.0;
        }
    }
    else trace(3,"pos var too high for accel term: %.4f,%.4f\n", var,rtk->opt.thresar[1]);
    for (i=0;i<nx;i++) {
        x[i]=rtk->x[ix[i]];
        for (j=0;j<nx;j++) {
            P[i+j*nx]=rtk->P[ix[i]+ix[j]*rtk->nx];
        }
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",nx,1,nx,F,x,xp);
    matmul("NN",nx,nx,nx,F,P,FP);
    matmul("NT",nx,nx,nx,FP,F,P);

    for (i=0;i<nx;i++) {
        rtk->x[ix[i]]=xp[i];
        for (j=0;j<nx;j++) {
            rtk->P[ix[i]+ix[j]*rtk->nx]=P[i+j*nx];
        }
    }
    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3])*fabs(rtk->tt);
    Q[8]=SQR(rtk->opt.prn[4])*fabs(rtk->tt);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        rtk->P[i+6+(j+6)*rtk->nx]+=Qv[i+j*3];
    }
    free(ix); free(F); free(P); free(FP); free(x); free(xp);
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk)
{
    double dtr;
	int i, numofsys;

	trace(3,"udclk_ppp:\n");

#ifdef BDS2BDS3
	numofsys = NSYS + 1;
#else
	numofsys = NSYS;
#endif

    /* initialize every epoch for clock (white noise) */
	for (i=0;i<numofsys;i++) {
        if (rtk->opt.sateph==EPHOPT_PREC) {
            /* time of prec ephemeris is based gpst */
            /* neglect receiver inter-system bias  */
			if (i == 0) dtr = rtk->sol.dtr[0];
            else dtr = rtk->sol.dtr[0] + rtk->sol.dtr[i];
        }
        else {
			dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i];
        }
		initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt));
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i=IT(&rtk->opt),j;

    trace(3,"udtrop_ppp:\n");

    if (rtk->x[i]==0.0) {
        ecef2pos(rtk->sol.rr,pos);
        ztd=sbstropcorr(rtk->sol.time,pos,azel,&var);
        initx(rtk,ztd,var,i);

        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) initx(rtk,1E-6,VAR_GRA,j);
        }
    }
    else {
        rtk->P[i+i*rtk->nx]+=SQR(rtk->opt.prn[2])*fabs(rtk->tt);

        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) {
                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt);
            }
        }
    }
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double freq1,freq2,ion,sinel,pos[3],*azel;
    char *p;
    int i,j,f2,gap_resion=GAP_RESION,sat;

    trace(3,"udiono_ppp:\n");

    if ((p=strstr(rtk->opt.pppopt,"-GAP_RESION="))) {
        sscanf(p,"-GAP_RESION=%d",&gap_resion);
    }
    /* reset ionosphere delay estimate if outage too long */
    for (i=0;i<MAXSAT;i++) {
        j=II(i+1,&rtk->opt);
        if (rtk->x[j]!=0.0&&(int)rtk->ssat[i].outc[0]>gap_resion) {
            rtk->x[j]=0.0;
        }
    }
    for (i=0;i<n;i++) {
        sat=obs[i].sat;
        j=II(sat,&rtk->opt);
        if (rtk->x[j]==0.0) {
            /* initialize ionosphere delay estimates if zero */
            f2=seliflc(rtk->opt.nf,satsys(sat,NULL));
            freq1=sat2freq(sat,obs[i].code[0],nav);
            freq2=sat2freq(sat,obs[i].code[f2],nav);
            if (obs[i].P[0]==0.0||obs[i].P[f2]==0.0||freq1==0.0||freq2==0.0) {
                continue;
            }
            /* use pseudorange difference adjusted by freq for initial estimate */
            ion=(obs[i].P[0]-obs[i].P[f2])/(SQR(FREQL1/freq1)-SQR(FREQL1/freq2));
            ecef2pos(rtk->sol.rr,pos);
            azel=rtk->ssat[sat-1].azel;
            /* adjust delay estimate by path length */
            ion/=ionmapf(pos,azel);
            initx(rtk,ion,VAR_IONO,j);
            trace(4,"ion init: sat=%d ion=%.4f\n",sat,ion);
        }
        else {
            sinel=sin(MAX(rtk->ssat[sat-1].azel[1],5.0*D2R));
            /* update variance of delay state */
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]/sinel)*fabs(rtk->tt);
        }
    }
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t *rtk)
{
    int i=ID(&rtk->opt);

    if (rtk->opt.if_model) {
        if (rtk->opt.if_model==1) {
            double q=rtk->opt.ifcb_prn>0?rtk->opt.ifcb_prn:0.001;
            int s;
            for (s=1;s<=NSATGPS;s++) {
                i=IK(s,&rtk->opt);
                if (rtk->x[i]!=0.0) rtk->P[i+i*rtk->nx]+=q*q*fabs(rtk->tt);
            }
            i=IFB(&rtk->opt);
            if (!rtk->x[i]) initx(rtk,1E-6,VAR_DCB,i);
            if (rtk->opt.navsys&SYS_GAL) for (s=1;s<=NSATGAL;s++) {
                i=IG(satno(SYS_GAL,s),&rtk->opt);
                if (rtk->x[i]) rtk->P[i+i*rtk->nx]+=q*q*fabs(rtk->tt);
            }
        }
        return;
    }

    trace(3,"uddcb_ppp:\n");

    if (rtk->x[i]==0.0) {
        initx(rtk,1E-6,VAR_DCB,i);
    }
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double Lc,Pc,bias[MAXOBS],offset=0.0,pos[3]={0};
    double freq1,freq2,ion;
    int i,j,k,f,sat,slip[MAXOBS]={0},clk_jump=0;
    meas_t ic;double *L=ic.raw_L,*P=ic.raw_P;
    char satid[4];
    //trace(3, "udbias  : n=%d\n", n);

    /* handle day-boundary clock jump */
    if (rtk->opt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) {
        rtk->ssat[i].slip[j]=0;
        if (j<2) rtk->ssat[i].if_slip[j]=0;
    }
    /* detect cycle slip by LLI */
    detslp_ll(rtk,obs,n);

    if (rtk->opt.if_model) detslp_if(rtk,obs,n,nav);
    else {
        detslp_gf(rtk,obs,n,nav);
        detslp_mw(rtk,obs,n,nav);
    }

    ecef2pos(rtk->sol.rr,pos);

    for (f=0;f<NF(&rtk->opt);f++) {
        offset=0.0;
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            satno2id(i, satid);
            trace(3, "udbias  : %s %s rtk->ssat[%d].outc[%d] = %d\n", time_str(obs[0].time, 2), satid, i, f, rtk->ssat[i].outc[f] + 1);
            if (++rtk->ssat[i].outc[f]>60||
                rtk->opt.modear==ARMODE_INST||clk_jump) {
                initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&rtk->opt);
            corr_meas(obs+i,nav,rtk->ssat[sat-1].azel,&rtk->opt,NULL,NULL,NULL,0,0,&ic);
            Lc=ic.L[rtk->opt.if_model?f:0];Pc=ic.P[rtk->opt.if_model?f:0];

            bias[i]=0.0;

            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                bias[i]=Lc!=0.0&&Pc!=0.0?Lc-Pc:0.0;
                slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[rtk->opt.if_model?f+1:1];
                if (rtk->opt.if_model) {
                    slip[i]|=rtk->ssat[sat-1].if_slip[f];
                    if (bias[i]&&f==1) {
                        if (satsys(sat,NULL)==SYS_GPS) {
                            int ik=IK(sat,&rtk->opt);
                            if (!rtk->x[ik]) {
                                double v=ic.P[0]!=0.0?ic.P[1]-ic.P[0]:0.0;
                                initx(rtk,v==0.0?1E-6:v,VAR_DCB,ik);
                            }
                            /* P13 = geometry+K, L13 = geometry-K+N13. */
                            bias[i]+=2.0*rtk->x[ik];
                        }
                        else if (satsys(sat,NULL)==SYS_GAL) {
                            int ig=IG(sat,&rtk->opt);
                            if (!rtk->x[ig]) {
                                double v=ic.P[0]?ic.P[1]-ic.P[0]:0.0;
                                initx(rtk,v==0?1E-6:v,VAR_DCB,ig);
                            }
                            bias[i]+=rtk->x[ig];
                        }
                        else bias[i]+=rtk->x[IFB(&rtk->opt)];
                    }
                }

            }
            else if (L[f]!=0.0&&P[f]!=0.0) {
                freq1=sat2freq(sat,obs[i].code[0],nav);
                freq2=sat2freq(sat,obs[i].code[f],nav);
                slip[i]=rtk->ssat[sat-1].slip[f];
                if (f==0||obs[i].P[0]==0.0||obs[i].P[f]==0.0||freq1==0.0||freq2==0.0)
                    ion=0;
                else
                    ion=(obs[i].P[0]-obs[i].P[f])/(1.0-SQR(freq1/freq2));
                bias[i]=L[f]-P[f]+2.0*ion*SQR(freq1/freq2);
            }
            if (rtk->x[j]==0.0||slip[i]||bias[i]==0.0) continue;

            offset+=bias[i]-rtk->x[j];
            k++;
        }
        /* correct phase-code jump to ensure phase-code coherence */
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
            for (i=0;i<MAXSAT;i++) {
                j=IB(i+1,f,&rtk->opt);
                if (rtk->x[j]!=0.0) rtk->x[j]+=offset/k;
            }
            trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                  time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
        }
        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&rtk->opt);

            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);

            if (bias[i]==0.0||(rtk->x[j]!=0.0&&!slip[i])) continue;

            /* reinitialize phase-bias if detecting cycle slip */
            initx(rtk,bias[i],VAR_BIAS,IB(sat,f,&rtk->opt));

            /* reset fix flags */
            for (k=0;k<MAXSAT;k++) rtk->ambc[sat-1].flags[k]=0;

            satno2id(sat, satid);
            trace(3,"udbias_ppp: sat=%2d %s bias=%.3f\n",sat,satid,bias[i]);
        }
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    int curweek = 0;
	double cursecond = 0.;

	cursecond = time2gpst(obs[0].time, &curweek);

    if(no1_flag == 1){
        preweek = curweek;
		presecond = cursecond;
        no1_flag = 0;
    }

	if ((cursecond - presecond) + (curweek - preweek)*604800. >= 8*3600.)
	{
		printf("Reset beginning time...\n");

		// for (int ix = 0; ix < rtk->nx; ix++)
		// 	initx(rtk, rtk->x[ix], VAR_BIAS, ix);

		preweek = curweek;
		presecond = cursecond;
	}


    time_t times = (time_t)1742515228;
    BOOL flag = times == obs[0].time.time ? TRUE : FALSE;

    trace(3,"udstate_ppp: n=%d\n",n);

    /* temporal update of position */
    udpos_ppp(rtk);
    /* temporal update of clock */
    udclk_ppp(rtk);
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        udtrop_ppp(rtk);
    }
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt==IONOOPT_EST) {
        udiono_ppp(rtk,obs,n,nav);
    }
    /* temporal update of L5-receiver-dcb parameters */
    if (rtk->opt.nf>=3) {
        uddcb_ppp(rtk);
    }
    /* temporal update of phase-bias */
    udbias_ppp(rtk,obs,n,nav);
}
/* satellite antenna phase center variation ----------------------------------*/
static int satantpcv(int sat,const double *rs, const double *rr, const pcv_t *pcv,
                      const uint8_t *code, double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return 0;

    cosa=dot3(eu,ez);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    if (code) return antmodel_s_bds(pcv,code,nadir,dant);
	antmodel_s(sat,pcv,nadir,dant);
    return 1;
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var)
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;

    /* zenith hydrostatic delay */
	zhd=tropmodel(time,pos,zazel,0.0);

    /* mapping function */
	m_h=tropmapf(time,pos,azel,&m_w);

    if (azel[1]>0.0) {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*(x[0]-zhd);
        dtdx[2]=grad_e*(x[0]-zhd);
    }
    dtdx[0]=m_w;
    *var=SQR(0.01);
    return m_h*zhd+m_w*(x[0]-zhd);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp, double *var)
{
    double trp[3]={0};

    if (opt->tropopt==TROPOPT_SAAS) {
		*dtrp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt==TROPOPT_SBAS) {
        *dtrp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        matcpy(trp,x+IT(opt),opt->tropopt==TROPOPT_EST?1:3,1);
		*dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);
        return 1;
    }
    return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var)
{
    if (opt->ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,dion,var);
    }
    if (opt->ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,dion,var);
    }
    if (opt->ionoopt==IONOOPT_BRDC) {
        *dion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*dion*ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt==IONOOPT_EST) {
        /* Estimated delay is a vertical delay, apply the mapping function. */
        *dion=x[II(sat,opt)]*ionmapf(pos,azel);
        *var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        *dion=*var=0.0;
        return 1;
    }
    return 0;
}
/* phase and code residuals --------------------------------------------------*/
/* PRNs are identifiers, not permanent BDS generations. The ANTEX record has
 * already been selected for the observation epoch. Modern B1C/B2a tracking
 * supplies a fallback if there is no descriptive satellite antenna record. */
static int ppp_bds3(const obsd_t *obs, const nav_t *nav)
{
    const pcv_t *pcv=nav->pcvs+obs->sat-1;
    int i,prn=0;
    if (pcv->sat==obs->sat) {
        if (!strncmp(pcv->type,"BEIDOU-3",8)) return 1;
        if (!strncmp(pcv->type,"BEIDOU-2",8)) return 0;
    }
    for (i=0;i<NFREQ;i++) {
        uint8_t c=obs->code[i];
        if (c==CODE_L1D||c==CODE_L1P||c==CODE_L1X||
            c==CODE_L5D||c==CODE_L5P||c==CODE_L5X) return 1;
    }
    satsys(obs->sat,&prn);
    return prn>16; /* compatibility for legacy observations without metadata */
}

static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *var_rs, const int *svh,
                   const double *dr, int *exc, const nav_t *nav,
                   const double *x, rtk_t *rtk, double *v, double *H, double *R,
                   double *azel)
{
    prcopt_t *opt=&rtk->opt;
    double y,r,cdtr,bias,rr[3],pos[3],e[3],dtdx[3];
    double var[MAXOBS*2*NFREQ],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0,dcb,freq;
    double crossvar[MAXOBS*2*NFREQ];
    double commonvar[MAXOBS*2*NFREQ];
    int rowsat[MAXOBS*2*NFREQ],rowcode[MAXOBS*2*NFREQ],rowpair[MAXOBS*2*NFREQ];
    meas_t ic;
	double ve[MAXOBS * 2 * NFREQ] = { 0 }, vr[MAXOBS * 2 * NFREQ] = { 0 }, vmax = 0;
    ppp_eq_trace_t eq[MAXOBS*2*NFREQ];
    char str[32];char id[8];
    int ne=0,obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ],maxobs,maxfrq,rej;
	int i,j,k,sat,sys,nv=0,nx=rtk->nx,stat=1,frq,code;
	double res = 0.0;
	double gravitationalDelayModel = 0.;
	int curprn = 0;
	int ambpos = 0;
    double debug_epoch[6];
    double stop_epoch[6] = {2025,3,7,5,54,00};

	time2str(obs[0].time,str,2);

    time2epoch(obs[0].time,debug_epoch);

    if (fabs(debug_epoch[0] - stop_epoch[0]) < 1e-6 &&  // 年
    fabs(debug_epoch[1] - stop_epoch[1]) < 1e-6 &&  // 月
    fabs(debug_epoch[2] - stop_epoch[2]) < 1e-6 &&  // 日
    fabs(debug_epoch[3] - stop_epoch[3]) < 1e-6 &&  // 时
    fabs(debug_epoch[4] - stop_epoch[4]) < 1e-6 &&  // 分
    fabs(debug_epoch[5] - stop_epoch[5]) < 1e-3) {  // 秒
    int a = 0;
}


    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;

    for (i=0;i<3;i++) rr[i]=x[i]+dr[i];
    ecef2pos(rr,pos);

    time_t times = (time_t)1742515228;
    BOOL flag = times == obs[0].time.time ? TRUE : FALSE;

    for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat; satno2id(sat,id);
        int nv_sat=nv;

        if ((r=geodist(rs+i*6,rr,e))<=0.0||
            satazel(pos,e,azel+i*2)<opt->elmin) {
			exc[i]=1;
            continue;
        }
		/* Configured IF uses its own physical measurements and B2b DCB.
           A signal unavailable to SPP (e.g. CNAV1 without TGD) can still
           be valid PPP data. Geometry, health, DCB and residual checks follow. */
        if (!(sys=satsys(sat,/*NULL*/&curprn))||
            (!opt->if_model&&!rtk->ssat[sat-1].vs)||
            satexclude(sat,var_rs[i],svh[i],opt)||exc[i]) {
            exc[i]=1;
            continue;
		}

        /* tropospheric and ionospheric model */
		if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,&vart)||
            !model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
            continue;
        }
        /* phase windup model */
        if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {
            continue;
		}

		//Gravitational delay correction */
		gravitationalDelayModel=gravitationalDelayCorrection(sys,rr,rs+i*6);

        /* corrected phase and code measurements */
        corr_meas(obs+i,nav,azel+i*2,opt,rs+i*6,rr,e,r/CLIGHT,rtk->ssat[sat-1].phw,&ic);

        /* stack phase and code residuals {L1,P1,L2,P2,...} */
        for (j=0;j<2*NF(opt);j++) {
            double C=0.0;

            dcb=bias=0.0;
            code=j%2; /* 0=phase, 1=code */
            frq=j/2;

            if (opt->ionoopt==IONOOPT_IFLC) {
                y=code?ic.P[frq]:ic.L[frq];
                if (y==0.0) continue;
            }
            else {
                if ((y=code==0?ic.raw_L[frq]:ic.raw_P[frq])==0.0) continue;

                if ((freq=sat2freq(sat,obs[i].code[frq],nav))==0.0) continue;
                /* The iono paths have already applied a slant factor. */
                C=SQR(FREQL1/freq)*(code==0?-1.0:1.0);
            }
            if (H) {
                for (k=0;k<nx;k++) H[k+nx*nv]=0.0;
                for (k=0;k<3;k++) H[k+nx*nv]=-e[k];
            }

            /* receiver clock */
			switch (sys) {
                case SYS_GLO: k=1; break;
                case SYS_GAL: k=2; break;
                case SYS_CMP: k=3; break;
                case SYS_IRN: k=4; break;
                default:      k=0; break;
			}

#ifdef	BDS2BDS3
			if(sys==SYS_CMP && (opt->if_model?ppp_bds3(obs+i,nav):curprn>16))
				k = NSYS;
#endif
			cdtr=x[IC(k,opt)];
			if (H) {
				H[IC(k,opt)+nx*nv]=1.0;

                if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                    for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                        H[IT(opt)+k+nx*nv]=dtdx[k];
                    }
                }
            }
            if (opt->ionoopt==IONOOPT_EST) {
                if (rtk->x[II(sat,opt)]==0.0) continue;
                /* The vertical iono delay is estimated, but the residual is
                 * in the direction of the slant, so apply the slat factor
                 * mapping function. */
                if (H) H[II(sat,opt)+nx*nv]=C*ionmapf(pos,azel+i*2);
            }
            if (opt->if_model==1&&frq==1) {
                if (sys==SYS_GPS) {
                    int ik=IK(sat,opt);
                    double sign=code?1.0:-1.0;
                    if (!x[ik]) continue;
                    dcb=sign*x[ik];
                    if (H) H[ik+nx*nv]=sign;
                }
                else if (code) {
                    int idx=sys==SYS_GAL?IG(sat,opt):IFB(opt);
                    dcb=x[idx];
                    if (H) H[idx+nx*nv]=1.0;
                }
            }
            else if (!opt->if_model&&frq==2&&code==1) { /* legacy L5-receiver-dcb */
                dcb+=rtk->x[ID(opt)];
                if (H) H[ID(opt)+nx*nv]=1.0;
            }
            if (code==0) { /* phase bias */
				if ((bias=x[IB(sat,frq,opt)])==0.0) continue;

				if (H) H[IB(sat,frq,opt)+nx*nv]=1.0;
            }
            /* residual */
			res=y-(r+cdtr-CLIGHT*dts[i*2]+dtrp+C*dion+dcb+bias-gravitationalDelayModel);
            if (v) v[nv]=res;

            if (code==0) rtk->ssat[sat-1].resc[frq]=res;  /* carrier phase */
            else         rtk->ssat[sat-1].resp[frq]=res;  /* pseudorange */

            /* variance */
            var[nv]=opt->if_model?ic.cov[code][frq][frq]:varerr(sat,sys,azel[1+i*2],
                    SNR_UNIT*rtk->ssat[sat-1].snr_rover[frq],
					j,opt,obs+i);
            commonvar[nv]=opt->if_model&&sys==SYS_GAL?var_rs[i]:0.0;
            var[nv]+=commonvar[nv]; /* also used by the postfit outlier test */
			if (sys==SYS_GLO&&code==1) var[nv]+=VAR_GLO_IFB;
			// if(sys == SYS_CMP  && code==1 && (curprn<=5||curprn>=59/*||curprn==31||curprn==38||curprn==39||curprn==40*/))
			// 	var[nv] *= 100; 
            // if(sys == SYS_CMP  && (curprn==31||curprn==38||curprn==39||curprn==40)){
            //     var[nv] *= 2; 
            // }
				
			
            if(!post)
                trace(3,"ppp_Res(%d)%s sat=%2d %s%d res=%9.4f y=%.4f r=%.4f cdtr=%.4f cdts=%.4f dtrp=%.4f dion=%.4f bias=%.4f dcb=%.4f gravity=%.4f el=%4.1f\n",post,str,sat,
				  code?"P":"L",frq+1,res,y,r, cdtr, CLIGHT * dts[i * 2], dtrp,dion, bias,dcb, gravitationalDelayModel,azel[1+i*2]*R2D);
            
            /* reject satellite by pre-fit residuals */
            if (!post&&opt->maxinno[code]>0.0&&fabs(res)>opt->maxinno[code]) {
				trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
                      post,str,sat,code?"P":"L",frq+1,res,azel[1+i*2]*R2D);
                exc[i]=1; rtk->ssat[sat-1].rejc[frq]++;
                if (opt->if_model) {
                    /* Discard earlier rows of this satellite too. */
                    nv=nv_sat;
                    rtk->ssat[sat-1].vsat[0]=rtk->ssat[sat-1].vsat[1]=0;
                    break;
                }
                continue;
            }
            /* record large post-fit residuals */
			if (post&&fabs(res)>sqrt(var[nv])*THRES_REJECT) {
				obsi[ne] = i; frqi[ne] = j; ve[ne] = res; vr[ne] = sqrt(var[nv]); ne++;
            }
            if (code==0) rtk->ssat[sat-1].vsat[frq]=1;
            rowsat[nv]=sat;rowcode[nv]=code;rowpair[nv]=frq;
            crossvar[nv]=opt->if_model?ic.cov[code][0][1]:0.0;
            if (!post) {
                eq[nv].sat=sat; eq[nv].frq=frq; eq[nv].code=code;
                eq[nv].elevation=azel[1+i*2]*R2D;
                eq[nv].y=y; eq[nv].model=y-res; eq[nv].residual=res;
                eq[nv].variance=var[nv];
                if (opt->ionoopt==IONOOPT_IFLC) {
                    snprintf(eq[nv].signal,sizeof(eq[nv].signal),"%s/%s",
                             code2obs(obs[i].code[0]),code2obs(obs[i].code[1]));
                    if (opt->if_model) snprintf(eq[nv].signal,sizeof(eq[nv].signal),"%s/%s",
                        code2obs(obs[i].code[0]),code2obs(obs[i].code[frq+1]));
                }
                else snprintf(eq[nv].signal,sizeof(eq[nv].signal),"%s",
                              code2obs(obs[i].code[frq]));
            }
            nv++;
        }
    }
    /* reject satellite with large and max post-fit residual */
    if (post&&ne>0) {
        vmax=ve[0]/vr[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
        for (j=1;j<ne;j++) {
			if (fabs(vmax) >= fabs(ve[j] / vr[j])) continue;
			vmax = ve[j] / vr[j]; maxobs = obsi[j]; maxfrq = frqi[j]; rej = j;
        }
        sat=obs[maxobs].sat;
        trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
            post,str,sat,maxfrq%2?"P":"L",maxfrq/2+1,vmax,azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; rtk->ssat[sat-1].rejc[maxfrq/2]++; stat=0;
        ve[rej]=0;
    }
    if (R) {
        for (j=0;j<nv;j++) for (i=0;i<nv;i++) R[i+j*nv]=0.0;
        for (i=0;i<nv;i++) R[i+i*nv]=var[i];
        if (opt->if_model) for (i=0;i<nv;i++) for (j=i+1;j<nv;j++) {
            if (rowsat[i]==rowsat[j]) R[i+j*nv]=R[j+i*nv]=commonvar[i];
            if (rowsat[i]==rowsat[j]&&rowcode[i]==rowcode[j]&&rowpair[i]!=rowpair[j])
                R[i+j*nv]=R[j+i*nv]=commonvar[i]+crossvar[i];
        }
    }
    if (!post&&H&&nv>0) trace_ppp_equations(obs[0].time,eq,nv,H,nx);
    return post?stat:nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
    const prcopt_t *opt=&rtk->opt;
    int i,j;
    char satid[4];
    /* test # of valid satellites */
    rtk->sol.ns=0;
    for (i=0;i<n&&i<MAXOBS;i++) {
        int used=0;
        for (j=0;j<NF(opt);j++) {
            if (!rtk->ssat[obs[i].sat-1].vsat[j]) continue;
            rtk->ssat[obs[i].sat-1].lock[j]++;
            rtk->ssat[obs[i].sat-1].outc[j]=0;
            satno2id(obs[i].sat - 1, satid);
            trace(3, "update_stat: %s %s rtk->ssat[%d].outc[%d] = %d\n", time_str(obs[0].time, 2), satid, obs[i].sat - 1, j, rtk->ssat[obs[i].sat - 1].outc[j]);
            used=1;
        }
        if (used) rtk->sol.ns++;
    }
    rtk->sol.stat=rtk->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;

    if (rtk->sol.stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];

        if (rtk->opt.dynamics) { /* velocity and covariance */
            for (i=3;i<6;i++) {
                rtk->sol.rr[i]=rtk->x[i];
                rtk->sol.qv[i-3]=(float)rtk->P[i+i*rtk->nx];
            }
            rtk->sol.qv[3]=(float)rtk->P[4+3*rtk->nx];
            rtk->sol.qv[4]=(float)rtk->P[5+4*rtk->nx];
            rtk->sol.qv[5]=(float)rtk->P[5+3*rtk->nx];
        }
    }
		rtk->sol.dtr[0]=rtk->x[IC(0,opt)]; /* GPS */
        rtk->sol.dtr[1]=rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)]; /* GLO-GPS */
        rtk->sol.dtr[2]=rtk->x[IC(2,opt)]-rtk->x[IC(0,opt)]; /* GAL-GPS */
        rtk->sol.dtr[3]=rtk->x[IC(3,opt)]-rtk->x[IC(0,opt)]; /* BDS-GPS */

#ifdef	BDS2BDS3
		rtk->sol.dtr[NSYS] = rtk->x[IC(NSYS, opt)] - rtk->x[IC(0, opt)]; /* BDS-GPS */
#endif

    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        rtk->ssat[obs[i].sat-1].snr_rover[j]=obs[i].SNR[j];
        rtk->ssat[obs[i].sat-1].snr_base[j] =0;
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (rtk->ssat[i].slip[j]&3) rtk->ssat[i].slipc[j]++;
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
    }
}
/* test hold ambiguity -------------------------------------------------------*/
static int test_hold_amb(rtk_t *rtk)
{
    int i,j,stat=0;

    /* no fix-and-hold mode */
	if (rtk->opt.modear!=ARMODE_FIXHOLD) return 0;

    /* reset # of continuous fixed if new ambiguity introduced */
    for (i=0;i<MAXSAT;i++) {
        if (rtk->ssat[i].fix[0]!=2&&rtk->ssat[i].fix[1]!=2) continue;
        for (j=0;j<MAXSAT;j++) {
            if (rtk->ssat[j].fix[0]!=2&&rtk->ssat[j].fix[1]!=2) continue;
            if (!rtk->ambc[j].flags[i]||!rtk->ambc[i].flags[j]) stat=1;
            rtk->ambc[j].flags[i]=rtk->ambc[i].flags[j]=1;
        }
    }
    if (stat) {
        rtk->nfix=0;
        return 0;
    }
    /* test # of continuous fixed */
    return ++rtk->nfix>=rtk->opt.minfix;
}


/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const prcopt_t *opt=&rtk->opt;
	double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp,*xb,*sb,dr[3]={0},std[3];
	char str[32];
	int i,j,nv,info,svh[MAXOBS],brdc[MAXOBS],exc[MAXOBS]={0},stat=SOLQ_SINGLE;
	double dpos[3] = {0x00};

	time2str(obs[0].time,str,2);
	trace(2,"pppos   : time=%s nx=%d n=%d\n",str,rtk->nx,n);

	rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);

	for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;
	for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
		rtk->ssat[obs[i].sat-1].snr_rover[j]=obs[i].SNR[j];
		rtk->ssat[obs[i].sat-1].snr_base[j] =0;
	}

    time_t times = (time_t)1742515228;
    BOOL flag = times == obs[0].time.time ? TRUE : FALSE;
	/* temporal update of ekf states */
	udstate_ppp(rtk,obs,n,nav);

	/* satellite positions and clocks */
    satposs_if(obs[0].time,obs,n,nav,&rtk->opt,rs,dts,var,svh);
    /* GAL may not participate in the SPP initializer (external F/NAV only).
     * Its availability for PPP is determined by its own valid orbit/clock. */
    if (rtk->opt.if_model&&(rtk->opt.navsys&SYS_GAL)) for (i=0;i<n;i++)
        if (satsys(obs[i].sat,NULL)==SYS_GAL)
            rtk->ssat[obs[i].sat-1].vs=norm(rs+6*i,3)>0&&svh[i]==0;

	/* exclude measurements of eclipsing satellite (block IIA) */
	if (rtk->opt.posopt[3]) {
		testeclipse(obs,n,nav,rs);
	}
	/* earth tides correction */
	if (opt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time),rtk->x,opt->tidecorr==1?1:3,&nav->erp,
				 opt->odisp[0],dr);
	}
	nv=n*rtk->opt.nf*2+MAXSAT+3;
	xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);
	xb=mat(rtk->nx,1); sb=mat(rtk->nx,1);
	v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv);

	for (i=0;i<MAX_ITER;i++) {

		matcpy(xp,rtk->x,rtk->nx,1);
		matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

		/* prefit residuals */
		if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel))) {
			trace(2,"%s ppp (%d) no valid obs data\n",str,i+1);
			break;
		}

		/* measurement update of ekf states */
		for (j=0;j<rtk->nx;j++) {
			xb[j]=xp[j];
			sb[j]=Pp[j+j*rtk->nx]>0.0?sqrt(Pp[j+j*rtk->nx]):0.0;
		}
		if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) {
			trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
			break;
		}
		trace_ppp_filter_states(obs[0].time,opt,xb,sb,xp,Pp,rtk->nx);

		/* postfit residuals */
		if (ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) {
			matcpy(rtk->x,xp,rtk->nx,1);
			matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
			stat=SOLQ_PPP;
			break;
		}
	}

	if (i>=MAX_ITER) {
		trace(2,"%s ppp (%d) iteration overflows\n",str,i);
	}

	if (stat==SOLQ_PPP) {

		if (ppp_res(9,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) {

			matcpy(rtk->xa,xp,rtk->nx,1);
			matcpy(rtk->Pa,Pp,rtk->nx,rtk->nx);

			/*ambiguity resolution in ppp*/
			if (opt->modear!=ARMODE_OFF)
			{
				// if (isapplypppar && pppamb(rtk,obs,n,nav,azel,exc)) {
                if (isapplypppar) {
                    if(pppamb(rtk,obs,n,nav,azel,exc)){
                        stat=SOLQ_FIX;
                    }
                    else stat=SOLQ_PPP;
				}
				else
				{
					stat=SOLQ_PPP;
					rtk->sol.ratio = 0.;
				}
			}
			else
			{
				stat = SOLQ_PPP;
			}

		}
		else {
			rtk->nfix=0;
        }

        /* update solution status */
        update_stat(rtk,obs,n,stat);

		if (stat==SOLQ_FIX&&test_hold_amb(rtk)) {
			matcpy(rtk->x,xp,rtk->nx,1);
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
            // matcpy(rtk->x,rtk->xa,rtk->nx,1);
            // matcpy(rtk->P,rtk->Pa,rtk->nx,rtk->nx);
            trace(2,"%s hold ambiguity\n",str);
            rtk->nfix=0;
        }
    }
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(xb); free(sb); free(v); free(H); free(R);
}
