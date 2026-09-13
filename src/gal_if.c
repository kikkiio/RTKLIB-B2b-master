/* Experimental GAL broadcast assistance, not precise orbit/clock products.
 * Galileo OS SIS ICD 2.2 sec 5.1.3--5: F/NAV clock datum is E1/E5a.
 * Broadcast IF APC is converted to COM before physical signal corrections. */
#include "rtklib.h"
#define SQR(x) ((x)*(x))
#define MAX(a,b) ((a)>(b)?(a):(b))

extern const eph_t *gal_select_eph(const nav_t *nav,int sat,gtime_t time,int fnav)
{
    const eph_t *best=NULL;double age,minage=MAXDTOE_GAL+1;int i;
    for (i=0;i<nav->ngal_eph;i++) {
        const eph_t *p=nav->gal_eph+i;
        if (p->sat!=sat||!(p->code&(1<<(fnav?8:9)))||p->svh||p->A<=0) continue;
        age=timediff(time,p->toe);
        if (age<0||age>MAXDTOE_GAL||timediff(time,p->ttr)<0) continue;
        if (age<minage) {best=p;minage=age;}
    }
    return best;
}

extern int gal_load_nav(nav_t *nav,const char *file,char *msg)
{
    nav_t *data=calloc(1,sizeof(*data));obs_t obs={0};sta_t sta={0};gtime_t z={0};
    int i,n=0,fn=0,in=0;eph_t *records;char path[MAXSTRPATH];
    if (!data) {strcpy(msg,"GAL navigation allocation failed");return 0;}
    if (strlen(file)>=sizeof(path)) {strcpy(msg,"GAL navigation path too long");free(data);return 0;}
    strcpy(path,file);
#ifdef _WIN32
    /* expath() recognises native separators when rebuilding absolute paths. */
    for (i=0;path[i];i++) if (path[i]=='/') path[i]='\\';
#endif
    if (!path[0]||readrnxt(path,1,z,z,0,"-SYS=E",&obs,data,&sta)<=0) {
        strcpy(msg,"cannot read GAL RINEX navigation file");free(obs.data);freenav(data,0xFF);free(data);return 0;
    }
    for (i=0;i<data->n;i++) if (satsys(data->eph[i].sat,NULL)==SYS_GAL) {
        n++;fn+=(data->eph[i].code&(1<<8))!=0;in+=(data->eph[i].code&(1<<9))!=0;
    }
    records=n?malloc(n*sizeof(*records)):NULL;
    if (!records||!fn||!in) {
        free(records);strcpy(msg,"GAL requires F/NAV clocks and I/NAV BGD records");
        free(obs.data);freenav(data,0xFF);free(data);return 0;
    }
    n=0;for (i=0;i<data->n;i++) if (satsys(data->eph[i].sat,NULL)==SYS_GAL) records[n++]=data->eph[i];
    free(nav->gal_eph);nav->gal_eph=records;nav->ngal_eph=n;
    trace(1,"GAL NAV: records=%d F/NAV=%d I/NAV=%d file=%s\n",n,fn,in,file);
    free(obs.data);freenav(data,0xFF);free(data);return 1;
}

/* Biases relative to F/NAV: zero IF12 BGD, IF13 c*(BGD_a-BGD_b).
 * No B2b GAL correction is fabricated. Phase constants remain ambiguities. */
extern int gal_code_bias(const nav_t *nav,int sat,gtime_t time,uint8_t code,double *bias)
{
    const eph_t *a=gal_select_eph(nav,sat,time,1),*b;
    char band=code2obs(code)[0];
    if (!a||!isfinite(a->tgd[0])) return 0;
    if (band=='1') *bias=CLIGHT*a->tgd[0];
    else if (band=='5') *bias=CLIGHT*SQR(FREQL1/FREQL5)*a->tgd[0];
    else if (band=='7') {
        b=gal_select_eph(nav,sat,time,0);
        if (!b||!isfinite(b->tgd[1])) return 0;
        *bias=CLIGHT*(a->tgd[0]+(SQR(FREQL1/FREQE5b)-1)*b->tgd[1]);
    }
    else return 0;
    return 1;
}

extern void satposs_if(gtime_t teph,const obsd_t *obs,int n,const nav_t *nav,
                      const prcopt_t *opt,double *rs,double *dts,double *var,int *svh)
{
    int i,j;double pr,dt,r2[3],dt2,v2,off1[3],off2[3],a,b;
    const eph_t *eph;gtime_t tx;
    satposs(teph,obs,n,nav,opt->sateph,rs,dts,var,svh);
    if (!opt->if_model||opt->gal_ephemeris!=1) return;
    for (i=0;i<n;i++) {
        if (satsys(obs[i].sat,NULL)!=SYS_GAL) continue;
        memset(rs+6*i,0,6*sizeof(double));memset(dts+2*i,0,2*sizeof(double));svh[i]=-1;
        for (j=0,pr=0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0) break;
        if (!pr) continue;
        tx=timeadd(obs[i].time,-pr/CLIGHT);
        if (!(eph=gal_select_eph(nav,obs[i].sat,tx,1))) continue;
        dt=eph2clk(tx,eph);tx=timeadd(tx,-dt);
        eph2pos(tx,eph,rs+6*i,dts+2*i,var+i);
        eph2pos(timeadd(tx,0.001),eph,r2,&dt2,&v2);
        for (j=0;j<3;j++) rs[6*i+j+3]=(r2[j]-rs[6*i+j])/0.001;
        dts[2*i+1]=(dt2-dts[2*i])/0.001;
        if (!signal_satantoff(tx,rs+6*i,obs[i].sat,nav,CODE_L1C,off1)||
            !signal_satantoff(tx,rs+6*i,obs[i].sat,nav,CODE_L5Q,off2)) {
            memset(rs+6*i,0,6*sizeof(double));continue;
        }
        if_coefficients(FREQL1,FREQL5,&a,&b);
        for (j=0;j<3;j++) rs[6*i+j]-=a*off1[j]+b*off2[j];
        var[i]=MAX(var[i],SQR(opt->gal_brdc_sigma>0?opt->gal_brdc_sigma:3.0));
        svh[i]=eph->svh;
    }
}
