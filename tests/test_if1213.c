/* White-box regression tests exercise the actual PPP measurement/state code. */
#include "rtklib.h"
#include "../src/ppp.c"
#include "fixtures/measurement_before_unification.h"
PPPGlobal_t PPP_Glo={0};
static int checks=0,failures=0;
#define CHECK(x) do {checks++;if (!(x)) {fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x);failures++;}} while (0)
#define NEAR(a,b,t) CHECK(fabs((a)-(b))<(t))

static prcopt_t options(void)
{
    prcopt_t o=prcopt_default;
    o.mode=PMODE_PPP_KINEMA;o.modear=ARMODE_OFF;o.nf=3;
    o.navsys=SYS_GPS|SYS_CMP;o.sateph=EPHOPT_B2b;o.ionoopt=IONOOPT_IFLC;
    o.if_model=1;o.bds_if=2;
    return o;
}

static void coefficients(void)
{
    const double f[3][3]={{FREQL1,FREQL2,FREQL5},{FREQL5,FREQL1,FREQ1_CMP},
                          {FREQ1_CMP,FREQ3_CMP,FREQL5}};
    double a[2][3]={{0}},v[3]={1,4,9},raw[3],combined;
    int s,k;
    for (s=0;s<3;s++) for (k=0;k<2;k++) {
        memset(a[k],0,sizeof(a[k]));
        CHECK(if_coefficients(f[s][0],f[s][k+1],a[k],a[k]+k+1));
        NEAR(a[k][0]+a[k][k+1],1.0,1E-12);
        raw[0]=10;raw[k+1]=10*SQR(f[s][0]/f[s][k+1]);
        combined=a[k][0]*raw[0]+a[k][k+1]*raw[k+1];NEAR(combined,0,1E-12);
        CHECK((s!=1)==(a[k][0]>0));
    }
    NEAR(if_covariance(a[0],a[1],v),a[0][0]*a[1][0],1E-12);
    CHECK(if_covariance(a[0],a[0],v)*if_covariance(a[1],a[1],v)>
          SQR(if_covariance(a[0],a[1],v)));
    CHECK(!if_coefficients(0,FREQL1,raw,raw+1));
    CHECK(!if_coefficients(FREQL1,FREQL1,raw,raw+1));
}

static void selection(void)
{
    prcopt_t o=options();char msg[128],opt[64];
    obsd_t a={0},b,c;
    int i;
    a.sat=satno(SYS_CMP,23);
    a.code[0]=CODE_L2I;a.code[1]=CODE_L6I;a.code[3]=CODE_L1P;a.code[4]=CODE_L5P;
    for (i=0;i<NFREQ+NEXOBS;i++) {a.P[i]=2E7+i;a.L[i]=1E8+i;a.LLI[i]=i%2;a.SNR[i]=40000+i;}
    bds_select_obs(&a,&b,1,2);bds_select_obs(&b,&c,1,2);
    CHECK(b.code[0]==CODE_L5P&&b.code[1]==CODE_L1P&&b.code[2]==CODE_L2I);
    CHECK(b.P[0]==a.P[4]&&b.L[2]==a.L[0]&&b.LLI[1]==a.LLI[3]);
    CHECK(!memcmp(&b,&c,sizeof(b)));
    bds_decode_options(opt,sizeof(opt),"",2);
    CHECK(code2obsidx(SYS_CMP,CODE_L5P,opt)==0);
    CHECK(code2obsidx(SYS_CMP,CODE_L1P,opt)==1);
    CHECK(code2obsidx(SYS_CMP,CODE_L2I,opt)==2);
    CHECK(code2idx(SYS_CMP,CODE_L5P)==3);
    CHECK(bds_options_valid(&o,msg));CHECK(NF(&o)==2);CHECK(ND(&o)==NSATGPS+1);
    CHECK(IB(1,1,&o)-IB(1,0,&o)==MAXSAT);
    o.modear=ARMODE_CONT;CHECK(!bds_options_valid(&o,msg));o.modear=ARMODE_OFF;
    o.if_model=2;CHECK(NF(&o)==1&&ND(&o)==0);CHECK(bds_options_valid(&o,msg));
    o.if_model=0;CHECK(!bds_options_valid(&o,msg));
    o.bds_if=1;o.nf=2;CHECK(bds_options_valid(&o,msg));CHECK(NF(&o)==1);
}

static void antenna(void)
{
    pcv_t p={0};double d,az[2]={0,PI/2},del[3]={0};
    p.var=calloc(NANTFREQ,sizeof(*p.var));strcpy(p.type,"TEST");
    p.valid[ANT_GPS_L5]=1;p.off[ANT_GPS_L5][2]=0.12;
    p.valid[1]=1;p.off[1][2]=9.0;
    CHECK(signal_ant_index(&p,SYS_GPS,CODE_L5Q,0)==ANT_GPS_L5);
    CHECK(signal_antmodel(&p,SYS_GPS,CODE_L5Q,del,az,0,0,0,&d));NEAR(d,-0.12,1E-12);
    p.valid[ANT_GPS_L5]=0;CHECK(signal_ant_index(&p,SYS_GPS,CODE_L5Q,0)<0);
    p.valid[3*NFREQ+1]=1;CHECK(signal_ant_index(&p,SYS_GPS,CODE_L5Q,1)==3*NFREQ+1);
    p.sat=satno(SYS_GPS,3);CHECK(signal_ant_index(&p,SYS_GPS,CODE_L5Q,1)<0);
    free(p.var);
}

static void equations(void)
{
    nav_t *nav=calloc(1,sizeof(*nav));rtk_t rtk={0};
    prcopt_t o=options();obsd_t obs={0};
    double ep[6]={2025,3,21,0,0,30},rs[6]={26378137,12000000,15000000},rr[3]={6378137,0,0};
    double dts[2]={0},vr[1]={0},dr[3]={0},az[2],v[4],R[16],*H,dist,e[3];
    double f,ion,a,b,K=2.5;int svh[1]={0},exc[1]={0},n,i;
    o.tropopt=TROPOPT_SAAS;o.posopt[2]=0;o.maxinno[0]=o.maxinno[1]=0;
    o.pcvr[0].var=calloc(NANTFREQ,sizeof(*o.pcvr[0].var));strcpy(o.pcvr[0].type,"TEST");
    o.pcvr[0].valid[0]=o.pcvr[0].valid[1]=o.pcvr[0].valid[ANT_GPS_L5]=1;
    obs.time=epoch2time(ep);obs.sat=satno(SYS_GPS,4);
    nav->pcvs[obs.sat-1]=o.pcvr[0];nav->pcvs[obs.sat-1].sat=obs.sat;
    obs.code[0]=CODE_L1C;obs.code[1]=CODE_L2W;obs.code[2]=CODE_L5Q;
    dist=geodist(rs,rr,e);if_coefficients(FREQL1,FREQL5,&a,&b);
    for (i=0;i<3;i++) {
        f=sat2freq(obs.sat,obs.code[i],nav);ion=10*SQR(FREQL1/f);
        obs.P[i]=dist+ion+(i==2?K/b:0);obs.L[i]=(dist-ion+5-(i==2?K/b:0))*f/CLIGHT;
        obs.SNR[i]=45000;
    }
    rtkinit(&rtk,&o);rtk.sol.time=obs.time;memcpy(rtk.sol.rr,rr,sizeof(rr));
    for (i=0;i<3;i++) initx(&rtk,rr[i],1,i);
    initx(&rtk,K,1,IK(obs.sat,&o));
    initx(&rtk,5,1,IB(obs.sat,0,&o));initx(&rtk,5,1,IB(obs.sat,1,&o));
    rtk.ssat[obs.sat-1].vs=0; /* independent PPP must not require SPP participation */
    H=calloc(rtk.nx*4,sizeof(double));
    n=ppp_res(0,&obs,1,rs,dts,vr,svh,dr,exc,nav,rtk.x,&rtk,v,H,R,az);
    CHECK(n==4);
    NEAR(H[IK(obs.sat,&o)+2*rtk.nx],-1,1E-12);
    NEAR(H[IK(obs.sat,&o)+3*rtk.nx],1,1E-12);
    NEAR(H[IB(obs.sat,0,&o)],1,1E-12);NEAR(H[IB(obs.sat,1,&o)+2*rtk.nx],1,1E-12);
    CHECK(R[0+2*4]>0&&R[1+3*4]>0);NEAR(R[0+2*4],R[2+0*4],1E-14);
    CHECK(R[0]*R[2+2*4]>SQR(R[0+2*4]));NEAR(R[0+1*4],0,1E-14);
    NEAR(v[0],v[2],1E-7);NEAR(v[1],v[3],1E-7);
    /* Capacity-only stress: duplicate synthetic satellites are NOT filtered.
     * This deliberately exceeds the old MAXOBS*2 local variance buffer. */
    {
        obsd_t many[MAXOBS];double rss[MAXOBS*6],dt[MAXOBS*2]={0},vvv[MAXOBS]={0},azs[MAXOBS*2];
        int ex[MAXOBS]={0},health[MAXOBS]={0},cap=MAXOBS*4;
        double *vv=calloc(cap+1,sizeof(double)),*hh=calloc(rtk.nx*cap,sizeof(double)),*rrr=calloc(cap*cap,sizeof(double));
        for (i=0;i<MAXOBS;i++) {many[i]=obs;memcpy(rss+6*i,rs,sizeof(rs));}
        vv[cap]=12345;
        n=ppp_res(0,many,MAXOBS,rss,dt,vvv,health,dr,ex,nav,rtk.x,&rtk,vv,hh,rrr,azs);
        CHECK(n==cap);NEAR(vv[cap],12345,1E-12);
        CHECK(isfinite(rrr[cap*cap-1])&&rrr[cap*cap-1]>0);
        free(vv);free(hh);free(rrr);
    }
    obs.L[2]=0;exc[0]=0;
    n=ppp_res(0,&obs,1,rs,dts,vr,svh,dr,exc,nav,rtk.x,&rtk,v,H,R,az);CHECK(n==2);
    /* A code outlier must also remove the phase row already stacked. */
    rtk.opt.maxinno[1]=1;obs.P[0]+=10000;exc[0]=0;
    n=ppp_res(0,&obs,1,rs,dts,vr,svh,dr,exc,nav,rtk.x,&rtk,v,H,R,az);CHECK(n==0);
    CHECK(!rtk.ssat[obs.sat-1].vsat[0]);
    free(H);rtkfree(&rtk);free(o.pcvr[0].var);free(nav);
}

static void measurement_and_slip(const char *freqs)
{
    nav_t *nav=calloc(1,sizeof(*nav));rtk_t rtk={0};
    prcopt_t o=options();obsd_t obs={0};meas_t x,y;
    double az[2]={0,PI/3},f,range=24000000.0,ion;
    double ep[6]={2025,3,21,0,0,30};int i;
    obs.time=epoch2time(ep);obs.sat=satno(SYS_CMP,23);
    strcpy(o.bds_freqs,freqs);CHECK(bds_parse_freqs(freqs,obs.code));
    o.pcvr[0].valid[ANT_B2A]=o.pcvr[0].valid[ANT_B1C]=o.pcvr[0].valid[2*NFREQ]=1;
    o.pcvr[0].valid[2*NFREQ+1]=1;
    nav->pcvs[obs.sat-1]=o.pcvr[0];nav->pcvs[obs.sat-1].sat=obs.sat;
    nav->B2bssr[obs.sat].t0[1]=obs.time;
    nav->B2bssr[obs.sat].t0[0]=nav->B2bssr[obs.sat].t0[2]=obs.time;
    for (i=0;i<3;i++) {
        f=sat2freq(obs.sat,obs.code[i],nav);ion=10*SQR(FREQL5/f);
        double cb=obs.code[i]==CODE_L6I?0:i;
        obs.P[i]=range+ion+cb;obs.L[i]=(range-ion+5)*f/CLIGHT;obs.SNR[i]=45000;
        nav->B2bssr[obs.sat].cbias_valid[obs.code[i]]=obs.code[i]!=CODE_L6I;
        nav->B2bssr[obs.sat].cbias[obs.code[i]]=(float)cb;
    }
    corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&x);
    for (i=0;i<2;i++) {NEAR(x.P[i],range,1E-7);NEAR(x.L[i],range+5,1E-7);CHECK(x.cov[0][i][i]>0);}
    CHECK(x.cov[0][0][1]>0);NEAR(x.cov[0][1][0],x.cov[0][0][1],1E-15);
    {
        pcv_t saved=o.pcvr[0];int pcvopt=o.posopt[1];
        memset(o.pcvr,0,sizeof(o.pcvr));o.posopt[1]=0;
        corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&y);
        NEAR(y.P[0],x.P[0],1E-8);NEAR(y.P[1],x.P[1],1E-8);
        o.posopt[1]=1;corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&y);
        CHECK(!y.P[0]&&!y.P[1]); /* enabling unknown PCV must not fabricate it */
        o.posopt[1]=0;strcpy(o.anttype[0],"MISSING CALIBRATION");
        corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&y);CHECK(!y.P[0]&&!y.P[1]);
        o.anttype[0][0]=0;o.pcvr[0]=saved;o.posopt[1]=pcvopt;
    }
    nav->B2bssr[obs.sat].cbias_valid[obs.code[2]]=0;
    corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&y);CHECK(y.P[0]!=0&&y.P[1]==0);
    nav->B2bssr[obs.sat].cbias_valid[obs.code[2]]=1;
    {
        double saved=obs.L[1];obs.L[1]=0;
        corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&y);CHECK(y.P[0]==0&&y.P[1]!=0);
        obs.L[1]=saved;
    }
    rtkinit(&rtk,&o);rtk.sol.time=obs.time;rtk.sol.rr[0]=6378137;
    memcpy(rtk.ssat[obs.sat-1].azel,az,sizeof(az));
    udbias_ppp(&rtk,&obs,1,nav);
    NEAR(rtk.x[IB(obs.sat,0,&o)],5,1E-7);NEAR(rtk.x[IB(obs.sat,1,&o)],5,1E-7);
    rtk.x[IB(obs.sat,0,&o)]=11;rtk.x[IB(obs.sat,1,&o)]=12;
    obs.LLI[2]=1;udbias_ppp(&rtk,&obs,1,nav);
    NEAR(rtk.x[IB(obs.sat,0,&o)],11,1E-9);NEAR(rtk.x[IB(obs.sat,1,&o)],5,1E-7);
    obs.LLI[2]=0;obs.LLI[0]=1;udbias_ppp(&rtk,&obs,1,nav);
    NEAR(rtk.x[IB(obs.sat,0,&o)],5,1E-7);NEAR(rtk.x[IB(obs.sat,1,&o)],5,1E-7);
    rtkfree(&rtk);free(nav);
}

static void configurable_selection(void)
{
    const char *names[4]={"B1I","B3I","B1C","B2a"};
    uint8_t codes[4]={CODE_L2I,CODE_L6I,CODE_L1P,CODE_L5P},parsed[3];
    prcopt_t o=options();obsd_t raw={0},out,again,spp;nav_t *nav=calloc(1,sizeof(*nav));
    eph_t eph={0};char text[64],decoder[128],msg[128];int i,j,k,q;
    double ep[6]={2025,3,21,0,0,30},bias;
    raw.time=epoch2time(ep);raw.sat=satno(SYS_CMP,23);
    for (i=0;i<4;i++) {
        raw.code[i]=codes[i];raw.P[i]=24000000+i;raw.L[i]=120000000+i;
        raw.LLI[i]=i;raw.SNR[i]=44000+i;raw.Pstd[i]=i;raw.D[i]=(float)i;
    }
    for (i=0;i<4;i++) for (j=0;j<4;j++) for (k=0;k<4;k++) {
        if (i==j||i==k||j==k) continue;
        snprintf(text,sizeof(text),"%s,%s,%s",names[i],names[j],names[k]);strcpy(o.bds_freqs,text);
        CHECK(bds_parse_freqs(text,parsed));CHECK(bds_options_valid(&o,msg));
        CHECK(parsed[0]==codes[i]&&parsed[1]==codes[j]&&parsed[2]==codes[k]);
        bds_decode_options_ex(decoder,sizeof(decoder),"-EPHALL -CBDSFREQ=512",&o);
        for (q=0;q<3;q++) CHECK(code2obsidx(SYS_CMP,parsed[q],decoder)==q);
        bds_select_obs_ex(&raw,&out,1,&o);bds_select_obs_ex(&out,&again,1,&o);
        CHECK(!memcmp(&out,&again,sizeof(out)));
        CHECK(out.code[0]==codes[i]&&out.code[1]==codes[j]&&out.code[2]==codes[k]);
        CHECK(out.P[0]==raw.P[i]&&out.L[1]==raw.L[j]&&out.LLI[2]==raw.LLI[k]);
        for (q=0;q<4;q++) {
            int z,found=0;for (z=0;z<NFREQ+NEXOBS;z++) found+=out.code[z]==codes[q];
            CHECK(found==1); /* preserve non-PPP B1C/B3I exactly once for SPP */
        }
        eph.sat=raw.sat;eph.toe=raw.time;eph.tgd_valid=1<<2;nav->eph=&eph;nav->n=1;
        bds_select_spp_obs(&out,&spp,1,nav);CHECK(spp.code[0]==CODE_L1P);
        eph.tgd_valid=0;bds_select_spp_obs(&out,&spp,1,nav);CHECK(spp.code[0]==CODE_L6I);
    }
    CHECK(bds_parse_freqs(" b1i , b3i , b2a ",parsed));
    CHECK(!bds_parse_freqs("B1I,B1I,B2a",parsed));
    CHECK(bds_parse_freqs("B1I,B3I",parsed));CHECK(!parsed[2]);
    CHECK(!bds_parse_freqs("B1I,B3I,B2a,",parsed));
    CHECK(!bds_parse_freqs("B1I,B3I,B2a,B1C",parsed));
    CHECK(!bds_parse_freqs("B1I,B3I,B2b",parsed));
    strcpy(o.bds_freqs,"B1I,B3I,B2a");o.bds_if=1;CHECK(!bds_options_valid(&o,msg));
    CHECK(code2obsidx(SYS_CMP,CODE_L2I,"-CBDSIF=2 -CBDSFREQ=22")<0);
    CHECK(signal_noise_index(SYS_CMP,CODE_L5P)==2);
    CHECK(signal_noise_index(SYS_CMP,CODE_L2I)==0);
    CHECK(signal_noise_index(SYS_CMP,CODE_L6I)==1);
    CHECK(signal_noise_index(SYS_CMP,CODE_L1P)==0);
    nav->B2bssr[raw.sat].t0[0]=nav->B2bssr[raw.sat].t0[2]=raw.time;
    CHECK(bds_code_bias(raw.time,nav->B2bssr+raw.sat,CODE_L6I,&bias));NEAR(bias,0,1E-12);
    CHECK(!bds_code_bias(timeadd(raw.time,43),nav->B2bssr+raw.sat,CODE_L6I,&bias));
    nav->B2bssr[raw.sat].iodssr[2]=1;
    CHECK(!bds_code_bias(raw.time,nav->B2bssr+raw.sat,CODE_L6I,&bias));
    {
        double az[2]={0,PI/3},newvar,oldvar;
        o=options();newvar=varerr(raw.sat,SYS_CMP,az[1],45,0,&o,&raw);
        o.if_model=0;oldvar=varerr(raw.sat,SYS_CMP,az[1],45,0,&o,&raw);
        NEAR(oldvar,newvar*9,1E-12);
    }
    free(nav);
}

static void pair_config_tests(void)
{
    const char *names[4]={"B1I","B3I","B1C","B2a"};
    prcopt_t o;uint8_t codes[3];char msg[128],decoder[128],pairs[64];
    obsd_t raw={0},selected,again;int i,j,k,q;
    FILE *fp;solopt_t sol;filopt_t fil;rtk_t rtk;
    raw.sat=satno(SYS_CMP,23);
    raw.code[0]=CODE_L2I;raw.code[1]=CODE_L6I;raw.code[2]=CODE_L1P;raw.code[3]=CODE_L5P;
    for (q=0;q<4;q++) {raw.P[q]=24000000+q;raw.L[q]=120000000+q;}
    for (i=0;i<4;i++) for (j=0;j<4;j++) {
        if (i==j) continue;
        o=options();strcpy(o.gps_if_pairs,"L1/L2");
        snprintf(o.bds_if_pairs,sizeof(o.bds_if_pairs),"%s/%s",names[i],names[j]);
        CHECK(if_parse_pairs(SYS_CMP,o.bds_if_pairs,codes)==1&&!codes[2]);
        CHECK(if_options_normalize(&o,msg));CHECK(o.nf==2&&o.if_model==2);
        CHECK(bds_options_valid(&o,msg));CHECK(NF(&o)==1&&ND(&o)==0);
        bds_decode_options_ex(decoder,sizeof(decoder),"",&o);
        CHECK(code2obsidx(SYS_CMP,codes[0],decoder)==0);
        CHECK(code2obsidx(SYS_CMP,codes[1],decoder)==1);
        bds_select_obs_ex(&raw,&selected,1,&o);bds_select_obs_ex(&selected,&again,1,&o);
        CHECK(!memcmp(&selected,&again,sizeof(selected)));
        CHECK(selected.code[0]==codes[0]&&selected.code[1]==codes[1]);
        for (k=0;k<4;k++) {
            if (k==i||k==j) continue;
            snprintf(pairs,sizeof(pairs),"%s/%s,%s/%s",names[i],names[j],names[i],names[k]);
            CHECK(if_parse_pairs(SYS_CMP,pairs,codes)==2);
        }
    }
    CHECK(!if_parse_pairs(SYS_CMP,"B1I/B1I",codes));
    CHECK(!if_parse_pairs(SYS_CMP,"B1I/B3I,B2a/B1C",codes));
    CHECK(!if_parse_pairs(SYS_CMP,"B1I/B3I,B1I/B3I",codes));
    CHECK(!if_parse_pairs(SYS_CMP,"B1I/B3I,",codes));
    CHECK(!if_parse_pairs(SYS_CMP,"B1I/B2b",codes));
    CHECK(!if_parse_pairs(SYS_GPS,"L1/L5",codes)); /* unimplemented clock datum */
    CHECK(if_parse_pairs(SYS_GAL,"E1/E5a,E1/E5b",codes)==2);
    CHECK(if_parse_pairs(SYS_CMP," b2a / b1c , b2a / b1i ",codes)==2);
    o=options();strcpy(o.gps_if_pairs,"L1/L2");strcpy(o.bds_if_pairs,"B2a/B1C,B2a/B1I");
    CHECK(!if_options_normalize(&o,msg)); /* different pair counts */
    strcpy(o.gps_if_pairs,"L1/L2,L1/L5");CHECK(if_options_normalize(&o,msg));
    CHECK(o.nf==3&&o.if_model==1);
    o.navsys|=SYS_GAL;CHECK(!if_options_normalize(&o,msg));CHECK(strstr(msg,"GAL")!=NULL);
    o=options();strcpy(o.gps_if_pairs,"L1/L2");strcpy(o.bds_if_pairs,"B1I/B3I");
    rtkinit(&rtk,&o);CHECK(rtk.opt.nf==2&&rtk.opt.if_model==2);
    CHECK(rtk.nx==pppnx(&rtk.opt));rtkfree(&rtk); /* normalize before state allocation */

    fp=fopen("unified_if_test.conf","w");CHECK(fp!=NULL);
    fprintf(fp,"prcopt.mode=7\nprcopt.sateph=5\nprcopt.ionoopt=3\nprcopt.modear=0\nprcopt.navsys=33\nprcopt.gps_if_pairs=L1/L2\nprcopt.bds_if_pairs=B1I/B3I\n");
    fclose(fp);load_config("unified_if_test.conf",&o,&sol,&fil);
    CHECK(o.nf==2&&o.if_model==2&&!strcmp(o.bds_freqs,"B1I,B3I"));
    CHECK(bds_options_valid(&o,msg));remove("unified_if_test.conf");
    resetsysopts();
    CHECK(str2opt(searchopt("pos1-gpsifpairs",sysopts),"L1/L2,L1/L5"));
    CHECK(str2opt(searchopt("pos1-bdsifpairs",sysopts),"B1I/B3I,B1I/B2a"));
    CHECK(str2opt(searchopt("pos1-navsys",sysopts),"33"));
    getsysopts(&o,NULL,NULL);CHECK(o.nf==3&&o.if_model==1);
    CHECK(!strcmp(o.bds_freqs,"B1I,B3I,B2a"));
    memset(pairs,'X',63);pairs[63]=0;CHECK(!if_parse_pairs(SYS_CMP,pairs,codes));
    {char longtext[80];memset(longtext,'X',79);longtext[79]=0;
     CHECK(!str2opt(searchopt("pos1-bdsifpairs",sysopts),longtext));}
    resetsysopts();getsysopts(&o,NULL,NULL);CHECK(!o.gps_if_pairs[0]&&!o.bds_if_pairs[0]);
}

static void measurement_regression(void)
{
    nav_t *nav=calloc(1,sizeof(*nav));prcopt_t o;obsd_t obs={0};
    meas_t actual;reference_ifobs_t ref;
    double rs[3]={26378137,12000000,15000000},rr[3]={6378137,0,0},e[3],az[2]={0.4,1.0};
    double ep[6]={2025,3,21,0,0,30},f,dantr[NFREQ],dants[NFREQ],L[NFREQ],P[NFREQ],Lc,Pc;
    int profile,geometry,i,j,c,sat;
    for (profile=0;profile<10;profile++) {
        o=options();o.posopt[0]=o.posopt[1]=0;
        if (profile<3) {o.if_model=0;o.nf=2;o.bds_if=profile==2?1:0;}
        if (profile==6) o.if_model=2;
        obs.time=epoch2time(ep);obs.sat=sat=satno(profile==0||profile==3?SYS_GPS:SYS_CMP,23);
        if (profile>=7) {
            o.if_model=0;o.bds_if=0;o.nf=3;
            o.sateph=profile==7?EPHOPT_PREC:profile==8?EPHOPT_SSRAPC:EPHOPT_BRDC;
            if (profile==9) o.ionoopt=IONOOPT_EST;
            obs.sat=sat=satno(profile==7?SYS_GAL:SYS_GPS,23);
        }
        memset(obs.code,0,sizeof(obs.code));memset(obs.P,0,sizeof(obs.P));memset(obs.L,0,sizeof(obs.L));
        if (profile==0||profile==3) {obs.code[0]=CODE_L1C;obs.code[1]=CODE_L2W;obs.code[2]=CODE_L5Q;}
        else if (profile==1) {obs.code[0]=CODE_L2I;obs.code[1]=CODE_L6I;}
        else if (profile==2) {obs.code[0]=CODE_L1P;obs.code[1]=CODE_L5P;}
        else {strcpy(o.bds_freqs,profile==5?"B1I,B3I,B2a":"B2a,B1C,B1I");CHECK(bds_parse_freqs(o.bds_freqs,obs.code));}
        if (profile>=7) {
            obs.code[0]=CODE_L1C;obs.code[1]=profile==7?CODE_L7Q:CODE_L2W;obs.code[2]=CODE_L5Q;
        }
        o.pcvr[0].var=calloc(NANTFREQ,sizeof(*o.pcvr[0].var));strcpy(o.pcvr[0].type,"TEST");
        for (i=0;i<NANTFREQ;i++) {
            o.pcvr[0].valid[i]=1;o.pcvr[0].off[i][0]=0.002*(i+1);
            o.pcvr[0].off[i][2]=0.03*(i+1);
        }
        nav->pcvs[sat-1]=o.pcvr[0];nav->pcvs[sat-1].sat=sat;
        nav->B2bssr[sat].t0[0]=nav->B2bssr[sat].t0[1]=nav->B2bssr[sat].t0[2]=obs.time;
        for (i=0;i<3;i++) if (obs.code[i]) {
            f=sat2freq(sat,obs.code[i],nav);
            obs.P[i]=24000000.0+10*SQR(FREQL1/f);
            obs.L[i]=(24000005.0-10*SQR(FREQL1/f))*f/CLIGHT;obs.SNR[i]=45000-i*1000;
            nav->B2bssr[sat].cbias_valid[obs.code[i]]=1;
            nav->B2bssr[sat].cbias[obs.code[i]]=obs.code[i]==CODE_L6I?0:0.23f*i;
        }
        geodist(rs,rr,e);
        for (geometry=0;geometry<2;geometry++) {
            corr_meas(&obs,nav,az,&o,geometry?rs:NULL,rr,e,0.08,0.17,&actual);
            if (o.if_model) {
                reference_configured_meas(&obs,nav,az,&o,geometry?rs:NULL,rr,e,0.08,0.17,&ref);
                for (i=0;i<NF(&o);i++) {
                    NEAR(actual.L[i],ref.L[i],1E-7);NEAR(actual.P[i],ref.P[i],1E-7);
                    for (j=0;j<NF(&o);j++) for (c=0;c<2;c++) NEAR(actual.cov[c][i][j],ref.cov[c][i][j],1E-12);
                }
            }
            else {
                memset(dantr,0,sizeof(dantr));memset(dants,0,sizeof(dants));
                if (geometry) {
                    if (o.bds_if) CHECK(antmodel_bds(o.pcvr,obs.code,o.antdel[0],az,0,0,dantr));
                    else antmodel(sat,o.pcvr,o.antdel[0],az,0,dantr);
                }
                reference_legacy_meas(&obs,nav,az,&o,dantr,geometry?rs:NULL,e,0.08,dants,0.17,L,P,&Lc,&Pc);
                for (i=0;i<o.nf;i++) {NEAR(actual.raw_L[i],L[i],1E-7);NEAR(actual.raw_P[i],P[i],1E-7);}
                NEAR(actual.L[0],Lc,1E-7);NEAR(actual.P[0],Pc,1E-7);
            }
        }
        if (profile==5) {
            meas_t triple=actual;
            strcpy(o.gps_if_pairs,"L1/L2");strcpy(o.bds_if_pairs,"B1I/B3I");
            CHECK(if_options_normalize(&o,(char[128]){0}));
            obs.code[2]=0;obs.P[2]=obs.L[2]=0;
            corr_meas(&obs,nav,az,&o,rs,rr,e,0.08,0.17,&actual);
            CHECK(actual.npair==1&&actual.P[1]==0&&actual.L[1]==0);
            NEAR(actual.P[0],triple.P[0],1E-7);NEAR(actual.L[0],triple.L[0],1E-7);
            NEAR(actual.cov[0][0][0],triple.cov[0][0][0],1E-12);
        }
        free(o.pcvr[0].var);
    }
    free(nav);
}

static void gal_model_tests(void)
{
    nav_t *nav=calloc(1,sizeof(*nav));prcopt_t o=options();rtk_t rtk;obsd_t obs={0};
    eph_t eph[2]={{0}};meas_t m;char msg[128];double ep[6]={2025,3,21,0,2,0};
    double rs[6]={26378137,12000000,15000000},rr[3]={6378137,0,0},dir[3],az[2]={0,1};
    double clocks[2]={0},satvar[1]={9},dr[3]={0},v[4],R[16],*H,a,b,f,range,bias;
    int i,j,n,svh[1]={0},exc[1]={0},sat=satno(SYS_GAL,4);
    o.navsys|=SYS_GAL;strcpy(o.gps_if_pairs,"L1/L2,L1/L5");
    strcpy(o.bds_if_pairs,"B2a/B1C,B2a/B1I");strcpy(o.gal_if_pairs,"E1/E5a,E1/E5b");
    strcpy(o.gal_navfile,"synthetic.rnx");
    CHECK(if_options_normalize(&o,msg));CHECK(!bds_options_valid(&o,msg));
    o.gal_ephemeris=1;CHECK(bds_options_valid(&o,msg));
    CHECK(ND(&o)==NSATGPS+1+NSATGAL);
    CHECK(IG(sat,&o)>IFB(&o)&&IG(sat,&o)<NR(&o));
    CHECK(signal_noise_index(SYS_GAL,CODE_L5Q)==2);
    CHECK(signal_noise_index(SYS_GAL,CODE_L7Q)==1);
    obs.time=epoch2time(ep);obs.sat=sat;obs.code[0]=CODE_L1C;obs.code[1]=CODE_L5Q;obs.code[2]=CODE_L7Q;
    nav->gal_eph=eph;nav->ngal_eph=2;
    for (i=0;i<2;i++) {
        eph[i].sat=sat;eph[i].A=29600000;eph[i].e=0.002;eph[i].i0=0.96;
        eph[i].toe=eph[i].toc=eph[i].ttr=timeadd(obs.time,-120);
        eph[i].code=1<<(i?9:8);
    }
    eph[0].tgd[0]=2E-9;eph[1].tgd[1]=-1E-9;
    CHECK(gal_select_eph(nav,sat,obs.time,1)==eph);
    CHECK(!gal_select_eph(nav,sat,timeadd(obs.time,MAXDTOE_GAL),1));
    eph[0].ttr=timeadd(obs.time,1);CHECK(!gal_select_eph(nav,sat,obs.time,1));eph[0].ttr=eph[0].toe;
    o.pcvr[0].var=calloc(NANTFREQ,sizeof(*o.pcvr[0].var));strcpy(o.pcvr[0].type,"TEST");
    o.pcvr[0].valid[3*NFREQ]=o.pcvr[0].valid[3*NFREQ+1]=o.pcvr[0].valid[3*NFREQ+2]=1;
    CHECK(signal_ant_index(o.pcvr,SYS_GAL,CODE_L7Q,0)==3*NFREQ+2);
    nav->pcvs[sat-1]=o.pcvr[0];nav->pcvs[sat-1].sat=sat;
    range=geodist(rs,rr,dir);if_coefficients(FREQL1,FREQE5b,&a,&b);
    for (i=0;i<3;i++) {
        CHECK(gal_code_bias(nav,sat,obs.time,obs.code[i],&bias));
        f=sat2freq(sat,obs.code[i],nav);
        obs.P[i]=range+10*SQR(FREQL1/f)+bias+(i==2?2.0/b:0);
        obs.L[i]=(range-10*SQR(FREQL1/f)+5)*f/CLIGHT;obs.SNR[i]=45000;
    }
    corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&m);
    NEAR(m.P[0],range,1E-7);NEAR(m.P[1],range+2,1E-7);
    NEAR(m.L[0],range+5,1E-7);NEAR(m.L[1],range+5,1E-7);
    eph[1].svh=1;corr_meas(&obs,nav,az,&o,NULL,NULL,NULL,0,0,&m);
    CHECK(m.P[0]!=0&&m.P[1]==0);eph[1].svh=0;
    o.tropopt=TROPOPT_SAAS;o.posopt[0]=o.posopt[1]=o.posopt[2]=0;o.maxinno[0]=o.maxinno[1]=0;
    rtkinit(&rtk,&o);rtk.sol.time=obs.time;memcpy(rtk.sol.rr,rr,sizeof(rr));
    for (i=0;i<3;i++) initx(&rtk,rr[i],1,i);
    initx(&rtk,2,1,IG(sat,&o));
    initx(&rtk,5,1,IB(sat,0,&o));initx(&rtk,5,1,IB(sat,1,&o));
    rtk.ssat[sat-1].vs=1;H=calloc(rtk.nx*4,sizeof(double));
    n=ppp_res(0,&obs,1,rs,clocks,satvar,svh,dr,exc,nav,rtk.x,&rtk,v,H,R,az);
    CHECK(n==4);NEAR(H[IG(sat,&o)+3*rtk.nx],1,1E-12);
    NEAR(H[IG(sat,&o)+2*rtk.nx],0,1E-12); /* code-only, unlike GPS link IFCB */
    NEAR(H[IFB(&o)+3*rtk.nx],0,1E-12);
    NEAR(v[0],v[2],1E-7);NEAR(v[1],v[3],1E-7);
    NEAR(R[1],9,1E-10);CHECK(R[2]>9); /* common orbit/clock also correlates phase/code */
    for (i=0;i<4;i++) for (j=0;j<4;j++) NEAR(R[i+4*j],R[j+4*i],1E-10);
    obs.L[2]=0;exc[0]=0;
    n=ppp_res(0,&obs,1,rs,clocks,satvar,svh,dr,exc,nav,rtk.x,&rtk,v,H,R,az);CHECK(n==2);
    /* No GAL B2b messages and no ordinary nav.eph: explicit product still works. */
    satposs_if(obs.time,&obs,1,nav,&o,rs,clocks,satvar,svh);
    CHECK(norm(rs,3)>1E7&&svh[0]==0&&satvar[0]>=9);
    o.gal_ephemeris=0;satposs_if(obs.time,&obs,1,nav,&o,rs,clocks,satvar,svh);CHECK(norm(rs,3)==0);
    free(H);rtkfree(&rtk);free(o.pcvr[0].var);nav->gal_eph=NULL;free(nav);
}

static void gal_navigation_tests(void)
{
    nav_t *nav=calloc(1,sizeof(*nav));rnxopt_t ro={0};prcopt_t o=options();
    eph_t eph={0};FILE *fp;char msg[128],path[MAXSTRPATH],text[MAXSTRPATH+1];
    double ep[6]={2025,3,21,0,0,0};gtime_t t=epoch2time(ep);int i;
    resetsysopts();
    CHECK(str2opt(searchopt("pos1-galifpairs",sysopts),"E1/E5a,E1/E5b"));
    CHECK(str2opt(searchopt("pos1-gpsifpairs",sysopts),"L1/L2,L1/L5"));
    CHECK(str2opt(searchopt("pos1-navsys",sysopts),"9"));
    CHECK(str2opt(searchopt("pos1-galeph",sysopts),"1"));
    CHECK(str2opt(searchopt("pos1-galbrdcsig",sysopts),"4.5"));
    CHECK(str2opt(searchopt("file-galnav",sysopts),"test.rnx"));
    getsysopts(&o,NULL,NULL);CHECK(o.nf==3&&o.if_model==1&&o.gal_ephemeris==1);
    NEAR(o.gal_brdc_sigma,4.5,1E-12);CHECK(!strcmp(o.gal_navfile,"test.rnx"));
    memset(text,'x',MAXSTRPATH);text[MAXSTRPATH]=0;
    CHECK(!str2opt(searchopt("file-galnav",sysopts),text));
    getsysopts(&o,NULL,NULL);CHECK(!o.gal_navfile[0]);resetsysopts();
    strcpy(o.gal_if_pairs,"E1/E5a");CHECK(!if_options_normalize(&o,msg));
    strcpy(o.gps_if_pairs,"L1/L2");CHECK(if_options_normalize(&o,msg));
    CHECK(o.nf==2&&o.if_model==2&&ND(&o)==0);
    strcpy(o.gal_navfile,"test.rnx");o.mode=PMODE_PPP_KINEMA;o.sateph=EPHOPT_B2b;
    o.ionoopt=IONOOPT_IFLC;o.modear=ARMODE_OFF;CHECK(bds_options_valid(&o,msg));
    o.navsys=SYS_GAL;CHECK(!bds_options_valid(&o,msg));
    /* Own a test-generated file; never overwrite an existing file. */
    strcpy(path,"if1213_synthetic_nav.rnx");
    if ((fp=fopen(path,"r"))) {fclose(fp);CHECK(0);free(nav);return;}
    fp=fopen(path,"w");CHECK(fp!=NULL);if (!fp) {free(nav);return;}
    ro.rnxver=304;ro.navsys=SYS_GAL;CHECK(outrnxnavh(fp,&ro,nav));
    eph.sat=satno(SYS_GAL,4);eph.toe=eph.toc=eph.ttr=t;
    eph.toes=time2gpst(t,&eph.week);eph.A=29600000;eph.i0=0.96;eph.e=0.002;
    eph.code=258;eph.tgd[0]=2E-9;CHECK(outrnxnavb(fp,&ro,&eph));
    eph.code=513;eph.tgd[1]=-1E-9;CHECK(outrnxnavb(fp,&ro,&eph));fclose(fp);
#ifdef _WIN32
    CHECK(_fullpath(text,path,sizeof(text))!=NULL);
    strcpy(path,text);for (i=0;path[i];i++) if (path[i]=='\\') path[i]='/';
#endif
    CHECK(gal_load_nav(nav,path,msg));CHECK(nav->ngal_eph==2);
    CHECK(gal_select_eph(nav,eph.sat,timeadd(t,1),1)!=NULL);
    CHECK(gal_select_eph(nav,eph.sat,timeadd(t,1),0)!=NULL);
    CHECK(!gal_select_eph(nav,eph.sat,timeadd(t,-1),1));
    CHECK(!gal_load_nav(nav,"nonexistent_GAL_test_input.rnx",msg));
    CHECK(nav->ngal_eph==2); /* failed reload leaves current products intact */
    CHECK(remove(path)==0);freenav(nav,1);CHECK(!nav->gal_eph&&!nav->ngal_eph);free(nav);
}

static void bds_prn_reassignment_tests(void)
{
    nav_t *nav=calloc(1,sizeof(*nav));
    obsd_t obs={0};prcopt_t o=options();
    double oldvar,newvar;
    obs.sat=satno(SYS_CMP,6);obs.code[0]=CODE_L2I;
    nav->pcvs[obs.sat-1].sat=obs.sat;
    strcpy(nav->pcvs[obs.sat-1].type,"BEIDOU-3I");
    CHECK(ppp_bds3(&obs,nav)); /* C06 after April 2026 */
    strcpy(nav->pcvs[obs.sat-1].type,"BEIDOU-2I");
    CHECK(!ppp_bds3(&obs,nav)); /* historical C06 */
    nav->pcvs[obs.sat-1].type[0]=0;obs.code[1]=CODE_L5P;
    CHECK(ppp_bds3(&obs,nav));
    obs.code[1]=0;CHECK(!ppp_bds3(&obs,nav));
    obs.sat=satno(SYS_CMP,39);
    newvar=varerr(obs.sat,SYS_CMP,PI/4,45,0,&o,&obs);
    NEAR(newvar,varerr(satno(SYS_CMP,21),SYS_CMP,PI/4,45,0,&o,&obs),1E-12);
    o.if_model=0;
    oldvar=varerr(obs.sat,SYS_CMP,PI/4,45,0,&o,&obs);
    CHECK(oldvar>newvar*35); /* legacy mode retains its historical weighting */
    free(nav);
}

int main(void)
{
    bds_prn_reassignment_tests();
    coefficients();selection();antenna();configurable_selection();pair_config_tests();measurement_regression();gal_model_tests();gal_navigation_tests();
    measurement_and_slip("B2a,B1C,B1I");measurement_and_slip("B1I,B3I,B2a");equations();
    printf("IF1213: %d checks, %d failures\n",checks,failures);
    return failures?EXIT_FAILURE:EXIT_SUCCESS;
}
