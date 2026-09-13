#include "rtklib.h"

PPPGlobal_t PPP_Glo={0};
static int checks=0;
#define CHECK(x) do {checks++; if (!(x)) {fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x); exit(1);}} while (0)
#define NEAR(a,b,e) CHECK(fabs((a)-(b))<(e))

static void test_selection(void)
{
    obsd_t in={0},out,again;
    prcopt_t opt=prcopt_default;
    char options[256],msg[128];
    int i;
    in.sat=satno(SYS_CMP,19);
    in.code[0]=CODE_L2I; in.code[1]=CODE_L6I;
    in.code[3]=CODE_L1P; in.code[4]=CODE_L5P;
    for (i=0;i<NFREQ+NEXOBS;i++) {
        in.P[i]=24000000+i;in.L[i]=120000000+i;in.D[i]=i;
        in.SNR[i]=40000+i;in.LLI[i]=i%2;in.Lstd[i]=i;in.Pstd[i]=i+1;
    }
    bds_select_obs(&in,&out,1,0);
    CHECK(!memcmp(&in,&out,sizeof(in)));
    bds_select_obs(&in,&out,1,1);
    CHECK(out.code[0]==CODE_L1P&&out.code[1]==CODE_L5P);
    CHECK(out.P[0]==in.P[3]&&out.L[1]==in.L[4]);
    CHECK(out.SNR[1]==in.SNR[4]&&out.LLI[0]==in.LLI[3]);
    CHECK(out.D[0]==in.D[3]&&out.Pstd[1]==in.Pstd[4]);
    bds_select_obs(&out,&again,1,1);
    CHECK(!memcmp(&out,&again,sizeof(out)));
    in.code[4]=0;
    bds_select_obs(&in,&out,1,1);
    CHECK(!out.code[1]&&!out.P[1]&&!out.L[1]);
    in.code[4]=CODE_L5X;
    bds_select_obs(&in,&out,1,1);
    CHECK(!out.code[1]);
    in.sat=satno(SYS_GPS,1);
    bds_select_obs(&in,&out,1,1);
    CHECK(!memcmp(&in,&out,sizeof(in)));
    bds_decode_options(options,sizeof(options),"-EPHALL -CBDSIF=0",1);
    CHECK(code2obsidx(SYS_CMP,CODE_L1P,options)==0);
    CHECK(code2obsidx(SYS_CMP,CODE_L5P,options)==1);
    CHECK(code2obsidx(SYS_CMP,CODE_L2I,options)>=NFREQ);
    CHECK(code2obsidx(SYS_CMP,CODE_L6I,options)>=NFREQ);
    CHECK(code2obsidx(SYS_CMP,CODE_L2I,"")==0);
    CHECK(code2obsidx(SYS_CMP,CODE_L6I,"")==1);
    CHECK(code2obsidx(SYS_GPS,CODE_L2W,options)==1);
    CHECK(code2idx(SYS_CMP,CODE_L1P)==5); /* canonical frequency IDs unchanged */
    NEAR(code2freq(SYS_CMP,CODE_L1P,0),1575420000.0,0.01);
    NEAR(code2freq(SYS_CMP,CODE_L5P,0),1176450000.0,0.01);
    CHECK(getcodepri(SYS_CMP,CODE_L1P,"")>getcodepri(SYS_CMP,CODE_L1D,""));
    CHECK(getcodepri(SYS_CMP,CODE_L5P,"")>getcodepri(SYS_CMP,CODE_L5D,""));
    CHECK(bds_options_valid(&opt,msg));
    opt.bds_if=1; CHECK(!bds_options_valid(&opt,msg));
    opt.mode=PMODE_PPP_KINEMA; opt.sateph=EPHOPT_B2b; opt.ionoopt=IONOOPT_IFLC;
    CHECK(bds_options_valid(&opt,msg));
    opt.nf=3; CHECK(!bds_options_valid(&opt,msg));
    opt.nf=2;opt.bds_if=2;CHECK(!bds_options_valid(&opt,msg));
}

static void test_bias(void)
{
    B2bssr_t ssr={0};
    nav_t *nav=(nav_t *)calloc(1,sizeof(nav_t));
    eph_t eph[2]={{0}};
    double ep[]={2025,3,21,0,0,0},bias;
    gtime_t t=epoch2time(ep);
    CHECK(nav!=NULL);
    ssr.t0[1]=t;
    ssr.iodssr[0]=ssr.iodssr[1]=ssr.iodssr[2]=1;
    ssr.cbias[CODE_L1P]=0.0f;
    CHECK(!bds_code_bias(t,&ssr,CODE_L1P,&bias));
    ssr.cbias_valid[CODE_L1P]=1;
    CHECK(bds_code_bias(t,&ssr,CODE_L1P,&bias)); NEAR(bias,0.0,1e-12);
    ssr.cbias[CODE_L5P]=-0.561f;ssr.cbias_valid[CODE_L5P]=1;
    CHECK(bds_code_bias(t,&ssr,CODE_L5P,&bias)); NEAR(bias,-0.561,1e-6);
    CHECK(!bds_code_bias(timeadd(t,86401),&ssr,CODE_L5P,&bias));
    CHECK(!bds_code_bias(timeadd(t,-2),&ssr,CODE_L5P,&bias));
    ssr.iodssr[2]=0;CHECK(!bds_code_bias(t,&ssr,CODE_L5P,&bias));
    nav->eph=eph;nav->n=2;
    eph[0].sat=eph[1].sat=satno(SYS_CMP,19);
    eph[0].toe=timeadd(t,-7200);eph[1].toe=t;
    eph[0].tgd[2]=1e-8;eph[1].tgd[2]=0;
    eph[0].tgd_valid=eph[1].tgd_valid=(1<<2);
    CHECK(bds_tgd_bias(t,eph[0].sat,nav,CODE_L1P,&bias));NEAR(bias,0.0,1e-12);
    CHECK(!bds_tgd_bias(t,eph[0].sat,nav,CODE_L5P,&bias));
    CHECK(!bds_tgd_bias(t,eph[0].sat,nav,CODE_L1D,&bias));
    CHECK(!bds_tgd_bias(t,eph[0].sat,nav,CODE_L1X,&bias));
    eph[1].tgd_valid|=(1<<4);eph[1].tgd[4]=1e-9;
    CHECK(bds_tgd_bias(t,eph[0].sat,nav,CODE_L1D,&bias));NEAR(bias,CLIGHT*1e-9,1e-12);
    free(nav);
}

static void test_ant(void)
{
    nav_t *nav=(nav_t *)calloc(1,sizeof(nav_t));
    pcv_t pcv={0};
    uint8_t modern[]={CODE_L1P,CODE_L5P},legacy[]={CODE_L2I,CODE_L6I};
    double azel[]={0,PI/2},del[3]={0},dant[3],v1[3],v2[3],rs[]={-15600000,20100000,21000000};
    double ep[]={2025,3,21,0,0,0};
    int sat=satno(SYS_CMP,19);
    CHECK(nav!=NULL);
    pcv.var=(double (*)[1600])calloc(NANTFREQ*1600,sizeof(double));CHECK(pcv.var!=NULL);
    CHECK(bds_ant_index(&pcv,CODE_L1P,0)<0);
    pcv.valid[3*NFREQ]=pcv.valid[3*NFREQ+1]=1;
    CHECK(bds_ant_index(&pcv,CODE_L1P,0)<0);
    CHECK(bds_ant_index(&pcv,CODE_L1P,1)==3*NFREQ);
    CHECK(bds_ant_index(&pcv,CODE_L5P,1)==3*NFREQ+1);
    pcv.valid[ANT_B1C]=pcv.valid[ANT_B2A]=1;
    pcv.off[ANT_B1C][2]=1;pcv.off[ANT_B2A][2]=2;
    CHECK(antmodel_bds(&pcv,modern,del,azel,1,0,dant));
    NEAR(dant[0],-1.0,1e-12);NEAR(dant[1],-2.0,1e-12);
    pcv.sat=sat;pcv.off[2*NFREQ][2]=3;pcv.off[2*NFREQ+1][2]=4;
    nav->pcvs[sat-1]=pcv;
    satantoff1(epoch2time(ep),rs,sat,nav,modern,v1,v2);
    NEAR(norm(v1,3),1.0,1e-10);NEAR(norm(v2,3),2.0,1e-10);
    satantoff1(epoch2time(ep),rs,sat,nav,legacy,v1,v2);
    NEAR(norm(v1,3),3.0,1e-10);NEAR(norm(v2,3),4.0,1e-10);
    satantoff2(epoch2time(ep),rs,sat,nav,dant,1);
    NEAR(norm(dant,3),4.0,1e-10); /* B3 orbit reference MUST stay B3 */
    pcv.valid[ANT_B1C]=0;CHECK(bds_ant_index(&pcv,CODE_L1P,1)<0); /* no satellite fallback */
    free(pcv.var);free(nav);
}

static void test_antex(void)
{
    const char *file="bds_antenna_test.atx";
    const char *freqs[]={"   C01","   C05","   C06","   E05","   G05"};
    pcvs_t pcvs={0};FILE *fp=fopen(file,"w");int i;
    CHECK(fp!=NULL);
    fprintf(fp,"%-60sSTART OF ANTENNA\n","");
    fprintf(fp,"%-60sTYPE / SERIAL NO\n","TESTANT         NONE");
    fprintf(fp,"%8.1f%52sDAZI\n",0.0,"");
    fprintf(fp,"  %6.1f%6.1f%6.1f%40sZEN1 / ZEN2 / DZEN\n",0.0,90.0,45.0,"");
    for (i=0;i<5;i++) {
        fprintf(fp,"%-60sSTART OF FREQUENCY\n",freqs[i]);
        fprintf(fp,"%10.2f%10.2f%10.2f%30sNORTH / EAST / UP\n",0.0,0.0,1000.0*(i+1),"");
        fprintf(fp,"   NOAZI%8.2f%8.2f%8.2f%28s\n",10.0,20.0,30.0,"");
        fprintf(fp,"%-60sEND OF FREQUENCY\n",freqs[i]);
    }
    fprintf(fp,"%-60sEND OF ANTENNA\n","");fclose(fp);
    CHECK(readpcv(file,&pcvs));CHECK(pcvs.n==1);
    CHECK(pcvs.pcv[0].valid[ANT_B1C]&&pcvs.pcv[0].valid[ANT_B2A]);
    CHECK(pcvs.pcv[0].valid[ANT_GPS_L5]);
    NEAR(pcvs.pcv[0].off[ANT_GPS_L5][2],5.0,1e-12);
    CHECK(!pcvs.pcv[0].valid[NFREQ+1]); /* G05 must not overwrite R02 */
    NEAR(pcvs.pcv[0].off[ANT_B1C][2],1.0,1e-12);
    NEAR(pcvs.pcv[0].off[ANT_B2A][2],2.0,1e-12);
    NEAR(pcvs.pcv[0].off[2*NFREQ+1][2],3.0,1e-12);
    NEAR(pcvs.pcv[0].off[3*NFREQ+1][2],4.0,1e-12);
    NEAR(pcvs.pcv[0].var[ANT_B1C][0],0.01,1e-12);
    NEAR(pcvs.pcv[0].var[ANT_B2A][2],0.03,1e-12);
    free(pcvs.pcv[0].var);free(pcvs.pcv);remove(file);
}

static void test_config(void)
{
    const char *file="bds_options_test.conf";
    FILE *fp=fopen(file,"w");
    prcopt_t opt;solopt_t sol;filopt_t fil;char msg[128];
    CHECK(fp!=NULL);
    fprintf(fp,"prcopt.mode=7\nprcopt.sateph=5\nprcopt.ionoopt=3\nprcopt.nf=2\nprcopt.bds_if=1\nprcopt.bds_ant_fallback=1\n");
    fclose(fp);load_config(file,&opt,&sol,&fil);
    CHECK(opt.bds_if==1&&opt.bds_ant_fallback==1);CHECK(bds_options_valid(&opt,msg));
    fp=fopen(file,"w");CHECK(fp!=NULL);fprintf(fp,"prcopt.mode=7\n");fclose(fp);
    load_config(file,&opt,&sol,&fil);CHECK(opt.bds_if==0&&opt.bds_ant_fallback==0);
    resetsysopts();CHECK(str2opt(searchopt("pos1-bdsif",sysopts),"b1c-b2a"));
    getsysopts(&opt,NULL,NULL);CHECK(opt.bds_if==1);
    fp=fopen(file,"w");CHECK(fp!=NULL);
    fprintf(fp,"prcopt.mode=7\nprcopt.sateph=5\nprcopt.ionoopt=3\nprcopt.nf=3\nprcopt.modear=0\nprcopt.navsys=33\nprcopt.bds_if=2\nprcopt.if_model=1\nprcopt.bds_freqs=B1I,B3I,B2a\n");
    fclose(fp);load_config(file,&opt,&sol,&fil);
    CHECK(!strcmp(opt.bds_freqs,"B1I,B3I,B2a"));CHECK(bds_options_valid(&opt,msg));
    resetsysopts();CHECK(str2opt(searchopt("pos1-bdsfreqs",sysopts),"B1I,B3I,B2a"));
    getsysopts(&opt,NULL,NULL);CHECK(!strcmp(opt.bds_freqs,"B1I,B3I,B2a"));
    resetsysopts();getsysopts(&opt,NULL,NULL);CHECK(!opt.bds_freqs[0]);
    remove(file);
}

static void test_rtcm(void)
{
    rtcm_t *encoder=(rtcm_t *)calloc(1,sizeof(rtcm_t));
    rtcm_t *decoder=(rtcm_t *)calloc(1,sizeof(rtcm_t));
    double ep[]={2025,3,21,0,0,0};
    uint8_t codes[]={CODE_L2I,CODE_L6I,CODE_L1P,CODE_L5P};
    int i,j,pair,ret=0;
    CHECK(encoder&&decoder);CHECK(init_rtcm(encoder));
    encoder->time=epoch2time(ep);encoder->staid=1;encoder->obs.n=1;
    memset(encoder->obs.data,0,sizeof(obsd_t));
    encoder->obs.data[0].time=encoder->time;
    encoder->obs.data[0].sat=satno(SYS_CMP,19);
    for (i=0;i<4;i++) {
        double f=code2freq(SYS_CMP,codes[i],0);
        encoder->obs.data[0].code[i]=codes[i];
        encoder->obs.data[0].P[i]=24000000.0+i;
        encoder->obs.data[0].L[i]=(24000000.0+i)/CLIGHT*f;
        encoder->obs.data[0].SNR[i]=45000;
    }
    CHECK(gen_rtcm3(encoder,1127,0,0));
    for (pair=0;pair<4;pair++) {
        CHECK(init_rtcm(decoder));decoder->time=encoder->time;
        if (pair==3) {
            prcopt_t o=prcopt_default;o.bds_if=2;strcpy(o.bds_freqs,"B1I,B3I,B2a");
            bds_decode_options_ex(decoder->opt,sizeof(decoder->opt),"",&o);
        }
        else bds_decode_options(decoder->opt,sizeof(decoder->opt),"",pair);
        for (j=0;j<encoder->nbyte;j++) ret=input_rtcm3(decoder,encoder->buff[j]);
        CHECK(ret==1&&decoder->obs.n==1);
        CHECK(strstr(decoder->opt,"-CBDSIF=")!=NULL);
        CHECK(decoder->obs.data[0].code[0]==(pair==2?CODE_L5P:pair==1?CODE_L1P:CODE_L2I));
        CHECK(decoder->obs.data[0].code[1]==(pair==2?CODE_L1P:pair==1?CODE_L5P:CODE_L6I));
        NEAR(decoder->obs.data[0].P[0],24000000.0+(pair==2?3:pair==1?2:0),0.002);
        if (pair==2) CHECK(decoder->obs.data[0].code[2]==CODE_L2I);
        if (pair==3) CHECK(decoder->obs.data[0].code[2]==CODE_L5P);
        free_rtcm(decoder);
    }
    free_rtcm(encoder);free(encoder);free(decoder);
}

static void test_rinex(void)
{
    const char *file="bds_selection_test.obs";
    FILE *fp=fopen(file,"w");
    obs_t obs={0};nav_t *nav=(nav_t *)calloc(1,sizeof(nav_t));sta_t sta={0};
    gtime_t zero={0};int i;
    CHECK(fp&&nav);
    fprintf(fp,"%9.2f           OBSERVATION DATA    M                   RINEX VERSION / TYPE\n",3.04);
    fprintf(fp,"%-60sSYS / # / OBS TYPES\n","C    8 C2I L2I C6I L6I C1P L1P C5P L5P");
    fprintf(fp,"%-60sEND OF HEADER\n","");
    fprintf(fp,"> 2025 03 21 00 00  0.0000000  0  1\nC19");
    for (i=0;i<8;i++) fprintf(fp,"%14.3f  ",24000000.0+i);
    fprintf(fp,"\n");fclose(fp);
    CHECK(readrnxt(file,1,zero,zero,0,"-CBDSIF=1",&obs,nav,&sta)>0);
    CHECK(obs.n==1);CHECK(obs.data[0].code[0]==CODE_L1P&&obs.data[0].code[1]==CODE_L5P);
    NEAR(obs.data[0].P[0],24000004.0,1e-6);NEAR(obs.data[0].L[1],24000007.0,1e-6);
    free(obs.data);memset(&obs,0,sizeof(obs));
    CHECK(readrnxt(file,1,zero,zero,0,"-CBDSIF=2",&obs,nav,&sta)>0);
    CHECK(obs.n==1);CHECK(obs.data[0].code[0]==CODE_L5P&&obs.data[0].code[1]==CODE_L1P&&obs.data[0].code[2]==CODE_L2I);
    NEAR(obs.data[0].P[0],24000006.0,1e-6);NEAR(obs.data[0].L[2],24000001.0,1e-6);
    free(obs.data);memset(&obs,0,sizeof(obs));
    CHECK(readrnxt(file,1,zero,zero,0,"-CBDSIF=2 -CBDSFREQ=261",&obs,nav,&sta)>0);
    CHECK(obs.n==1);CHECK(obs.data[0].code[0]==CODE_L2I&&obs.data[0].code[1]==CODE_L6I&&obs.data[0].code[2]==CODE_L1P);
    free(obs.data);memset(&obs,0,sizeof(obs));
    CHECK(readrnxt(file,1,zero,zero,0,"-CBDSIF=2 -CBDSFREQ=265",&obs,nav,&sta)>0);
    CHECK(obs.n==1);CHECK(obs.data[0].code[0]==CODE_L2I&&obs.data[0].code[1]==CODE_L6I&&obs.data[0].code[2]==CODE_L5P);
    NEAR(obs.data[0].P[2],24000006.0,1e-6);NEAR(obs.data[0].L[1],24000003.0,1e-6);
    free(obs.data);memset(&obs,0,sizeof(obs));
    CHECK(readrnxt(file,1,zero,zero,0,"",&obs,nav,&sta)>0);
    CHECK(obs.n==1);CHECK(obs.data[0].code[0]==CODE_L2I&&obs.data[0].code[1]==CODE_L6I);
    free(obs.data);free(nav);remove(file);
}

int main(void)
{
    test_selection();test_bias();test_ant();test_antex();test_config();test_rinex();test_rtcm();
    printf("PASS: %d BDS signal/bias/antenna/RINEX checks\n",checks);
    return 0;
}
