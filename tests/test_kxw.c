#include "rtklib.h"
PPPGlobal_t PPP_Glo={0};
static int checks=0;
#define CHECK(x) do {checks++;if (!(x)) {fprintf(stderr,"FAIL %d: %s\n",__LINE__,#x);exit(1);}} while(0)
#define CLOSE(a,b) CHECK(fabs((a)-(b))<1E-6)
static int packet(raw_t *r,int id,int info,const uint8_t *p,int n,int corrupt)
{
    uint8_t b[1030]={0};int i,ret=0;unsigned crc;
    b[0]=0xd3;b[1]=(n+4)>>8;b[2]=(n+4)&255;b[3]=0xfc;b[4]=0xf0;b[5]=id;b[6]=info;
    memcpy(b+7,p,n);crc=rtk_crc24q(b,n+7);setbitu(b,(n+7)*8,24,crc^(corrupt?1:0));
    for (i=0;i<n+10;i++) ret=input_raw(r,STRFMT_KXW,b[i]);
    return ret;
}
static void header(uint8_t *b,int type)
{
    memset(b,0,64);setbitu(b,0,6,type);setbitu(b,6,17,14000);setbitu(b,27,2,2);
}
static void payload(const uint8_t *b,uint8_t *p)
{
    int i,j;p[0]=62;p[1]=0;p[2]=getbitu(b,0,6);
    for (i=0;i<16;i++) for (j=0;j<4;j++) p[3+4*i+j]=b[4*i+3-j];
}
static int frame(raw_t *r,const uint8_t *b)
{
    uint8_t p[67];payload(b,p);return packet(r,64,0,p,67,0);
}
static void put2(uint8_t *p,int v) {p[0]=v&255;p[1]=(v>>8)&255;}
static void put4(uint8_t *p,int v) {int i;for(i=0;i<4;i++)p[i]=(v>>(8*i))&255;}
static void put8(uint8_t *p,double v) {uint64_t u;int i;memcpy(&u,&v,8);for(i=0;i<8;i++)p[i]=(uint8_t)(u>>(8*i));}
static void tag_tests(void)
{
    int width,poswidth,i;FILE *f;uint8_t h[64]={0},out[8];uint64_t sec=1788925237,pos;
    double fraction=.015;uint32_t tick;stream_t stream;char name[64],tag[80],url[128];gtime_t t;
    for(width=4;width<=8;width+=4) for(poswidth=4;poswidth<=8;poswidth+=4) {
        sprintf(name,"test_kxw_tag_%d_%d.bin",width,poswidth);sprintf(tag,"%s.tag",name);
        f=fopen(name,"wb");CHECK(f!=NULL);fwrite("ABCD",1,4,f);fclose(f);
        f=fopen(tag,"wb");CHECK(f!=NULL);memcpy(h,"TIMETAG RTKLIB 2.4.2",19);
        fwrite(h,1,64,f);fwrite(&sec,1,width,f);fwrite(&fraction,1,8,f);
        for(i=0;i<2;i++) {tick=i*1000;pos=2+2*i;fwrite(&tick,1,4,f);fwrite(&pos,1,poswidth,f);}
        fclose(f);strinit(&stream);sprintf(url,"%s::T::P=%d::+2",name,poswidth);
        CHECK(stropen(&stream,STR_FILE,STR_MODE_R,url));t=strgettime(&stream);
        CHECK(t.time==sec+2);CLOSE(t.sec,fraction);
        CHECK(strread(&stream,out,8)==4);CHECK(!memcmp(out,"ABCD",4));strclose(&stream);
        CHECK(remove(name)==0);CHECK(remove(tag)==0);
    }
    {
        pcv_t pcv={0};double az[2]={0,PI/2},del[3]={0,0,1.2},d;
        CHECK(signal_antmodel(&pcv,SYS_GPS,CODE_L5Q,del,az,0,0,0,&d));CLOSE(d,-1.2);
        CHECK(!signal_antmodel(&pcv,SYS_GPS,CODE_L5Q,del,az,0,1,0,&d));
    }
}
int main(void)
{
    raw_t *r=calloc(1,sizeof(*r)),*other=calloc(1,sizeof(*r));
    uint8_t b[64],p[178]={0};double ep[6]={2026,9,9,3,40,0},v[20]={0};
    int sat=satno(SYS_CMP,19),k,i;gtime_t saved;
    CHECK(init_raw(r,STRFMT_KXW));CHECK(init_raw(other,STRFMT_KXW));
    r->time=other->time=epoch2time(ep);
    /* A valid mask carried by C62 must not be hard-coded away. */
    header(b,1);setbitu(b,29,4,3);setbitu(b,33+18,1,1);setbitu(b,33+63,1,1);
    CHECK(frame(r,b)==20);CHECK(r->geoprn==62);
    header(b,2);k=29+69; /* empty first slot, followed by signed corrections */
    setbitu(b,k,9,19);setbitu(b,k+9,10,55);setbitu(b,k+19,3,5);
    setbits(b,k+22,15,-1250);setbits(b,k+37,13,100);setbits(b,k+50,13,-200);setbitu(b,k+63,6,9);
    CHECK(frame(other,b)==0); /* each stream needs its own mask */
    CHECK(frame(r,b)==20);CLOSE(r->nav.B2bssr[sat].deph[0],-2.0);
    CLOSE(r->nav.B2bssr[sat].deph[1],.64);CLOSE(r->nav.B2bssr[sat].deph[2],-1.28);
    CHECK(r->nav.B2bssr[sat].ura==9);CHECK(r->nav.B2bssr[sat].iodn==55);
    CHECK(!r->nav.B2bssr[0].update);
    header(b,3);setbitu(b,29,5,1);setbitu(b,34,9,19);setbitu(b,43,4,2);
    setbitu(b,47,4,0);setbits(b,51,12,-100);setbitu(b,63,4,1);setbits(b,67,12,0);
    CHECK(frame(r,b)==20);CLOSE(r->nav.B2bssr[sat].cbias[b2b_bds_codebias_mode[0]],-1.7);
    CHECK(r->nav.B2bssr[sat].cbias_valid[b2b_bds_codebias_mode[0]]);
    if(b2b_bds_codebias_mode[1]) CHECK(r->nav.B2bssr[sat].cbias_valid[b2b_bds_codebias_mode[1]]);
    header(b,4);setbitu(b,29,4,3);setbitu(b,38,3,5);setbits(b,41,15,-625);
    CHECK(frame(r,b)==20);CLOSE(r->nav.B2bssr[sat].dclk[0],-1.0);
    saved=r->time;header(b,63);setbitu(b,6,17,131071);CHECK(frame(r,b)==0);CLOSE(timediff(saved,r->time),0);
    header(b,4);setbitu(b,29,4,3);setbitu(b,33,5,31);CHECK(frame(r,b)==0);
    header(b,3);setbitu(b,29,5,31);for(i=0;i<2;i++)setbitu(b,43+i*253,4,15);
    /* Do not write beyond the synthetic frame: first two records already overflow. */
    CHECK(frame(r,b)==-1);
    header(b,1);setbitu(b,29,4,3);setbitu(b,33+18,1,1);payload(b,p);
    CHECK(packet(r,64,0,p,67,1)==-1); /* transport CRC */
    p[1]=0x20;CHECK(packet(r,64,0,p,67,0)==0);p[1]=0;
    CHECK(packet(r,64,1,p,20,0)==0);CHECK(packet(r,64,0x11,p+20,47,0)==20);
    CHECK(packet(r,64,0x11,p+20,47,0)==-1); /* orphan */
    CHECK(packet(r,64,2,p,20,0)==0);CHECK(packet(r,64,0x22,p+20,47,0)==-1); /* missing middle */
    CHECK(packet(r,64,1,p,20,0)==0);CHECK(packet(r,1,0,p,67,0)==-1); /* interleaved */
    CHECK(packet(r,64,0,p,67,0)==20); /* resync */
    /* CNAV1 unaligned fields, BDT epoch and absence of TGD. */
    memset(p,0,sizeof(p));p[0]=19;p[1]=7;p[2]=3;put2(p+4,55);put2(p+6,1079);put4(p+10,270000);put4(p+14,270000);
    v[0]=.0001;v[3]=.002;v[8]=.95;v[16]=400;v[17]=.01;v[19]=1E-14;
    for(i=0;i<20;i++)put8(p+18+8*i,v[i]);
    CHECK(packet(r,1,0,p,178,0)==2);CHECK(r->ephsat==sat);CHECK(r->nav.eph[sat-1].flag==BDS_NAV_KXW_CNV1);
    CHECK(r->nav.eph[sat-1].iode==7&&r->nav.eph[sat-1].iodc==55);CHECK(!r->nav.eph[sat-1].tgd_valid);
    CLOSE(r->nav.eph[sat-1].A,27906500.0);CLOSE(r->nav.eph[sat-1].Adot,.01);
    CLOSE(timediff(r->nav.eph[sat-1].toe,bdt2gpst(bdt2time(1079,270000))),0);
    CHECK(!PPP_Glo.BDS_CNV1_flag);CHECK(packet(r,1,0,p,178,0)==0);
    put2(p+6,-1);CHECK(packet(r,1,0,p,178,0)==-1);
    tag_tests();free_raw(r);free_raw(other);free(r);free(other);printf("KXW: %d checks passed\n",checks);return 0;
}
