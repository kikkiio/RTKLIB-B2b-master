/* Audit tagged RTCM3 observations with the same decoder used by rtppp. */
#include "rtklib.h"
PPPGlobal_t PPP_Glo={0};
int main(int argc,char **argv)
{
    rtcm_t *r;stream_t st;FILE *fp;char path[MAXSTRPATH],first[64],last[64];
    int ret,i,j,k,errors=0,epochs=0,nav=0,counts[2][MAXCODE+1]={{0}};
    gtime_t t0={0},t1={0};double gap,maxgap=0;
    if(argc!=2){fprintf(stderr,"usage: inspect_rtcm_input ROVER (requires ROVER.tag)\n");return 2;}
    if(strlen(argv[1])+4>=sizeof(path))return 2;
    r=(rtcm_t *)calloc(1,sizeof(*r));if(!r||!init_rtcm(r))return 2;
    snprintf(path,sizeof(path),"%s::T",argv[1]);strinit(&st);
    if(!stropen(&st,STR_FILE,STR_MODE_R,path))return 2;
    r->time=strgettime(&st);strclose(&st);strcpy(r->opt,"-EPHALL");
    if(!(fp=fopen(argv[1],"rb")))return 2;
    while((ret=input_rtcm3f(r,fp))!=-2){
        if(ret<0)errors++;
        if(ret==2)nav++;
        if(ret!=1||!r->obs.n)continue;
        if(!epochs)t0=r->obs.data[0].time;
        if(epochs){gap=timediff(r->obs.data[0].time,t1);if(gap>maxgap)maxgap=gap;}
        t1=r->obs.data[0].time;epochs++;
        for(i=0;i<r->obs.n;i++){
            int sys=satsys(r->obs.data[i].sat,NULL);
            k=sys==SYS_GPS?0:sys==SYS_CMP?1:-1;if(k<0)continue;
            for(j=0;j<NFREQ+NEXOBS;j++)if(r->obs.data[i].P[j]&&r->obs.data[i].L[j])
                counts[k][r->obs.data[i].code[j]]++;
        }
    }
    time2str(t0,first,3);time2str(t1,last,3);
    printf("{\"epochs\":%d,\"decode_errors\":%d,\"navigation_messages\":%d,",epochs,errors,nav);
    printf("\"first_gpst\":\"%s\",\"last_gpst\":\"%s\",\"span_seconds\":%.3f,\"max_gap_seconds\":%.3f,\"signals\":{",first,last,timediff(t1,t0),maxgap);
    for(k=0;k<2;k++){
        int comma=0;printf("%s\"%s\":{",k?",":"",k?"BDS":"GPS");
        for(j=1;j<=MAXCODE;j++)if(counts[k][j]){printf("%s\"%s\":%d",comma?",":"",code2obs(j),counts[k][j]);comma=1;}
        printf("}");
    }
    printf("}}\n");fclose(fp);free_rtcm(r);free(r);return epochs?0:1;
}
