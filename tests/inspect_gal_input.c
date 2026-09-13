#include "rtklib.h"
PPPGlobal_t PPP_Glo={0};
int main(int argc,char **argv)
{
    rtcm_t *r=calloc(1,sizeof(*r));FILE *fp;int ret,i,j,counts[MAXCODE+1]={0},ep=0,nav=0;
    double date[6]={2025,3,21,0,0,0};gtime_t start=epoch2time(date);
    if (argc!=2||!(fp=fopen(argv[1],"rb"))||!init_rtcm(r)) return 1;
    r->time=start;strcpy(r->opt,"-EPHALL");
    while ((ret=input_rtcm3f(r,fp))!=-2) {
        if (ret==2&&satsys(r->ephsat,NULL)==SYS_GAL) {
            nav++;printf("GAL NAV sat=%d set=%d\n",r->ephsat,r->ephset);
        }
        if (ret==1) {
            ep++;
            for (i=0;i<r->obs.n;i++) if (satsys(r->obs.data[i].sat,NULL)==SYS_GAL)
                for (j=0;j<NFREQ+NEXOBS;j++) if (r->obs.data[i].P[j]&&r->obs.data[i].L[j])
                    counts[r->obs.data[i].code[j]]++;
            if (ep>=3600) break; /* 1 Hz rover stream; avoid decoder wall-clock week ambiguity */
        }
    }
    printf("epochs=%d GAL_eph=%d bytes=%ld\n",ep,nav,ftell(fp));
    for (j=1;j<=MAXCODE;j++) if (counts[j]) printf("GAL %s count=%d\n",code2obs(j),counts[j]);
    fclose(fp);free_rtcm(r);free(r);return 0;
}
