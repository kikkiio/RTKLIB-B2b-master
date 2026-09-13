/* Physical-signal helpers shared by IF1213 observations and unit tests.
 * No changes to code2freq(), code2idx() or the PPP-B2b product APC datum. */
#include "rtklib.h"

extern int if_coefficients(double f1, double f2, double *a, double *b)
{
    double d=f1*f1-f2*f2;
    *a=*b=0.0;
    if (!isfinite(f1)||!isfinite(f2)||f1<=0.0||f2<=0.0||fabs(f1-f2)<1.0) return 0;
    *a=f1*f1/d; *b=-f2*f2/d;
    return 1;
}

extern double if_covariance(const double *a, const double *b, const double *rawvar)
{
    int i;
    double v=0.0;
    for (i=0;i<3;i++) v+=a[i]*b[i]*rawvar[i];
    return v;
}

extern int signal_ant_index(const pcv_t *pcv, int sys, uint8_t code, int fallback)
{
    const char *s=code2obs(code);
    int idx=-1;
    if (sys==SYS_GPS) {
        idx=s[0]=='1'?0:s[0]=='2'?1:s[0]=='5'?ANT_GPS_L5:-1;
    }
    else if (sys==SYS_CMP) {
        idx=s[0]=='5'?ANT_B2A:s[0]=='1'?ANT_B1C:s[0]=='2'?2*NFREQ:
            s[0]=='6'?2*NFREQ+1:-1;
    }
    else if (sys==SYS_GAL) idx=s[0]=='1'?3*NFREQ:s[0]=='5'?3*NFREQ+1:s[0]=='7'?3*NFREQ+2:-1;
    if (idx<0) return -1;
    if (pcv->valid[idx]) return idx;
    if (!pcv->sat&&fallback&&(s[0]=='1'||s[0]=='5')) {
        idx=3*NFREQ+(s[0]=='1'?0:1); /* explicit same-frequency E01/E05 */
        if (pcv->valid[idx]) return idx;
    }
    return -1;
}

/* Weight/SNR groups remain physical when observation slots are permuted:
 * L1/B1I/B1C -> 0; L2/B3I -> 1; L5/B2a -> 2. */
extern int signal_noise_index(int sys, uint8_t code)
{
    const char *s=code2obs(code);
    if (sys==SYS_CMP) return s[0]=='2'||s[0]=='1'?0:s[0]=='6'?1:s[0]=='5'?2:-1;
    if (sys==SYS_GPS) return s[0]=='1'?0:s[0]=='2'?1:s[0]=='5'?2:-1;
    if (sys==SYS_GAL) return s[0]=='1'?0:s[0]=='5'?2:s[0]=='7'?1:-1;
    return -1;
}

extern int signal_satantoff(gtime_t time, const double *rs, int sat,
                           const nav_t *nav, uint8_t code, double *dant)
{
    const pcv_t *pcv=nav->pcvs+sat-1;
    int i,idx=signal_ant_index(pcv,satsys(sat,NULL),code,0);
    double ex[3],ey[3],ez[3],es[3],r[3],sun[3],erpv[5]={0};
    for (i=0;i<3;i++) dant[i]=0.0;
    if (idx<0) return 0;
    sunmoonpos(gpst2utc(time),erpv,sun,NULL,NULL);
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return 0;
    for (i=0;i<3;i++) r[i]=sun[i]-rs[i];
    if (!normv3(r,es)) return 0;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return 0;
    cross3(ey,ez,ex);
    for (i=0;i<3;i++) dant[i]=pcv->off[idx][0]*ex[i]+pcv->off[idx][1]*ey[i]+pcv->off[idx][2]*ez[i];
    return 1;
}
