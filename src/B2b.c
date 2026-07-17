#include "rtklib.h"
#include "B2b.h"

#define EPSILON 1e-3


/* mjd_time for B2b tod to gpst */ 

/* constants/macros ----------------------------------------------------------*/
#define MAXLEAPS    64       /* max number of leap seconds table */

/* GPS time reference (1980-01-06 00:00:00) */
static const double gpst0[] = { 1980, 1, 6, 0, 0, 0 };
/* BDS time reference (2006-01-01 00:00:00) */
static const double bdst0[] = { 2006, 1, 1, 0, 0, 0 };
/* Leap seconds table (y, m, d, h, m, s, UTC-GPST) */
static double leaps[MAXLEAPS + 1][7] = {
    { 2017, 1, 1, 0, 0, 0, -18 },
    { 2015, 7, 1, 0, 0, 0, -17 },
    { 2012, 7, 1, 0, 0, 0, -16 },
    { 2009, 1, 1, 0, 0, 0, -15 },
    { 2006, 1, 1, 0, 0, 0, -14 },
    { 1999, 1, 1, 0, 0, 0, -13 },
    { 1997, 7, 1, 0, 0, 0, -12 },
    { 1996, 1, 1, 0, 0, 0, -11 },
    { 1994, 7, 1, 0, 0, 0, -10 },
    { 1993, 7, 1, 0, 0, 0, -9 },
    { 1992, 7, 1, 0, 0, 0, -8 },
    { 1991, 1, 1, 0, 0, 0, -7 },
    { 1990, 1, 1, 0, 0, 0, -6 },
    { 1988, 1, 1, 0, 0, 0, -5 },
    { 1985, 7, 1, 0, 0, 0, -4 },
    { 1983, 7, 1, 0, 0, 0, -3 },
    { 1982, 7, 1, 0, 0, 0, -2 },
    { 1981, 7, 1, 0, 0, 0, -1 },
    { 0 }
};

/* Add a time difference (in seconds) to an MJD-based time */
static mjd_gtime_t mjd_time_add(mjd_gtime_t t0, double dt) {
    mjd_gtime_t t1 = t0;
    t1.sod = t0.sod + dt;

    while (1) {
        if (t1.sod >= 86400.0) {
            t1.sod -= 86400.0;
            t1.mjd += 1;
        } else if (t1.sod < 0.0) {
            t1.sod += 86400.0;
            t1.mjd -= 1;
        } else {
            break;
        }
    }
    return t1;
}

/* Compute the difference (in seconds) between two MJD-based times */
static double mjd_time_diff(mjd_gtime_t t1, mjd_gtime_t t0) {
    return (t1.mjd - t0.mjd) * 86400.0 + (t1.sod - t0.sod);
}

/* Convert seconds of day (SOD) to hours, minutes, and seconds */
static void sod_to_hms(double sod, int *hh, int *minu, double *sec) {
    int hh_tmp = (int)floor(sod / 3600.0 + 1.0e-9);
    if (hh) *hh = hh_tmp;

    double tmp = sod - hh_tmp * 3600.0;
    int minu_tmp = (int)floor(tmp / 60.0 + 1.0e-9);
    if (minu) *minu = minu_tmp;

    double sec_tmp = sod - hh_tmp * 3600.0 - minu_tmp * 60.0;
    if (sec) *sec = sec_tmp;
} /* end of sod_to_hms */

/* Convert hours, minutes, and seconds to seconds of day (SOD) */
static double hms_to_sod(int hh, int minu, double sec) {
    return hh * 3600.0 + minu * 60.0 + sec;
}

/* Convert a two-digit year to a four-digit year */
int yy_to_yyyy(int yy) {
    int yyyy = yy;
    if (yyyy <= 50) yyyy += 2000;
    else if (yy > 50 && yy < 1900) yyyy += 1900;
    return yyyy;
}

/* Convert year, month, day, hour, minute, and second to MJD-based time */
static mjd_gtime_t ymdhms_to_mjd_time(const double *date) {
    mjd_gtime_t tt = {0};
    tt.sod = hms_to_sod((int)floor(date[3]), (int)floor(date[4]), date[5]);

    int yyyy = yy_to_yyyy((int)floor(date[0]));
    int month = (int)floor(date[1]);
    int day = (int)floor(date[2]);
    if (month <= 2) {
        yyyy -= 1;
        month += 12;
    }

    double jd = (int)floor(365.25 * yyyy + 1.0e-9) + (int)floor(30.6001 * (month + 1) + 1.0e-9) + day + 1720981.5;
    tt.mjd = (int)floor(jd - 2400000.5);
    return tt;
}

/* Convert year, day of year (DOY), and seconds of day (SOD) to MJD-based time */
mjd_gtime_t yrdoy_to_mjd_time(int year, int doy, double sod) {
    /* Time at yyyy-01-01 00:00:00 */
    double date[6] = { 0 };
    date[0] = year; date[1] = 1; date[2] = 1;
    date[3] = 0;    date[4] = 0; date[5] = 0;

    /* Convert January 1st of the given year to MJD-based time */
    mjd_gtime_t t_jan1 = ymdhms_to_mjd_time(date);

    /* Calculate total seconds from January 1st to the specified DOY plus SOD */
    double dt = (doy - 1) * 86400.0 + sod;

    /* Add the calculated seconds to get the final MJD-based time */
    mjd_gtime_t tt = mjd_time_add(t_jan1, dt);
    return tt;
}

/* Convert GPS time to BDS time (subtract 14 seconds) */
mjd_gtime_t gpst_to_bdst(mjd_gtime_t tt_gps) {
    return mjd_time_add(tt_gps, -14.0);
}

/* Convert BDS time to GPS time (add 14 seconds) */
mjd_gtime_t bdst_to_gpst(mjd_gtime_t tt_bds) {
    return mjd_time_add(tt_bds, 14.0);
}

/* Convert MJD-based time to GPS week and seconds of week (SOW) */
int mjd_time_to_gpst(mjd_gtime_t tt, int *week, double *sow) {
    /* Compute the GPS epoch start time (1980-01-06 00:00:00) */
    mjd_gtime_t t_gpst0 = ymdhms_to_mjd_time(gpst0);

    /* Compute the time difference in days */
    double delta_day = mjd_time_diff(tt, t_gpst0) / 86400.0;

    /* Compute the GPS week */
    int week_tmp = (int)floor(delta_day / 7.0 + 1.0e-9);
    if (week) *week = week_tmp;

    /* Compute the day of week (days since the last GPS week start) */
    int dow = (int)floor(delta_day - week_tmp * 7 + 1.0e-9);

    /* Compute the seconds of week (SOW) */
    if (sow) *sow = dow * 24.0 * 3600.0 + tt.sod;

    return dow;
} /* end of mjd_time_to_gpst */

/* ---------------------------- */ 
static char *rtklib_codes[]={       /* observation code strings */

    ""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
    "1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
    "2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
    "6A","6B","6C","6X","6Z", "6S","6L","8I","8Q","8X", /* 30-39 */
    "2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
    "5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
    "6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
};

static void merge_mask_arrays(B2bmask_t *mask, int *merged_array) {

    int index = 0;
    
    memcpy(merged_array + index, mask->MASK_BD, sizeof(mask->MASK_BD));
    index += 63;
    
    memcpy(merged_array + index, mask->MASK_GPS, sizeof(mask->MASK_GPS));
    index += 37;
    
    memcpy(merged_array + index, mask->MASK_GALILEO, sizeof(mask->MASK_GALILEO));
    index += 37;
    
    memcpy(merged_array + index, mask->MASK_GLONASS, sizeof(mask->MASK_GLONASS));
    index += 37;
}

/* initialize B2b control -----------------------------------------------------
* initialize B2b control struct and reallocate memory for related buffers
* args   : B2b_t *b2b      IO  B2b control struct
* return : status (1:ok,0:memory allocation error)
*-----------------------------------------------------------------------------*/
extern int init_B2b(B2b_t *b2b)
{
    gtime_t time0 = {0}; // Initialize time to zero
    int i, j;

    trace(3, "init_B2b:\n");

    // Initialize main structure members
    b2b->geoprn = 0;
    b2b->time = b2b->time_s = time0;
    b2b->num_PPPB2BINF01 = b2b->num_PPPB2BINF02 = b2b->num_PPPB2BINF03 = b2b->num_PPPB2BINF04 = 0;
    b2b->outtype = 0;
    b2b->msgtype[0] = '\0';
    b2b->nbyte = b2b->nbit = b2b->len = 0;
    b2b->word = 0;

    // Initialize message buffers and counts
    for (i = 0; i < 100; i++) b2b->nmsg2[i] = 0;
    for (i = 0; i < 400; i++) b2b->nmsg3[i] = 0;
    for (i = 0; i < 1200; i++) b2b->buff[i] = 0;

    // Initialize B2bssr array
    for (i = 0; i < MAXSAT; i++) {
        b2b->B2bssr[i] = (B2bssr_t){0}; // Initialize B2bssr_t elements (you can modify based on structure)
    }

    // Optionally, reallocate memory for dynamic buffers if necessary (e.g. for observation data)
    // if (!(b2b->some_dynamic_buffer = malloc(sizeof(SomeType) * SIZE))) {
    //     free_B2b(b2b);
    //     return 0;
    // }

    return 1; // Return success
}

extern int slot2satno(int slot)
{
	int prn = 0, sys = 0;

	if (slot >= B2B_BDS_MINSAT && slot <= B2B_BDS_MAXSAT) {
		sys = SYS_CMP; prn = slot - B2B_BDS_MINSAT + 1;
	}
	else if (slot >= B2B_GPS_MINSAT && slot <= B2B_GPS_MAXSAT) {
		sys = SYS_GPS; prn = slot - B2B_GPS_MINSAT + 1;
	}
	else if (slot >= B2B_GAL_MINSAT && slot <= B2B_GAL_MAXSAT) {
		sys = SYS_GAL; prn = slot - B2B_GAL_MINSAT + 1;
	}
	else if (slot >= B2B_GLO_MINSAT && slot <= B2B_GLO_MAXSAT) {
		sys = SYS_GLO; prn = slot - B2B_GLO_MINSAT + 1;
	}
	else return 0;

	return satno(sys, prn);
}

extern int mask2satno(B2bmask_t *B2bmask)
{
    int i, satslot, satno, effect_num = 0;
    int total_mask[B2B_MAXSAT] = {0};

    // Merge mask arrays
    merge_mask_arrays(B2bmask, total_mask);

    // Find valid satellites and store their numbers in B2bmask->satno array
    for (i = 0; i < B2B_MAXSAT; i++) {
        if (total_mask[i] == 1) {
            satslot = i + 1;  // Assuming satellite numbers start from 1
            satno = slot2satno(satslot);
            B2bmask->satno[effect_num++] = satno;  // Store valid satellite number in B2bmask->satno
        }
    }

    // Update the number of valid satellites
    B2bmask->satnum = effect_num;

    return effect_num;
}

/* output Type 4 before a matching Type 1 mask is available -----------------
 * Type 4 addresses satellites by their ordinal position in the Type 1 mask.
 * Preserve the decoded values as INDEX records instead of inventing PRNs. */
static void output_B2bClockIndex(raw_t *raw, gtime_t epoch, double udi,
                                 int iodssr, int iodp, int subtype,
                                 const int *iodcorr, const int *c0)
{
    double ep[6],rp[6];
    int i,n=0;

    for (i=0;i<23;i++) {
        if (c0[i]>-16383&&(c0[i]!=0||iodcorr[i]!=0)) n++;
    }
    if (n<=0) return;
    time2epoch(epoch,ep);
    time2epoch(raw->time,rp);
    B2b_trace(22,"> CLOCK_INDEX %04d %02d %02d %02d %02d %.1f %d %d "
              "B2bSatPrn_%d IODSSR_%d IODP_%d SUBTYPE_%d @ "
              "%04d %02d %02d %02d %02d %.1f\n",
              (int)ep[0],(int)ep[1],(int)ep[2],(int)ep[3],(int)ep[4],ep[5],
              (int)udi,n,raw->geoprn,iodssr,iodp,subtype,
              (int)rp[0],(int)rp[1],(int)rp[2],(int)rp[3],(int)rp[4],rp[5]);
    for (i=0;i<23;i++) {
        if (c0[i]<=-16383||(c0[i]==0&&iodcorr[i]==0)) continue;
        B2b_trace(22,"INDEX_%03d %5d %10.4f\n",
                  subtype*23+i+1,iodcorr[i],c0[i]*0.0016);
    }
}

/* decode one standard PPP-B2b information block ----------------------------
 * The receiver-specific frame, preamble, CRC and LDPC symbols are deliberately
 * outside this function. Its caller has already checked the standard CRC24Q;
 * this function consumes PRN/status/type and the 456 information bits defined
 * by the PPP-B2b SIS ICD.
 *
 * return : message type (1-4), 0 if the page belongs to an obsolete mask,
 *         -1 on malformed input
 *---------------------------------------------------------------------------*/
extern int decode_B2b(raw_t *raw, B2bmask_t *mask, int bitpos, int endbit)
{
    const uint8_t *p;
    B2bssr_t *ssr;
    gtime_t epoch, old_mask_time;
    gtime_t last_time[4];
    double ep[6], udi;
    uint32_t prn, type, sow, iodssr, iodp=0, nsat, nsig, mode;
    int i, j, pos, body_end, slot, sat, sys, code, v1, v2, v3;
    int clock_iodcorr[23],clock_c0[23];
    int *codes=NULL;

    if (!raw||!mask||bitpos<0||endbit<bitpos+18+456) return -1;
    p=raw->buff;
    pos=bitpos;
    prn=getbitu(p,pos,6); pos+=6;
    pos+=6;                              /* PPP service status/reserved */
    type=getbitu(p,pos,6); pos+=6;
    body_end=pos+456;
    if (prn<1||prn>63||type<1||type>4||body_end>endbit) return -1;

    raw->geoprn=(int)prn;
    sow=getbitu(p,pos,17); pos+=17;
    pos+=4;                              /* SSR update interval index */
    iodssr=getbitu(p,pos,2); pos+=2;
    epoch=B2btod2time(raw->time,(double)sow);
    /* Keep the raw stream time aligned with the correction epoch. Recorded
     * 4047/64 frames have no receiver header time of their own; the server
     * uses raw->time to reject corrections older than the rover solution. */
    raw->time=epoch;
    time2epoch(epoch,ep);

    if (type==1) {                       /* satellite mask */
        old_mask_time=mask->time;
        memcpy(last_time,mask->last_time,sizeof(last_time));
        iodp=getbitu(p,pos,4); pos+=4;
        memset(mask,0,sizeof(*mask));
        memcpy(mask->last_time,last_time,sizeof(last_time));
        mask->m_time=raw->time;
        mask->time=epoch;
        mask->IOD_SSR=(int)iodssr;
        mask->IODP=(int)iodp;
        for (i=0;i<63;i++) {mask->MASK_BD[i]=(int)getbitu(p,pos,1); pos++;}
        for (i=0;i<37;i++) {mask->MASK_GPS[i]=(int)getbitu(p,pos,1); pos++;}
        for (i=0;i<37;i++) {mask->MASK_GALILEO[i]=(int)getbitu(p,pos,1); pos++;}
        for (i=0;i<37;i++) {mask->MASK_GLONASS[i]=(int)getbitu(p,pos,1); pos++;}
        mask2satno(mask);
        raw->num_PPPB2BINF01++;
        raw->raw_nmsg[0]++;
        output_B2bInfo1(raw,mask,old_mask_time.time?
                        (int)timediff(epoch,old_mask_time):0);
        mask->last_time[0]=epoch;
        return 1;
    }
    /* Type 2/3 contain explicit SatSlot values and can be inspected before
     * Type 1 arrives. Once a mask exists, enforce the IOD SSR association. */
    if (mask->IOD_SSR>=0&&(int)iodssr!=mask->IOD_SSR) return 0;

    if (type==2) {                       /* orbit and URA: six records */
        udi=mask->last_time[1].time?timediff(epoch,mask->last_time[1]):0.0;
        if (udi<0.0||udi>86400.0) udi=0.0;
        for (i=0;i<6;i++) {
            slot=(int)getbitu(p,pos,9); pos+=9;
            v1=(int)getbitu(p,pos,10); pos+=10;       /* IODN */
            v2=(int)getbitu(p,pos,3); pos+=3;         /* IOD Corr */
            v3=getbits(p,pos,15); pos+=15;            /* radial */
            j=getbits(p,pos,13); pos+=13;             /* along */
            code=getbits(p,pos,13); pos+=13;          /* cross */
            mode=(uint32_t)getbitu(p,pos,6); pos+=6;  /* URAI */
            sat=slot2satno(slot);
            if (sat<=0||sat>=MAXSAT||v3==-16384||j==-4096||code==-4096) continue;
            /* B2bSSR storage in this project is indexed directly by RTKLIB
             * satellite number (index 0 is unused), unlike nav.eph[]. */
            ssr=&raw->nav.B2bssr[sat];
            ssr->sow=(int)sow;
            ssr->verify_sow=(int)(ep[3]*3600.0+ep[4]*60.0+ep[5]);
            ssr->t0[0]=epoch;
            ssr->udi[0]=udi;
            ssr->iodssr[0]=(int)iodssr;
            ssr->iodn=v1;
            ssr->iodcorr[0]=(uint16_t)v2;
            ssr->deph[0]=v3*0.0016;
            ssr->deph[1]=j*0.0064;
            ssr->deph[2]=code*0.0064;
            ssr->ura=(int)mode;
            ssr->update|=B2B_UPD_ORBIT|B2B_UPD_URA;
        }
        pos=body_end;                    /* 19 reserved bits */
        raw->num_PPPB2BINF02++;
        raw->raw_nmsg[1]++;
        output_B2bInfo2(raw,&raw->nav);
        mask->last_time[1]=epoch;
        return 2;
    }
    if (type==3) {                       /* code bias */
        udi=mask->last_time[2].time?timediff(epoch,mask->last_time[2]):0.0;
        if (udi<0.0||udi>86400.0) udi=0.0;
        nsat=getbitu(p,pos,5); pos+=5;
        for (i=0;i<(int)nsat;i++) {
            if (pos+13>body_end) return -1;
            slot=(int)getbitu(p,pos,9); pos+=9;
            nsig=getbitu(p,pos,4); pos+=4;
            sat=slot2satno(slot);
            codes=NULL;
            if (sat>0&&sat<MAXSAT) {
                sys=satsys(sat,NULL);
                if      (sys==SYS_CMP) codes=b2b_bds_codebias_mode;
                else if (sys==SYS_GPS) codes=b2b_gps_codebias_mode;
                else if (sys==SYS_GAL) codes=b2b_gal_codebias_mode;
                else if (sys==SYS_GLO) codes=b2b_glo_codebias_mode;
                ssr=&raw->nav.B2bssr[sat];
                memset(ssr->cbias,0,sizeof(ssr->cbias));
            }
            for (j=0;j<(int)nsig;j++) {
                if (pos+16>body_end) return -1;
                mode=getbitu(p,pos,4); pos+=4;
                v1=getbits(p,pos,12); pos+=12;
                if (!codes||mode>=B2B_CodeBiasModeNum||v1==-2048) continue;
                code=codes[mode];
                if (code>CODE_NONE&&code<=MAXCODE)
                    ssr->cbias[code-1]=(float)(v1*0.017);
            }
            if (codes) {
                ssr->sow=(int)sow;
                ssr->verify_sow=(int)(ep[3]*3600.0+ep[4]*60.0+ep[5]);
                ssr->t0[1]=epoch;
                ssr->udi[1]=udi;
                ssr->iodssr[1]=(int)iodssr;
                ssr->update|=B2B_UPD_CBIAS;
            }
        }
        raw->num_PPPB2BINF03++;
        raw->raw_nmsg[2]++;
        output_B2bInfo3(raw,&raw->nav);
        mask->last_time[2]=epoch;
        return 3;
    }

    iodp=getbitu(p,pos,4); pos+=4;       /* clock */
    mode=getbitu(p,pos,5); pos+=5;       /* subtype */
    if (mode>11) return 0;
    udi=mask->last_time[3].time?timediff(epoch,mask->last_time[3]):0.0;
    if (udi<0.0||udi>86400.0) udi=0.0;
    if (mask->IOD_SSR<0) {
        for (i=0;i<23;i++) {
            clock_iodcorr[i]=(int)getbitu(p,pos,3); pos+=3;
            clock_c0[i]=getbits(p,pos,15); pos+=15;
        }
        output_B2bClockIndex(raw,epoch,udi,(int)iodssr,(int)iodp,(int)mode,
                             clock_iodcorr,clock_c0);
        raw->num_PPPB2BINF04++;
        raw->raw_nmsg[3]++;
        mask->last_time[3]=epoch;
        return 4;
    }
    if ((int)iodp!=mask->IODP) return 0;
    for (i=0;i<23;i++) {
        v1=(int)getbitu(p,pos,3); pos+=3;
        v2=getbits(p,pos,15); pos+=15;
        j=(int)mode*23+i;
        /* The live stream uses the lower negative boundary as the explicit
         * unavailable value (repeated -26.2128 m records). */
        if (j>=mask->satnum||v2<=-16383) continue;
        sat=mask->satno[j];
        if (sat<=0||sat>=MAXSAT) continue;
        ssr=&raw->nav.B2bssr[sat];
        ssr->sow=(int)sow;
        ssr->verify_sow=(int)(ep[3]*3600.0+ep[4]*60.0+ep[5]);
        ssr->t0[2]=epoch;
        ssr->udi[2]=udi;
        ssr->iodssr[2]=(int)iodssr;
        ssr->iodp[0]=(int)iodp;
        ssr->iodcorr[1]=(uint16_t)v1;
        ssr->dclk[0]=v2*0.0016;
        ssr->dclk[1]=ssr->dclk[2]=0.0;
        ssr->update|=B2B_UPD_CLOCK;
    }
    raw->num_PPPB2BINF04++;
    raw->raw_nmsg[3]++;
    output_B2bInfo4(raw,&raw->nav);
    mask->last_time[3]=epoch;
    return 4;
}

/* merge decoded PPP-B2b corrections ---------------------------------------*/
extern int update_B2b(nav_t *nav, raw_t *raw)
{
    B2bssr_t *src, *dst;
    double dt;
    int i, j, n=0, changed;
    uint8_t upd;

    if (!nav||!raw) return -1;
    for (i=0;i<MAXSAT;i++) {
        src=&raw->nav.B2bssr[i];
        dst=&nav->B2bssr[i];
        upd=src->update;
        changed=0;
        if (!upd) continue;

        if (upd&B2B_UPD_ORBIT) {
            dt=dst->t0[0].time?timediff(src->t0[0],dst->t0[0]):0.0;
            src->udi[0]=(dt>0.0&&dt<=86400.0)?dt:0.0;
            dst->t0[0]=src->t0[0]; dst->udi[0]=src->udi[0];
            dst->iodssr[0]=src->iodssr[0]; dst->iodn=src->iodn;
            dst->iodcorr[0]=src->iodcorr[0]; dst->sow=src->sow;
            dst->verify_sow=src->verify_sow;
            for (j=0;j<3;j++) {dst->deph[j]=src->deph[j]; dst->ddeph[j]=src->ddeph[j];}
            changed=1;
        }
        if (upd&B2B_UPD_URA) {dst->ura=src->ura; changed=1;}
        if (upd&B2B_UPD_CBIAS) {
            dt=dst->t0[1].time?timediff(src->t0[1],dst->t0[1]):0.0;
            src->udi[1]=(dt>0.0&&dt<=86400.0)?dt:0.0;
            dst->t0[1]=src->t0[1]; dst->udi[1]=src->udi[1];
            dst->iodssr[1]=src->iodssr[1];
            memcpy(dst->cbias,src->cbias,sizeof(dst->cbias));
            changed=1;
        }
        if (upd&B2B_UPD_CLOCK) {
            dt=dst->t0[2].time?timediff(src->t0[2],dst->t0[2]):0.0;
            src->udi[2]=(dt>0.0&&dt<=86400.0)?dt:0.0;
            dst->t0[2]=src->t0[2]; dst->udi[2]=src->udi[2];
            dst->iodssr[2]=src->iodssr[2]; dst->iodp[0]=src->iodp[0];
            dst->iodcorr[1]=src->iodcorr[1];
            for (j=0;j<3;j++) dst->dclk[j]=src->dclk[j];
            changed=1;
        }
        src->update=0;
        dst->update=0;
        if (changed) n++;
    }
    raw->num_PPPB2BINF01=raw->num_PPPB2BINF02=0;
    raw->num_PPPB2BINF03=raw->num_PPPB2BINF04=0;
    return n;
}

static int adjday_B2b(double H_sod, int H_doy, double D_sod){
    double tt= H_sod-D_sod;
    int D_doy;
    if (tt > 43200){
        printf("!!!D_doy = H_doy + 1!!!");
        D_doy = H_doy + 1;
        return D_doy;
    }
    if (tt < -43200){
        D_doy = H_doy - 1;
        return D_doy;
    }

    return H_doy;
}

static void adjyear_B2b(int H_year ,int *D_year, int *D_doy){

    int maxday;
    if((H_year % 4 == 0 && H_year % 100 != 0) || H_year % 400 == 0){
        maxday = 366;
    }
    else maxday = 365;

    if (*D_doy > maxday){
        printf("!!!D_year = H_year + 1 and D_doy = 1!!!");
        *D_year = H_year+1;
        *D_doy = 1;
    }
    else *D_year = H_year;
}



// Convert B2b time data: BDS: year + doy (provided by header) + tod (carried by data) --> BDS: mjd + tod --> GPS: mjd + tod --> rtklib_time (gpst)
extern gtime_t B2btod2time(gtime_t H_tt, double D_sod)
{
    double H_ep[6], D_sow;
    int H_year, H_doy, H_sod, D_year, D_doy, D_week;
    mjd_gtime_t mjdsod = {0};
    gtime_t H_bdt = {0};
    gtime_t D_time = {0};

    // Convert header time from GPST to BDST to determine the year and day of year in BeiDou time
    H_bdt = gpst2bdt(H_tt);
    time2epoch(H_bdt, H_ep);
    H_year = (int)H_ep[0];
    H_sod = (int)H_ep[3] * 3600 + (int)H_ep[4] * 60 + (int)H_ep[5];
    H_doy = (int)time2doy(H_bdt);

    // Calculate the year and day of year for the data (in BDST) based on the header's year and day of year
    D_doy = adjday_B2b(H_sod, H_doy, D_sod);
    adjyear_B2b(H_year, &D_year, &D_doy);

    // Convert the data's year and day of year to MJD, then to GPST
    mjdsod = yrdoy_to_mjd_time(D_year, D_doy, D_sod);
    mjdsod = bdst_to_gpst(mjdsod);
    mjd_time_to_gpst(mjdsod, &D_week, &D_sow);

    // Convert the data's GPST to RTKLIB's time structure
    D_time = gpst2time(D_week, D_sow);
    
    return D_time;
} /* end of time2bdst */

extern int checkout_B2beph(B2bssr_t *B2bssr0, B2bssr_t *B2bssr1) {

    if (B2bssr0->iodssr[0] != B2bssr1->iodssr[0]) {
        printf("Error: iodssr mismatch! B2bssr0->iodssr = %d, B2bssr1->iodssr = %d\n", 
               B2bssr0->iodssr[0], B2bssr1->iodssr[0]);
        return 0;
    }

    if (B2bssr0->iodn != B2bssr1->iodn) {
        printf("Error: iodn mismatch! B2bssr0->iodn = %d, B2bssr1->iodn = %d\n", 
               B2bssr0->iodn, B2bssr1->iodn);
        return 0;
    }

    if (B2bssr0->iodcorr[0] != B2bssr1->iodcorr[0]) {
        printf("Error: iodcorr mismatch! B2bssr0->iodcorr = %d, B2bssr1->iodcorr = %d\n", 
               B2bssr0->iodcorr[0], B2bssr1->iodcorr[0]);
        return 0;
    }

    if (B2bssr0->ura != B2bssr1->ura) {
        printf("Error: ura mismatch! B2bssr0->ura = %d, B2bssr1->ura = %d\n", 
               B2bssr0->ura, B2bssr1->ura);
        return 0;
    }

    if (fabs(B2bssr0->deph[0] - B2bssr1->deph[0]) > EPSILON || 
        fabs(B2bssr0->deph[1] - B2bssr1->deph[1]) > EPSILON || 
        fabs(B2bssr0->deph[2] - B2bssr1->deph[2]) > EPSILON) {
        printf("Error: deph mismatch! B2bssr0->deph = [%f, %f, %f], B2bssr1->deph = [%f, %f, %f]\n", 
               B2bssr0->deph[0], B2bssr0->deph[1], B2bssr0->deph[2],
               B2bssr1->deph[0], B2bssr1->deph[1], B2bssr1->deph[2]);
        return 0;
    }

    return 1;
}


extern int checkout_B2bcbia(B2bssr_t *B2bssr0, B2bssr_t *B2bssr1) {

    if (B2bssr0->iodssr[1] != B2bssr1->iodssr[1]) {
        printf("Error: iodssr mismatch! B2bssr0->iodssr = %d, B2bssr1->iodssr = %d\n", 
               B2bssr0->iodssr[1], B2bssr1->iodssr[1]);
        return 0;
    }

    for (int i = 0; i < MAXCODE; i++) {
        if (fabs(B2bssr0->cbias[i] - B2bssr1->cbias[i]) > EPSILON) {
            printf("Error: cbias mismatch! B2bssr0->cbias[%d] = %f, B2bssr1->cbias[%d] = %f\n", 
                   i, B2bssr0->cbias[i], i, B2bssr1->cbias[i]);
            return 0;
        }
    }

    return 1;
}


extern int checkout_B2bclk(B2bssr_t *B2bssr0, B2bssr_t *B2bssr1) {

    if (B2bssr0->iodssr[2] != B2bssr1->iodssr[2]) {
        printf("Error: iodssr mismatch! B2bssr0->iodssr = %d, B2bssr1->iodssr = %d\n", 
               B2bssr0->iodssr[2], B2bssr1->iodssr[2]);
        return 0;
    }
    if (B2bssr0->iodcorr[1] != B2bssr1->iodcorr[1]) {
        printf("Error: iodssr mismatch! B2bssr0->iodcorr = %d, B2bssr1->iodcorr = %d\n", 
               B2bssr0->iodcorr[1], B2bssr1->iodcorr[1]);
        return 0;
    }
    if (fabs(B2bssr0->dclk[0] - B2bssr1->dclk[0]) > EPSILON) {
        printf("Error: iodssr mismatch! B2bssr0->dclk = %f, B2bssr1->dclk = %f\n", 
               B2bssr0->dclk[0], B2bssr1->dclk[0]);
        return 0;
    }

    return 1;
}

/* Convert URAI to URA */
extern void B2b_urai2ura(gtime_t time, const char *satid, int URAI, FILE *UraFile) {
    char time_str[64];

    if (UraFile == NULL) {
        fprintf(stderr, "Error: File pointer is NULL.\n");
        return;
    }

    time2str(time, time_str, 3);

    /* If URAI is undefined or out of range */
    if (URAI == 0) {
        fprintf(UraFile, "[%s] Satellite: %s - URAI is undefined or unknown. Corrections are unreliable.\n", 
                time_str, satid);
    } 
    else if (URAI == 63) { /* 63 == 111111 binary */
        fprintf(UraFile, "[%s] Satellite: %s - URAI indicates URA > 5466.5 mm. Corrections are unreliable.\n", 
                time_str, satid);
    } 
    else {
        /* Parse URA_CLASS and URA_VALUE */
        int URA_CLASS = (URAI >> 3) & 0x07;  /* Top 3 bits */
        int URA_VALUE = URAI & 0x07;         /* Bottom 3 bits */

        /* Calculate URA using the formula */
        double URA = pow(3.0, URA_CLASS) * (1.0 + 0.25 * URA_VALUE) - 1.0;

        fprintf(UraFile, "[%s] Satellite: %s - URA_CLASS: %d, URA_VALUE: %d, URA: %.3f mm\n", 
                time_str, satid, URA_CLASS, URA_VALUE, URA);
    }

    fflush(UraFile);
}

// Convert mask array to ASCII string
static void MaskToStr(const int *mask, int len, char *str) {
    for (int i = 0; i < len; i++) {
        str[i] = mask[i] ? '1' : '0';
    }
    str[len] = '\0';
}

static void SatnoToPrn(const B2bmask_t *mask, char *bds_prn, char *gps_prn, char *gal_prn, char *glo_prn) {
    char prn_str[8];
    int sys, prn;

    bds_prn[0] = gps_prn[0] = gal_prn[0] = glo_prn[0] = '\0';


    for (int i = 0; i < mask->satnum; i++) {
        int sat = mask->satno[i];
        sys = satsys(sat, &prn); 
        switch (sys) {
            case SYS_CMP:  
                sprintf(prn_str, "C%02d", prn); 
                strcat(bds_prn, prn_str);
                strcat(bds_prn, " ");
                break;
            case SYS_GPS:  
                sprintf(prn_str, "G%02d", prn); 
                strcat(gps_prn, prn_str);
                strcat(gps_prn, " ");
                break;
            case SYS_GAL:  
                sprintf(prn_str, "E%02d", prn);
                strcat(gal_prn, prn_str);
                strcat(gal_prn, " ");
                break;
            case SYS_GLO:
                sprintf(prn_str, "R%02d", prn);  
                strcat(glo_prn, prn_str);
                strcat(glo_prn, " ");
                break;
            default:
                break;
        }
    }
}

// Output B2b mask information
extern void output_B2bInfo1(raw_t *raw, const B2bmask_t *mask, int udi) {
    double m_epoch[6], t_epoch[6];
    char str_MASK_BD[64], str_MASK_GPS[38], str_MASK_GALILEO[38], str_MASK_GLONASS[38];
    char str_BDSPRN[512], str_GPSPRN[512], str_GALPRN[512], str_GLOPRN[512];

    // Generate mask strings
    MaskToStr(mask->MASK_BD, 63, str_MASK_BD);
    MaskToStr(mask->MASK_GPS, 37, str_MASK_GPS);
    MaskToStr(mask->MASK_GALILEO, 37, str_MASK_GALILEO);
    MaskToStr(mask->MASK_GLONASS, 37, str_MASK_GLONASS);

    // Generate PRN lists
    SatnoToPrn(mask, str_BDSPRN, str_GPSPRN, str_GALPRN, str_GLOPRN);

    // Get times: reference time and receive time
    time2epoch(mask->time, t_epoch);
    time2epoch(raw->time, m_epoch);

    // Output epoch record
    B2b_trace(22, "> MASK %04d %02d %02d %02d %02d %.1f %d %d B2bSatPrn_%d @ %04d %02d %02d %02d %02d %.1f\n",
              (int)t_epoch[0], (int)t_epoch[1], (int)t_epoch[2], (int)t_epoch[3],
              (int)t_epoch[4], t_epoch[5], udi, mask->satnum, raw->geoprn,
              (int)m_epoch[0], (int)m_epoch[1], (int)m_epoch[2], (int)m_epoch[3],
              (int)m_epoch[4], m_epoch[5]);

    // Output IOD information
    B2b_trace(22, "IOD:      %d            %d\n", mask->IOD_SSR, mask->IODP);

    // Output mask information
    B2b_trace(22, "BDSMASK:  %s\n", str_MASK_BD);
    B2b_trace(22, "GPSMASK:  %s\n", str_MASK_GPS);
    B2b_trace(22, "GALMASK:  %s\n", str_MASK_GALILEO);
    B2b_trace(22, "GLOMASK:  %s\n", str_MASK_GLONASS);

    // Output PRN lists
    B2b_trace(22, "BDS PRN:  %s\n", str_BDSPRN);
    B2b_trace(22, "GPS PRN:  %s\n", str_GPSPRN);
    B2b_trace(22, "GAL PRN:  %s\n", str_GALPRN);
    B2b_trace(22, "GLO PRN:  %s\n", str_GLOPRN);
}

// Output B2b orbit and URAI information
extern void output_B2bInfo2(raw_t *raw, nav_t *nav) {
    char satid[16];
    double t_epoch[6] = {0}, m_epoch[6];
    int udi = 0, satnum = 0;
    gtime_t t_time;

    // Get current time
    time2epoch(raw->time, m_epoch);

    // Get udi and t0 of the first valid satellite (also update satnum)
    for (int i = 1; i < MAXSAT; i++) {
        const B2bssr_t *b2b = &nav->B2bssr[i];
        if (b2b->update == 0) continue;
        satnum++;
        udi = (int)b2b->udi[0];  // Update interval for orbit correction
        t_time = b2b->t0[0];  // Reference time for orbit correction
        time2epoch(t_time, t_epoch);
    }
    if (satnum == 0) return;

    // Output epoch record
    B2b_trace(22, "> ORBIT_URAI %04d %02d %02d %02d %02d %.1f %d %d B2bSatPrn_%d @ %04d %02d %02d %02d %02d %.1f\n",
              (int)t_epoch[0], (int)t_epoch[1], (int)t_epoch[2], (int)t_epoch[3], 
              (int)t_epoch[4], t_epoch[5], udi, satnum, raw->geoprn,
              (int)m_epoch[0], (int)m_epoch[1], (int)m_epoch[2], (int)m_epoch[3], 
              (int)m_epoch[4], m_epoch[5]);

    // Traverse all satellites
    for (int i = 1; i < MAXSAT; i++) {
        const B2bssr_t *b2b = &nav->B2bssr[i];
        if (b2b->update == 0) continue;  // Skip satellites without updates

        satno2id(i, satid);  // Convert satellite number to PRN, e.g., G01, C01

        // Satellite record (URAI not provided, set to 0 by default)
        B2b_trace(22, "%-4s %5d %5d %5d %5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                  satid, b2b->iodssr[0], b2b->iodn, b2b->iodcorr[0], b2b->ura,  // PRN, IODN, IOD Corr, URAI
                  b2b->deph[0], b2b->deph[1], b2b->deph[2],
                  b2b->ddeph[0] * 1000.0, b2b->ddeph[1] * 1000.0, b2b->ddeph[2] * 1000.0);  // Velocity (mm/s)
    }
}

// Output B2b differential code bias information
extern void output_B2bInfo3(raw_t *raw, nav_t *nav) {
    char satid[16];
    double t_epoch[6] = {0}, m_epoch[6];
    int udi = 0, satnum = 0;
    gtime_t t_time;

    // Get current time
    time2epoch(raw->time, m_epoch);

    // Get udi and t0 of the first valid satellite (also update satnum)
    for (int i = 1; i < MAXSAT; i++) {
        const B2bssr_t *b2b = &nav->B2bssr[i];
        if (b2b->update == 0) continue;
        satnum++;
        udi = (int)b2b->udi[1];  // Update interval for code bias
        t_time = b2b->t0[1];  // Reference time for code bias
        time2epoch(t_time, t_epoch);
    }
    if (satnum == 0) return;

    // Output epoch record
    B2b_trace(22, "> DIFF_CODE_BIAS %04d %02d %02d %02d %02d %.1f %d %d B2bSatPrn_%d @ %04d %02d %02d %02d %02d %.1f\n",
              (int)t_epoch[0], (int)t_epoch[1], (int)t_epoch[2], (int)t_epoch[3], 
              (int)t_epoch[4], t_epoch[5], udi, satnum, raw->geoprn,
              (int)m_epoch[0], (int)m_epoch[1], (int)m_epoch[2], (int)m_epoch[3], 
              (int)m_epoch[4], m_epoch[5]);

    // Traverse all satellites
    for (int i = 1; i < MAXSAT; i++) {
        const B2bssr_t *b2b = &nav->B2bssr[i];
        if (b2b->update == 0) continue;

        satno2id(i, satid);  // Get PRN, e.g., G01, C01

        // Calculate the number of code biases
        int bias_num = 0;
        for (int j = 0; j < MAXCODE; j++) {
            if (b2b->cbias[j] != 0.0) bias_num++;
        }

        // Output satellite PRN and number of biases
        B2b_trace(22, "%-4s %5d %5d", satid, b2b->iodssr[1], bias_num);

        // Output signal indicators and bias values
        for (int j = 0; j < MAXCODE; j++) {
            if (b2b->cbias[j] != 0.0) {
                char signal[4] = {0};
                strcpy(signal, rtklib_codes[j]);
                // rtklib_codes[j];
                // code_to_signal(j, signal);  // Convert code type to signal string
                B2b_trace(22, " %7s %7.4f", signal, b2b->cbias[j]);
            }
        }
        B2b_trace(22, "\n");
    }
}

// Output B2b clock correction information
extern void output_B2bInfo4(raw_t *raw, nav_t *nav) {
    char satid[16];
    double t_epoch[6] = {0}, m_epoch[6];
    int udi = 0, satnum = 0;
    gtime_t t_time;

    // Get current time
    time2epoch(raw->time, m_epoch);

    // Get udi and t0 of the first valid satellite (clock uses index 2)
    for (int i = 1; i < MAXSAT; i++) {
        const B2bssr_t *b2b = &nav->B2bssr[i];
        if (b2b->update == 0) continue;
        satnum++;
        udi = (int)b2b->udi[2];  // Update interval for clock correction
        t_time = b2b->t0[2];  // Reference time for clock correction
        time2epoch(t_time, t_epoch);
    }
    if (satnum == 0) return;

    // Output epoch record
    B2b_trace(22, "> CLOCK %04d %02d %02d %02d %02d %.1f %d %d B2bSatPrn_%d @ %04d %02d %02d %02d %02d %.1f\n",
              (int)t_epoch[0], (int)t_epoch[1], (int)t_epoch[2], (int)t_epoch[3], 
              (int)t_epoch[4], t_epoch[5], udi, satnum, raw->geoprn,
              (int)m_epoch[0], (int)m_epoch[1], (int)m_epoch[2], (int)m_epoch[3], 
              (int)m_epoch[4], m_epoch[5]);

    // Output satellite records
    for (int i = 1; i < MAXSAT; i++) {
        const B2bssr_t *b2b = &nav->B2bssr[i];
        if (b2b->update == 0) continue;

        satno2id(i, satid);

        // Output PRN, IOD_SSR, IODP, and C0, C1, C2
        B2b_trace(22, "%-4s %5d %5d %10.4f %10.4f %10.4f\n",
                  satid, 
                  b2b->iodssr[2],  // IOD_SSR
                  b2b->iodp[0],    // IODP
                  b2b->dclk[0],    // C0 (meters)
                  b2b->dclk[1] * 1000.0,  // C1 (mm/s)
                  b2b->dclk[2] * 1000.0); // C2 (mm/s²)
    }
}
