/* Test-only golden implementation from Git 7e77af5, before unification.
 * Not linked into production. Function names changed only to avoid collisions. */
/* antenna corrected measurements --------------------------------------------*/
static void reference_legacy_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *dantr, const double *rs,const double *e,double elapsedtime,
                      const double *dants, double phw, double *L, double *P,
                      double *Lc, double *Pc)
{
    double freq[NFREQ]={0},C1,C2,cbias[2]={0};
    int i,ix=0,frq,frq2,bias_ix,sys=satsys(obs->sat,NULL);

    double	dantf[NFREQ][3] = { 0x00 };
    *Lc=*Pc=0.0;
    for (i=0;i<NFREQ;i++) L[i]=P[i]=0.0;
    if (sys==SYS_CMP&&opt->bds_if) {
        for (i=0;i<2;i++) {
            if (bds_ant_index(nav->pcvs+obs->sat-1,obs->code[i],0)<0) return;
            if (opt->sateph==EPHOPT_B2b &&
                !bds_code_bias(obs->time,&nav->B2bssr[obs->sat],obs->code[i],cbias+i)) {
                trace(3,"bds_if: missing/stale/inconsistent DCB sat=%d code=%d\n",obs->sat,obs->code[i]);
                return;
            }
        }
    }

	if (rs != 0 && (opt->sateph == EPHOPT_PREC || opt->sateph == EPHOPT_SSRCOM || opt->sateph == EPHOPT_B2b))
	{
		double	rawdantf[NFREQ][3] = { 0x00 };
		satantoff1(obs->time, rs, obs->sat, nav, obs->code, rawdantf[0], rawdantf[1]);

		dantf[0][0] = rawdantf[0][0] + rawdantf[0][1] * OMGE * elapsedtime;
		dantf[0][1] = rawdantf[0][1] - rawdantf[0][0] * OMGE * elapsedtime;
		dantf[0][2] = rawdantf[0][2];

		dantf[1][0] = rawdantf[1][0] + rawdantf[1][1] * OMGE * elapsedtime;
		dantf[1][1] = rawdantf[1][1] - rawdantf[1][0] * OMGE * elapsedtime;
		dantf[1][2] = rawdantf[1][2];
	}

    for (i=0;i<opt->nf;i++) {
        L[i]=P[i]=0.0;
        /* skip if low SNR or missing observations */
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);
        if (freq[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) 
            continue;

        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*CLIGHT/freq[i]-dants[i]-dantr[i]-phw*CLIGHT/freq[i];
        P[i]=obs->P[i]               -dants[i]-dantr[i];

		// Correct PCO
		if (rs != 0)
		{
			L[i] = L[i] - (dantf[i][0] * e[0] + dantf[i][1] * e[1] + dantf[i][2] * e[2]);
			P[i] = P[i] - (dantf[i][0] * e[0] + dantf[i][1] * e[1] + dantf[i][2] * e[2]);
		}

        if (opt->sateph==EPHOPT_SSRAPC||opt->sateph==EPHOPT_SSRCOM) {
            /* select SSR code correction based on code */
			if (sys==SYS_GPS)
                ix=(i==0?CODE_L1W-1:CODE_L2W-1);
            else if (sys==SYS_GLO)
                ix=(i==0?CODE_L1P-1:CODE_L2P-1);
            else if (sys==SYS_GAL)
				ix=(i==0?CODE_L1X-1:CODE_L7X-1);
			else if (sys==SYS_CMP)
				ix=(i==0?CODE_L2I-1:CODE_L6I-1);

			/* apply SSR correction */
//			P[i]+=(nav->ssr[obs->sat-1].cbias[obs->code[i]-1]-nav->ssr[obs->sat-1].cbias[ix]);

			P[i]+=(nav->ssr[obs->sat-1].cbias[obs->code[i]-1]);

/*			if(nav->ssr[obs->sat-1].pbias[obs->code[i]-1] != 0)
			{
				L[i]+=nav->ssr[obs->sat-1].pbias[obs->code[i]-1];
            }*/
        }
        if (opt->sateph==EPHOPT_B2b) {
            /* select SSR code correction based on code */
			if (sys==SYS_GPS)
                ix=(i==0?CODE_L1W:CODE_L2W);
            else if (sys==SYS_GLO)
                ix=(i==0?CODE_L1P:CODE_L2P);
            else if (sys==SYS_GAL)
				ix=(i==0?CODE_L1X:CODE_L7X);
			else if (sys==SYS_CMP)
				ix=obs->code[i];

            P[i]-=(nav->B2bssr[obs->sat].cbias[ix]);
        }
        else {   /* apply code bias corrections from file */
            if (sys==SYS_GAL&&(i==1||i==2)) frq=3-i;  /* GAL biases are L1/L5 */
            else frq=i;  /* other biases are L1/L2 */
            if (frq>=MAX_CODE_BIAS_FREQS) continue;  /* only 2 freqs per system supported in code bias table */
            bias_ix=code2bias_ix(sys,obs->code[i]); /* look up bias index in table */
            if (bias_ix>0) {  /*  0=ref code */
                P[i]+=nav->cbias[obs->sat-1][frq][bias_ix-1]; /* code bias */
            }
        }
    }
    /* choose freqs for iono-free LC */
    *Lc=*Pc=0.0;
    frq2=(sys==SYS_CMP&&opt->bds_if)?1:(L[1]==0?2:1);
    if (freq[0]==0.0||freq[frq2]==0.0) return;
    if (fabs(freq[0]-freq[frq2])<1.0) return;
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[frq2]));
    C2=-SQR(freq[frq2])/(SQR(freq[0])-SQR(freq[frq2]));

    if (L[0]!=0.0&&L[frq2]!=0.0) *Lc=C1*L[0]+C2*L[frq2];
    if (P[0]!=0.0&&P[frq2]!=0.0) *Pc=C1*P[0]+C2*P[frq2];
    if (sys==SYS_CMP&&opt->bds_if) {
        trace(4,"bds_iflc: sat=%d code=%d/%d freq=%.3f/%.3f coef=%.9f/%.9f dcb=%.3f/%.3f Pc=%.3f\n",
              obs->sat,obs->code[0],obs->code[1],freq[0],freq[1],C1,C2,cbias[0],cbias[1],*Pc);
    }
}
typedef struct {
    double L[2],P[2],a[2][3],cov[2][2][2]; /* covariance [phase/code][pair][pair] */
} reference_ifobs_t;

/* GPS: TF-C link IFCB is estimated, not a zero-valued L5 DCB correction.
 * BDS: three actual code biases are tied to the B2b B3I clock datum first.
 * Model reference: Wang et al., Remote Sensing 14 (2022) 4509, TF-C/TF-F.
 * The small time-varying phase terms left in TF-C code are code-level noise. */
static void reference_configured_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                       const prcopt_t *opt, const double *rs, const double *rr,
                       const double *e, double elapsed, double phw, reference_ifobs_t *out)
{
    double L[3]={0},P[3]={0},freq[3]={0},rv[2][3]={{0}},bias;
    double off[3],rot[3],dr,ds,ru[3],rz[3],eu[3],ez[3],nadir=0.0,cosa;
    int i,j,k,c,sys=satsys(obs->sat,NULL),valid[3]={0},wi;
    uint8_t requested[3]={CODE_L1C,CODE_L2W,CODE_L5Q};
    prcopt_t noiseopt=*opt;
    memset(out,0,sizeof(*out));
    if (sys==SYS_CMP&&!bds_parse_freqs(opt->bds_freqs,requested)) return;
    if (rs) {
        for (j=0;j<3;j++) {ru[j]=rr[j]-rs[j];rz[j]=-rs[j];}
        if (!normv3(ru,eu)||!normv3(rz,ez)) return;
        cosa=dot(eu,ez,3); nadir=acos(MAX(-1.0,MIN(1.0,cosa)));
    }
    for (i=0;i<3;i++) {
        const char *s=code2obs(obs->code[i]);
        if (s[0]!=code2obs(requested[i])[0]) continue;
        if (sys==SYS_CMP&&(s[0]=='2'||s[0]=='6'?s[1]!='I':s[1]!='P'&&s[1]!='D')) continue;
        if ((wi=signal_noise_index(sys,obs->code[i]))<0) continue;
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);
        if (freq[i]<=0||!obs->P[i]||!obs->L[i]||
            testsnr(0,wi,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;
        if (signal_ant_index(nav->pcvs+obs->sat-1,sys,obs->code[i],0)<0||
            signal_ant_index(opt->pcvr,sys,obs->code[i],opt->bds_ant_fallback)<0) continue;
        bias=0.0;
        if (sys==SYS_CMP&&!bds_code_bias(obs->time,nav->B2bssr+obs->sat,obs->code[i],&bias)) {
            trace(3,"IF1213 missing BDS DCB sat=%d code=%s\n",obs->sat,code2obs(obs->code[i]));
            continue;
        }
        L[i]=obs->L[i]*CLIGHT/freq[i]-phw*CLIGHT/freq[i];
        P[i]=obs->P[i]-bias;
        if (rs) {
            if (!signal_antmodel(opt->pcvr,sys,obs->code[i],opt->antdel[0],azel,0,
                                 opt->posopt[1],opt->bds_ant_fallback,&dr)||
                !signal_antmodel(nav->pcvs+obs->sat-1,sys,obs->code[i],NULL,NULL,nadir,
                                 opt->posopt[0],0,&ds)||
                !signal_satantoff(obs->time,rs,obs->sat,nav,obs->code[i],off)) continue;
            rot[0]=off[0]+off[1]*OMGE*elapsed;
            rot[1]=off[1]-off[0]*OMGE*elapsed;rot[2]=off[2];
            L[i]-=dr+ds+dot(rot,e,3); P[i]-=dr+ds+dot(rot,e,3);
        }
        valid[i]=1;
        noiseopt.eratio[i]=opt->eratio[wi];
        for (c=0;c<2;c++) rv[c][i]=varerr(obs->sat,sys,azel[1],
            obs->SNR[i]*SNR_UNIT,i*2+c,&noiseopt,obs);
    }
    for (k=0;k<NF(opt);k++) {
        j=k+1;
        if (!valid[0]||!valid[j]||!if_coefficients(freq[0],freq[j],out->a[k],out->a[k]+j)) continue;
        out->L[k]=out->a[k][0]*L[0]+out->a[k][j]*L[j];
        out->P[k]=out->a[k][0]*P[0]+out->a[k][j]*P[j];
    }
    for (c=0;c<2;c++) for (i=0;i<NF(opt);i++) for (j=0;j<NF(opt);j++)
        out->cov[c][i][j]=if_covariance(out->a[i],out->a[j],rv[c]);
}
