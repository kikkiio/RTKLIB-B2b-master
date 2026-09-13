/* Configurable BeiDou observation slots. Physical code2freq/code2idx and
 * PPP-B2b's B3 reference antenna are deliberately left unchanged. */
#include "rtklib.h"
#include <ctype.h>

/* Public names, not RINEX frequency digits. No implicit P/D/X equivalence. */
extern int bds_parse_freqs(const char *text, uint8_t codes[3])
{
    const char *p=text?text:"";
    int i,j,n;char word[8];
    while (isspace((unsigned char)*p)) p++;
    if (!*p) {
        codes[0]=CODE_L5P;codes[1]=CODE_L1P;codes[2]=CODE_L2I;return 1;
    }
    for (i=0;i<3;i++) {
        n=0;while (isspace((unsigned char)*p)) p++;
        while (*p&&*p!=','&&!isspace((unsigned char)*p)) {
            if (n>=7) return 0;
            word[n++]=(char)toupper((unsigned char)*p++);
        }
        word[n]=0;while (isspace((unsigned char)*p)) p++;
        if      (!strcmp(word,"B1I")) codes[i]=CODE_L2I;
        else if (!strcmp(word,"B3I")) codes[i]=CODE_L6I;
        else if (!strcmp(word,"B1C")) codes[i]=CODE_L1P;
        else if (!strcmp(word,"B2A")) codes[i]=CODE_L5P;
        else return 0;
        for (j=0;j<i;j++) if (codes[j]==codes[i]) return 0;
        if (i==1&&!*p) {codes[2]=0;return 1;}
        if (i<2) {if (*p++!=',') return 0;}
        else if (*p) return 0;
    }
    return 1;
}

/* Pair syntax is the public interface; old numeric switches are only a
 * compatibility representation used by existing state/decoder code. */
extern int if_parse_pairs(int sys, const char *text, uint8_t codes[3])
{
    char compact[64],left[8],right[8],anchor[8]="",bands[32]="";
    const char *p=text?text:"";char *q,*end;int n=0,k=0,i;
    memset(codes,0,3*sizeof(*codes));
    while (*p) {
        if (!isspace((unsigned char)*p)) {
            if (n>=63) return 0;
            compact[n++]=(char)toupper((unsigned char)*p);
        }
        p++;
    }
    compact[n]=0;
    if (sys==SYS_GAL) {
        if (strcmp(compact,"E1/E5A")&&strcmp(compact,"E1/E5A,E1/E5B")) return 0;
        codes[0]=CODE_L1C;codes[1]=CODE_L5Q;
        if (strchr(compact,',')) codes[2]=CODE_L7Q;
        return codes[2]?2:1;
    }
    if (sys==SYS_GPS) {
        if (strcmp(compact,"L1/L2")&&strcmp(compact,"L1/L2,L1/L5")) return 0;
        codes[0]=CODE_L1C;codes[1]=CODE_L2W;
        if (strchr(compact,',')) codes[2]=CODE_L5Q;
        return codes[2]?2:1;
    }
    if (sys!=SYS_CMP) return 0;
    q=compact;
    while (*q) {
        if (k==2) return 0;
        end=strchr(q,',');if (end) *end=0;
        p=strchr(q,'/');if (!p||strchr(p+1,'/')) return 0;
        i=(int)(p-q);if (i<1||i>=8||!p[1]||strlen(p+1)>=8) return 0;
        memcpy(left,q,i);left[i]=0;strcpy(right,p+1);
        if (!k) {strcpy(anchor,left);snprintf(bands,sizeof(bands),"%s,%s",left,right);}
        else {if (strcmp(left,anchor)) return 0;strcat(bands,",");strcat(bands,right);}
        k++;
        if (!end) break;
        q=end+1;if (!*q) return 0;
    }
    return k&&bds_parse_freqs(bands,codes)?k:0;
}

extern int if_options_normalize(prcopt_t *opt, char *msg)
{
    uint8_t gps[3],bds[3],gal[3];int ng=0,nb=0,ne=0,n;
    const char *names[4]={"B1I","B3I","B1C","B2a"};
    const uint8_t values[4]={CODE_L2I,CODE_L6I,CODE_L1P,CODE_L5P};
    int i,j;
    if (opt->robust<0||opt->robust>1||!isfinite(opt->robust_k)||
        (opt->robust_k!=0.0&&(opt->robust_k<1.0||opt->robust_k>10.0))) {
        strcpy(msg,"PPP robust must be 0/1; robust_k must be 0(default) or 1..10 sigma");return 0;
    }
    if (!opt->gps_if_pairs[0]&&!opt->bds_if_pairs[0]&&!opt->gal_if_pairs[0]) return 1;
    if (opt->if_model<0) {strcpy(msg,"invalid configuration (including replay_end)");return 0;}
    if (opt->navsys&~(SYS_GPS|SYS_CMP|SYS_GAL)) {
        strcpy(msg,"configured IF supports GPS/BDS/GAL only");return 0;
    }
    if (opt->gps_if_pairs[0]) ng=if_parse_pairs(SYS_GPS,opt->gps_if_pairs,gps);
    if (opt->bds_if_pairs[0]) nb=if_parse_pairs(SYS_CMP,opt->bds_if_pairs,bds);
    if (opt->gal_if_pairs[0]) ne=if_parse_pairs(SYS_GAL,opt->gal_if_pairs,gal);
    if ((opt->gps_if_pairs[0]&&!ng)||(opt->bds_if_pairs[0]&&!nb)||
        ((opt->navsys&SYS_GPS)&&!ng)||((opt->navsys&SYS_CMP)&&!nb)||
        (opt->gal_if_pairs[0]&&!ne)||((opt->navsys&SYS_GAL)&&!ne)||
        (ng&&nb&&ng!=nb)||(ng&&ne&&ng!=ne)||(nb&&ne&&nb!=ne)||!(n=ng?ng:nb?nb:ne)) {
        strcpy(msg,"IF pairs: GPS L1/L2[,L1/L5]; GAL E1/E5a[,E1/E5b]; BDS common-anchor; same count required");return 0;
    }
    opt->nf=n+1;opt->if_model=n==2?1:2;opt->bds_if=2;
    opt->bds_freqs[0]=0;
    if (!nb) strcpy(opt->bds_freqs,n==2?"B2a,B1C,B1I":"B1I,B3I");
    else for (i=0;i<n+1;i++) {
        for (j=0;j<4&&values[j]!=bds[i];j++);
        if (i) strcat(opt->bds_freqs,",");
        strcat(opt->bds_freqs,names[j]);
    }
    return 1;
}

static void bds_slot_map(const uint8_t codes[3], int map[6])
{
    int i,j=codes[2]?3:2;
    for (i=0;i<6;i++) map[i]=-1;
    for (i=0;i<3&&codes[i];i++) map[code2idx(SYS_CMP,codes[i])]=i;
    for (i=0;i<6;i++) if (map[i]<0) map[i]=j++;
}

extern int bds_options_valid(const prcopt_t *opt, char *msg)
{
    uint8_t codes[3];
    if (opt->gal_ephemeris<0||opt->gal_ephemeris>1||!isfinite(opt->gal_brdc_sigma)||opt->gal_brdc_sigma<0) {
        strcpy(msg,"invalid GAL product or broadcast sigma option");return 0;
    }
    if (opt->if_model&&(opt->navsys&SYS_GAL)&&!(opt->navsys&(SYS_GPS|SYS_CMP))) {
        strcpy(msg,"GAL broadcast assistance currently requires GPS/BDS for SPP initialization");return 0;
    }
    if (opt->if_model&&(opt->navsys&SYS_GAL)&&
        (opt->gal_ephemeris!=1||!opt->gal_navfile[0]||
         if_parse_pairs(SYS_GAL,opt->gal_if_pairs,codes)!=(opt->if_model==1?2:1))) {
        strcpy(msg,"GAL IF requires explicit gal_ephemeris=1, gal_navfile and E1/E5a[,E1/E5b]");return 0;
    }
    if (!bds_parse_freqs(opt->bds_freqs,codes)||
        (opt->bds_freqs[0]&&(!opt->if_model||opt->bds_if!=2))) {
        strcpy(msg,"bds_freqs requires configured IF,bds_if=2 and distinct B1I/B3I/B1C/B2a bands");
        return 0;
    }
    if (opt->bds_if<0||opt->bds_if>2||opt->bds_ant_fallback<0||opt->bds_ant_fallback>1||
        opt->if_model<0||opt->if_model>2||!isfinite(opt->ifcb_prn)||opt->ifcb_prn<0) {
        strcpy(msg,"invalid BDS/IF model or IFCB process noise option");
        return 0;
    }
    if (opt->if_model) {
        if (NFREQ<3||(opt->if_model==1?opt->nf!=3||!codes[2]:opt->nf<2||opt->nf>3)||
            opt->bds_if!=2||opt->mode<PMODE_PPP_KINEMA||
            opt->mode>PMODE_PPP_FIXED||opt->sateph!=EPHOPT_B2b||
            opt->ionoopt!=IONOOPT_IFLC||opt->modear!=ARMODE_OFF||
            (opt->navsys&~(SYS_GPS|SYS_CMP|SYS_GAL))) {
            strcpy(msg,"IF requires GPS/BDS B2b plus explicit GAL broadcast, float, nf=2/3 matching pairs");
            return 0;
        }
        return 1;
    }
    if (opt->bds_if==2) {
        strcpy(msg,"bds_if=2 requires if_model=1 or 2"); return 0;
    }
    if (opt->bds_if==1 && (opt->nf!=2 ||
        (opt->mode!=PMODE_SINGLE &&
         (opt->mode<PMODE_PPP_KINEMA||opt->sateph!=EPHOPT_B2b||opt->ionoopt!=IONOOPT_IFLC)))) {
        strcpy(msg,"B1C/B2a requires nf=2 and SPP or PPP-B2b IF mode");
        return 0;
    }
    return 1;
}

/* Append the processing selection after receiver options (last token wins). */
extern void bds_decode_options(char *dst, size_t size, const char *src, int pair)
{
    snprintf(dst,size,"%.*s -CBDSIF=%d",(int)size-20,src?src:"",pair);
}

extern void bds_decode_options_ex(char *dst, size_t size, const char *src, const prcopt_t *opt)
{
    uint8_t codes[3];
    if (opt->bds_if!=2) {bds_decode_options(dst,size,src,opt->bds_if);return;}
    if (!bds_parse_freqs(opt->bds_freqs,codes)) {if (size) dst[0]=0;return;}
    snprintf(dst,size,"%.*s -CBDSIF=2 -CBDSFREQ=%c%c%s",
        size>40?(int)size-40:0,src?src:"",code2obs(codes[0])[0],
        code2obs(codes[1])[0],codes[2]?(char[2]){code2obs(codes[2])[0],0}:"");
}

extern int code2obsidx(int sys, uint8_t code, const char *opt)
{
    const char *p;
    int pair=0,value,idx=code2idx(sys,code);
    if (sys!=SYS_CMP) return idx;
    for (p=opt;p&&(p=strstr(p,"-CBDSIF="));p++) {
        if (sscanf(p,"-CBDSIF=%d",&value)==1) pair=value;
    }
    if (pair==2) {
        uint8_t codes[3]={CODE_L5P,CODE_L1P,CODE_L2I};int i,j,map[6];
        for (p=opt;p&&(p=strstr(p,"-CBDSFREQ="));p++) {
            const char *s=p+10;
            if (strlen(s)<2) return -1;
            for (i=0;i<3;i++) {
                if (i==2&&(!s[i]||isspace((unsigned char)s[i]))) {codes[i]=0;break;}
                codes[i]=s[i]=='2'?CODE_L2I:s[i]=='6'?CODE_L6I:
                         s[i]=='1'?CODE_L1P:s[i]=='5'?CODE_L5P:0;
                if (!codes[i]) return -1;
                for (j=0;j<i;j++) if (codes[i]==codes[j]) return -1;
            }
            if (i==3&&s[3]&&!isspace((unsigned char)s[3])) return -1;
        }
        bds_slot_map(codes,map);
        return idx>=0&&idx<6?map[idx]:-1;
    }
    if (pair!=1) return idx;
    switch (idx) {
        case 5: return 0; /* B1C */
        case 3: return 1; /* B2a */
        case 0: return 5; /* B1I: extended */
        case 1: return 3; /* B3I: extended */
        default: return idx;
    }
}

static void copy_signal(obsd_t *dst, int f, const obsd_t *src, int j)
{
    dst->code[f]=j<0?0:src->code[j];dst->P[f]=j<0?0:src->P[j];
    dst->L[f]=j<0?0:src->L[j];dst->D[f]=j<0?0:src->D[j];
    dst->SNR[f]=j<0?0:src->SNR[j];dst->LLI[f]=j<0?0:src->LLI[j];
    dst->Lstd[f]=j<0?0:src->Lstd[j];dst->Pstd[f]=j<0?0:src->Pstd[j];
}

extern void bds_select_obs_ex(const obsd_t *src, obsd_t *dst, int n, const prcopt_t *opt)
{
    uint8_t codes[3];int map[6],i,j,f,idx,best,score,maxscore;
    if (opt->bds_if!=2) {bds_select_obs(src,dst,n,opt->bds_if);return;}
    if (!bds_parse_freqs(opt->bds_freqs,codes)) {memset(dst,0,n*sizeof(*dst));return;}
    bds_slot_map(codes,map);
    for (i=0;i<n;i++) {
        obsd_t in=src[i];dst[i]=in;
        if (satsys(in.sat,NULL)!=SYS_CMP) continue;
        /* Rebuild extended slots too: preserve B1C for SPP even when PPP omits it.
         * This permutation is idempotent, and never duplicates a raw signal. */
        for (f=0;f<NFREQ+NEXOBS;f++) {
            best=-1;maxscore=-1;
            for (j=0;j<NFREQ+NEXOBS;j++) {
                const char *s=code2obs(in.code[j]);
                idx=code2idx(SYS_CMP,in.code[j]);
                if (idx<0||idx>=6||map[idx]!=f||!in.P[j]) continue;
                if (f<(codes[2]?3:2)&&(s[0]=='2'||s[0]=='6'?s[1]!='I':s[1]!='P'&&s[1]!='D')) continue;
                score=(in.L[j]!=0.0?100:0)+getcodepri(SYS_CMP,in.code[j],"");
                if (score>maxscore) {best=j;maxscore=score;}
            }
            copy_signal(dst+i,f,&in,best);
        }
    }
}

extern void bds_select_spp_obs(const obsd_t *src, obsd_t *dst, int n, const nav_t *nav)
{
    int i,j,best,score,maxscore;double bias;
    for (i=0;i<n;i++) {
        obsd_t in=src[i];dst[i]=in;
        if (satsys(in.sat,NULL)!=SYS_CMP) continue;
        best=-1;maxscore=-1;
        for (j=0;j<NFREQ+NEXOBS;j++) {
            if (!in.P[j]) continue;
            score=in.code[j]==CODE_L6I?1:0; /* B3I is the broadcast clock datum. */
            if ((in.code[j]==CODE_L1P||in.code[j]==CODE_L1D)&&
                bds_tgd_bias(in.time,in.sat,nav,in.code[j],&bias))
                score=in.code[j]==CODE_L1P?3:2;
            if (score>maxscore&&score>0) {best=j;maxscore=score;}
        }
        copy_signal(dst+i,0,&in,best);
        if (best>0) copy_signal(dst+i,best,&in,0);
    }
}

/* Normalize external/raw observations too. A missing selected signal remains
 * missing: never silently substitute a legacy signal in the IF combination. */
extern void bds_select_obs(const obsd_t *src, obsd_t *dst, int n, int pair)
{
    int i,j,f,best,score,maxscore;
    obsd_t in;
    for (i=0;i<n;i++) {
        in=src[i]; dst[i]=in;
        if (!pair||satsys(in.sat,NULL)!=SYS_CMP) continue;
        for (f=0;f<(pair==2?3:2);f++) {
            char band=pair==2?"512"[f]:"15"[f];
            best=-1; maxscore=-1;
            for (j=0;j<NFREQ+NEXOBS;j++) {
                const char *s=code2obs(in.code[j]);
                if (s[0]!=band||!in.P[j]) continue;
                if (band=='2'?s[1]!='I':(s[1]!='P'&&s[1]!='D')) continue;
                score=(s[1]=='P'?2:1)+(in.L[j]!=0.0?4:0);
                if (score>maxscore) {best=j;maxscore=score;}
            }
            dst[i].code[f]=best<0?0:in.code[best];
            dst[i].P[f]=best<0?0.0:in.P[best];
            dst[i].L[f]=best<0?0.0:in.L[best];
            dst[i].D[f]=best<0?0.0f:in.D[best];
            dst[i].SNR[f]=best<0?0:in.SNR[best];
            dst[i].LLI[f]=best<0?0:in.LLI[best];
            dst[i].Lstd[f]=best<0?0:in.Lstd[best];
            dst[i].Pstd[f]=best<0?0:in.Pstd[best];
        }
        trace(4,"bds_select: sat=%d pair=%d code=%d/%d/%d P=%.3f/%.3f\n",
              in.sat,pair,dst[i].code[0],dst[i].code[1],dst[i].code[2],dst[i].P[0],dst[i].P[1]);
    }
}

/* Independent physical calibration lookup. Receiver fallback is opt-in and
 * uses the same RF frequencies, not B1I/B3I's different frequencies. */
extern int bds_ant_index(const pcv_t *pcv, uint8_t code, int fallback)
{
    const char *s=code2obs(code);
    int idx=s[0]=='1'?ANT_B1C:s[0]=='5'?ANT_B2A:-1;
    if (idx<0) return -1;
    if (pcv->valid[idx]) return idx;
    if (!pcv->sat&&fallback) {
        idx=3*NFREQ+(s[0]=='1'?0:1); /* E01/E05 */
        if (pcv->valid[idx]) return idx;
    }
    return -1;
}

extern int bds_code_bias(gtime_t time, const B2bssr_t *ssr, uint8_t code, double *bias)
{
    double age=timediff(time,ssr->t0[1]);
    if (code==CODE_L6I) {
        /* Known reference-code zero, not a missing non-reference DCB.
         * Require a current, internally consistent PPP-B2b orbit/clock datum. */
        double orb=timediff(time,ssr->t0[0]),clk=timediff(time,ssr->t0[2]);
        if (!ssr->t0[0].time||!ssr->t0[2].time||orb < -1||orb>126||clk < -1||clk>42||
            ssr->iodssr[0]!=ssr->iodssr[2]||ssr->iodcorr[0]!=ssr->iodcorr[1]) return 0;
        *bias=0.0;return 1;
    }
    if (!code||code>=MAXCODE||!ssr->cbias_valid[code]||!ssr->t0[1].time||
        age < -1.0||age > 86400.0||ssr->iodssr[1]!=ssr->iodssr[0]||
        ssr->iodssr[1]!=ssr->iodssr[2]) return 0;
    *bias=ssr->cbias[code];
    return 1;
}

extern int bds_tgd_bias(gtime_t time, int sat, const nav_t *nav, uint8_t code, double *bias)
{
    const char *s=code2obs(code);
    const eph_t *eph=NULL;
    int i,t=s[0]=='1'?2:s[0]=='5'?3:-1,mask;
    double dt,best=MAXDTOE_CMP+1.0;
    if (t<0||(s[1]!='P'&&s[1]!='D')) return 0;
    mask=(1<<t)|(s[1]=='D'?(1<<(t+2)):0);
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat||(nav->eph[i].tgd_valid&mask)!=mask) continue;
        dt=fabs(timediff(time,nav->eph[i].toe));
        if (dt<=MAXDTOE_CMP && dt<best) {best=dt;eph=nav->eph+i;}
    }
    if (!eph) return 0;
    *bias=CLIGHT*(eph->tgd[t]+(s[1]=='D'?eph->tgd[t+2]:0.0));
    return 1;
}
