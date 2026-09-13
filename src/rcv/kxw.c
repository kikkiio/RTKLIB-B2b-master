/*------------------------------------------------------------------------------
* kxw.c : decoder for KXW KSRTCM satellite correction messages
*
* The transport is RTCM-like message 4047. KS message 0_40 carries one
* PPP-B2b frame as 16 little-endian UINT32 words. The words are converted to
* network bit order. V2.3 0_01 CNAV1 ephemerides and 0_40 B2b corrections
* are decoded with per-receiver assembly/mask state.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "B2b.h"

#define KXW_MSGTYPE_PPPB2B  0
#define KXW_MSGID_PPPB2B 0x40
#define KXW_MAX_PAYLOAD  16304
#define KXW_PPPB2B_LEN      67

typedef struct {
    B2bmask_t mask; /* per-stream mask: never share state with another receiver */
    uint8_t payload[KXW_MAX_PAYLOAD];
    int len;
    int msgtype;
    int msgid;
    int last_part;
    int next_part;
    int assembling;
} kxw_t;

/* KSRTCM scalar fields are little endian, including unaligned doubles. */
static unsigned int kxw_u2(const uint8_t *p) {return p[0]|((unsigned int)p[1]<<8);}
static int kxw_i4(const uint8_t *p)
{
    uint32_t u=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);
    int32_t v;memcpy(&v,&u,4);return v;
}
static double kxw_r8(const uint8_t *p)
{
    uint64_t u=0;double v;int i;
    for (i=0;i<8;i++) u|=(uint64_t)p[i]<<(8*i);
    memcpy(&v,&u,8);return v;
}
static int decode_cnav1(raw_t *raw,const uint8_t *p,int len)
{
    eph_t e={0};double v[20];int i,sat,week,toc,toe;
    if (len!=178) return -1;
    week=(int16_t)kxw_u2(p+6);toc=kxw_i4(p+10);toe=kxw_i4(p+14);
    if (!(sat=satno(SYS_CMP,p[0]))||p[2]<1||p[2]>3||week<0||
        toe<0||toe>=604800||toc<0||toc>=604800) return -1;
    for (i=0;i<20;i++) {v[i]=kxw_r8(p+18+8*i);if (!isfinite(v[i])) return -1;}
    e.sat=sat;e.iode=p[1];e.iodc=kxw_u2(p+4);e.week=week;
    e.flag=BDS_NAV_KXW_CNV1;e.toes=toe;
    e.toe=bdt2gpst(bdt2time(week,toe));e.toc=bdt2gpst(bdt2time(week,toc));
    e.ttr=raw->time;
    e.f0=v[0];e.f1=v[1];e.f2=v[2];e.e=v[3];e.omg=v[4];e.M0=v[5];
    e.OMG0=v[6];e.OMGd=v[7];e.i0=v[8];e.idot=v[9];
    e.cuc=v[10];e.cus=v[11];e.crc=v[12];e.crs=v[13];e.cic=v[14];e.cis=v[15];
    e.A=(p[2]==3?27906100.0:42162200.0)+v[16];
    e.Adot=v[17];e.deln=v[18];e.ndot=v[19];
    if (e.A<2E7||e.A>5E7||e.e<0.0||e.e>=1.0) return -1;
    /* V2.3 supplies no TGD/ISC, health or URA. Do not invent valid TGD.
       PPP uses the B2b URA; code bias uses actual B2b DCB presence. */
    e.sva=15;e.tgd_valid=0;
    raw->raw_nmsg[9]++;
    if (!strstr(raw->opt,"-EPHALL")&&raw->nav.eph[sat-1].sat==sat&&
        raw->nav.eph[sat-1].iode==e.iode&&raw->nav.eph[sat-1].iodc==e.iodc&&
        timediff(raw->nav.eph[sat-1].toe,e.toe)==0.0&&
        timediff(raw->nav.eph[sat-1].toc,e.toc)==0.0) return 0;
    raw->nav.eph[sat-1]=e;raw->ephsat=sat;raw->ephset=0;
    trace(3,"KXW CNAV1 C%02d IODE=%d IODC=%d toe=%s\n",p[0],e.iode,e.iodc,time_str(e.toe,0));
    if (raw->outtype) sprintf(raw->msgtype,"KXW 4047 0_01: C%02d IODC=%d",p[0],e.iodc);
    return 2;
}

static int decode_pppb2b(raw_t *raw, const uint8_t *payload, int len)
{
    kxw_t *state=(kxw_t *)raw->rcv_data;
    B2bmask_t *mask=&state->mask;
    uint8_t frame[64];B2bssr_t *ssr;
    int i,j,k=29,type,sod,iod,iodp,sub,slot,sat,sys,n,nsig,mode,bias,changed=0;
    int iodn,iodcorr,radial,along,cross,ura,*codes;
    gtime_t time;
    if (len!=KXW_PPPB2B_LEN) return -1;
    if (payload[1]&0x20) return 0;
    for (i=0;i<16;i++) for (j=0;j<4;j++) frame[4*i+j]=payload[3+4*i+3-j];
    type=payload[2]&63;
    if ((int)getbitu(frame,0,6)!=type) return -1;
    /* Null/unsupported messages have no usable time-of-day or observations. */
    if (type<1||type>4) {raw->raw_nmsg[7]++;return 0;}
    sod=getbitu(frame,6,17);
    if (sod>=86400||!raw->time.time) return -1;
    time=B2btod2time(raw->time,sod);iod=getbitu(frame,27,2);
    /* Validate variable-length DCB records before changing any state. */
    if (type==3) {
        n=getbitu(frame,k,5);k+=5;
        for (i=0;i<n;i++) {
            if (k+13>462) return -1;
            nsig=getbitu(frame,k+9,4);k+=13+16*nsig;
            if (k>462) return -1;
        }
    }
    raw->time=time;raw->geoprn=payload[0];k=29;
    raw->num_PPPB2BINF01=raw->num_PPPB2BINF02=0;
    raw->num_PPPB2BINF03=raw->num_PPPB2BINF04=0;
    raw->raw_nmsg[type-1]++;
    if (raw->outtype) sprintf(raw->msgtype,"KXW 4047 0_40: PRN=%d type=%d",payload[0],type);
    if (type==1) {
        iodp=getbitu(frame,k,4);k+=4;
        if (iod!=mask->IOD_SSR) memset(raw->nav.B2bssr,0,sizeof(raw->nav.B2bssr));
        else if (iodp!=mask->IODP) for (i=0;i<MAXSAT;i++) raw->nav.B2bssr[i].t0[2].time=0;
        memset(mask,0,sizeof(*mask));mask->IOD_SSR=iod;mask->IODP=iodp;
        mask->time=mask->m_time=time;
        for (i=0;i<63;i++) mask->MASK_BD[i]=getbitu(frame,k++,1);
        for (i=0;i<37;i++) mask->MASK_GPS[i]=getbitu(frame,k++,1);
        for (i=0;i<37;i++) mask->MASK_GALILEO[i]=getbitu(frame,k++,1);
        for (i=0;i<37;i++) mask->MASK_GLONASS[i]=getbitu(frame,k++,1);
        mask2satno(mask);raw->num_PPPB2BINF01=1;output_B2bInfo1(raw,mask,0);return 20;
    }
    if (iod!=mask->IOD_SSR) return 0;
    if (type==2) {
        for (i=0;i<6;i++) {
            slot=getbitu(frame,k,9);k+=9;iodn=getbitu(frame,k,10);k+=10;
            iodcorr=getbitu(frame,k,3);k+=3;radial=getbits(frame,k,15);k+=15;
            along=getbits(frame,k,13);k+=13;cross=getbits(frame,k,13);k+=13;
            ura=getbitu(frame,k,6);k+=6; /* consume even invalid/empty slots */
            sat=slot2satno(slot);
            if (sat<=0||sat>=MAXSAT||abs(radial)>=16383||abs(along)>=4095||abs(cross)>=4095) continue;
            ssr=raw->nav.B2bssr+sat;ssr->t0[0]=time;ssr->sow=sod;
            ssr->iodssr[0]=iod;ssr->iodn=iodn;ssr->iodcorr[0]=iodcorr;ssr->ura=ura;
            ssr->deph[0]=radial*0.0016;ssr->deph[1]=along*0.0064;ssr->deph[2]=cross*0.0064;
            ssr->update=1;changed++;
        }
        if (changed) {raw->num_PPPB2BINF02=1;output_B2bInfo2(raw,&raw->nav);}
    }
    if (type==3) {
        n=getbitu(frame,k,5);k+=5;
        for (i=0;i<n;i++) {
            slot=getbitu(frame,k,9);k+=9;nsig=getbitu(frame,k,4);k+=4;
            sat=slot2satno(slot);sys=satsys(sat,NULL);
            codes=sys==SYS_CMP?b2b_bds_codebias_mode:sys==SYS_GPS?b2b_gps_codebias_mode:
                  sys==SYS_GAL?b2b_gal_codebias_mode:sys==SYS_GLO?b2b_glo_codebias_mode:NULL;
            ssr=sat>0&&sat<MAXSAT&&codes?raw->nav.B2bssr+sat:NULL;
            if (ssr) {memset(ssr->cbias_valid,0,sizeof(ssr->cbias_valid));ssr->t0[1]=time;ssr->iodssr[1]=iod;}
            for (j=0;j<nsig;j++) {
                mode=getbitu(frame,k,4);k+=4;bias=getbits(frame,k,12);k+=12;
                if (!ssr||mode>=B2B_CodeBiasModeNum||!codes[mode]||abs(bias)>=2047) continue;
                ssr->cbias[codes[mode]]=bias*0.017;ssr->cbias_valid[codes[mode]]=1;
            }
            if (ssr) {ssr->update=1;changed++;}
        }
        if (changed) {raw->num_PPPB2BINF03=1;output_B2bInfo3(raw,&raw->nav);}
    }
    if (type==4) {
        iodp=getbitu(frame,k,4);k+=4;sub=getbitu(frame,k,5);k+=5;
        if (iodp!=mask->IODP||sub*23>=B2B_MAXSAT) return 0;
        for (i=0;i<23;i++) {
            iodcorr=getbitu(frame,k,3);k+=3;bias=getbits(frame,k,15);k+=15;
            j=sub*23+i;if (j>=B2B_MAXSAT) continue;
            sat=mask->satno[j];if (sat<=0||sat>=MAXSAT||abs(bias)>=16383) continue;
            ssr=raw->nav.B2bssr+sat;ssr->t0[2]=time;ssr->sow=sod;
            ssr->iodssr[2]=iod;ssr->iodp[0]=iodp;ssr->iodcorr[1]=iodcorr;
            ssr->dclk[0]=bias*0.0016;ssr->update=1;changed++;
        }
        if (changed) {raw->num_PPPB2BINF04=1;output_B2bInfo4(raw,&raw->nav);}
    }
    return changed?20:0;
}

static int decode_packet(raw_t *raw, int msgtype, int msgid,
                         const uint8_t *payload, int len)
{
    if (msgtype == KXW_MSGTYPE_PPPB2B && msgid == KXW_MSGID_PPPB2B) {
        return decode_pppb2b(raw, payload, len);
    }
    if (msgtype==0&&msgid==1) return decode_cnav1(raw,payload,len);
    trace(3, "KXW unsupported message: type=%d id=%d len=%d\n",
          msgtype, msgid, len);
    return 0;
}

static int decode_kxw(raw_t *raw)
{
    kxw_t *kxw = (kxw_t *)raw->rcv_data;
    uint32_t crc, crc_msg;
    int data_len, msgtype, msgid, part, last_part, payload_len;
    const uint8_t *payload;

    data_len = ((raw->buff[1] & 0x03) << 8) | raw->buff[2];
    crc = rtk_crc24q(raw->buff, data_len + 3);
    crc_msg = getbitu(raw->buff, (data_len + 3) * 8, 24);
    if (crc != crc_msg) {
        trace(2, "KXW CRC-24Q error: len=%d crc=%06X expected=%06X\n",
              data_len, crc, crc_msg);
        kxw->assembling = 0;
        return -1;
    }
    if (raw->buff[3] != 0xFC || (raw->buff[4] & 0xF0) != 0xF0) {
        trace(2, "KXW message ID error: %02X %02X\n",
              raw->buff[3], raw->buff[4]);
        kxw->assembling = 0;
        return -1;
    }

    msgtype = raw->buff[4] & 0x0F;
    msgid = raw->buff[5];
    part = raw->buff[6] >> 4;
    last_part = raw->buff[6] & 0x0F;
    payload = raw->buff + 7;
    payload_len = data_len - 4;

    if (kxw->assembling && (msgtype!=kxw->msgtype || msgid!=kxw->msgid ||
        last_part!=kxw->last_part || part!=kxw->next_part)) {
        kxw->assembling=0;return -1; /* V2.3: discard the offending packet */
    }
    if (last_part == 0) {
        if (part != 0) return -1;
        kxw->assembling = 0;
        return decode_packet(raw, msgtype, msgid, payload, payload_len);
    }

    if (part == 0 && !kxw->assembling) {
        kxw->len = 0;
        kxw->msgtype = msgtype;
        kxw->msgid = msgid;
        kxw->last_part = last_part;
        kxw->next_part = 0;
        kxw->assembling = 1;
    }
    if (!kxw->assembling || msgtype != kxw->msgtype || msgid != kxw->msgid ||
        last_part != kxw->last_part || part != kxw->next_part ||
        kxw->len + payload_len > KXW_MAX_PAYLOAD) {
        trace(2, "KXW packet sequence error: type=%d id=%d part=%d/%d\n",
              msgtype, msgid, part, last_part);
        kxw->assembling = 0;
        return -1;
    }
    memcpy(kxw->payload + kxw->len, payload, payload_len);
    kxw->len += payload_len;
    kxw->next_part++;

    if (part != last_part) return 0;

    kxw->assembling = 0;
    return decode_packet(raw, kxw->msgtype, kxw->msgid,
                         kxw->payload, kxw->len);
}

extern int init_kxw(raw_t *raw)
{
    raw->rcv_data = calloc(1, sizeof(kxw_t));
    if (!raw->rcv_data) return 0;
    ((kxw_t *)raw->rcv_data)->mask.IOD_SSR=-1;
    ((kxw_t *)raw->rcv_data)->mask.IODP=-1;
    return 1;
}

extern void free_kxw(raw_t *raw)
{
    free(raw->rcv_data);
    raw->rcv_data = NULL;
}

extern int input_kxw(raw_t *raw, uint8_t data)
{
    int data_len, total_len;

    if (raw->nbyte == 0) {
        if (data != 0xD3) return 0;
        raw->buff[raw->nbyte++] = data;
        return 0;
    }
    raw->buff[raw->nbyte++] = data;

    if (raw->nbyte == 3) {
        data_len = ((raw->buff[1] & 0x03) << 8) | raw->buff[2];
        if ((raw->buff[1] & 0xFC) || data_len < 4 || data_len > 1023) {
            trace(2, "KXW length error: len=%d\n", data_len);
            raw->nbyte = 0;
            ((kxw_t *)raw->rcv_data)->assembling=0;
            return -1;
        }
        raw->len = data_len + 3;
    }
    if (raw->nbyte < 3) return 0;

    total_len = raw->len + 3;
    if (raw->nbyte < total_len) return 0;
    raw->nbyte = 0;
    return decode_kxw(raw);
}

extern int input_kxwf(raw_t *raw, FILE *fp)
{
    int i, data, ret;

    if (!raw || !fp) return -1;
    for (i = 0; i < 4096; i++) {
        data = fgetc(fp);
        if (data == EOF) return -2;
        if ((ret = input_kxw(raw, (uint8_t)data))) return ret;
    }
    return 0;
}
