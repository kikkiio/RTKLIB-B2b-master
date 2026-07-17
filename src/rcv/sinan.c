/*------------------------------------------------------------------------------
* sinan.c : Decoder for Shanghai Sinan custom receiver format
*
*           Developed by Xu@CUGB and Liu@APM
*
* description : This file contains the implementation for parsing Shanghai
*               Sinan company receiver's custom binary format. The decoder
*               includes strategies for data reading and partial decoding.
* options : None
*
* version : $Revision:$ $Date:$
* history : 2024/10/30   Initial version by Xu@CUGB
*           2024/11/19   Improved reading strategies and decoding methods,
*                        changed cpp to c by Liu@APM
*
*-----------------------------------------------------------------------------*/


#include "rtklib.h"
#include "B2b.h"


#define SINAN_SYNC1       0xAA    /* message start sync code 1 */
#define SINAN_SYNC2       0x44    /* message start sync code 2 */
#define SINAN_SYNC3       0x12    /* message start sync code 3 */

#define SINAN_LEN         28      /* message header length (bytes) */
#define MAX_SINAN_LEN     16384     /* max message length:H(28)+D(127)+CRC(4) (bytes) */
// Liu@APM:B2b length is too short for another message


#define B2BRAWNAVSUBFRAMEb 1697   /*sinan B2b type */
#define GPSEPHEM 71
#define BD3EPHEM 72

/* KXW/SSR outer frame and standard PPP-B2b page layout. */
#define SSR3PREAMB               0xD3
#define KXW_PPPB2B_MSG_ID        0x40
#ifndef KXW_PPPB2B_PAGE_OFFSET
#define KXW_PPPB2B_PAGE_OFFSET   6
#endif
#define PPPB2B_PREAMBLE           0xEB90
#define PPPB2B_DATA_BITS          456
#define PPPB2B_CRC_BITS            24
#define PPPB2B_MESSAGE_BITS       (6+PPPB2B_DATA_BITS+PPPB2B_CRC_BITS)
#define PPPB2B_CRC_INPUT_BITS     (6+PPPB2B_DATA_BITS)
#define PPPB2B_PACKED_BITS        512
#define PPPB2B_PACKED_BYTES       (PPPB2B_PACKED_BITS/8)
#define PPPB2B_INFO_BITS          (6+6+6+PPPB2B_DATA_BITS)
#define RTCM_PRIVATE_TYPE         4047
#define RTCM_PPPB2B_SUBTYPE       64
#define RTCM_PPPB2B_PRN_OFFSET    7
#define RTCM_PPPB2B_RSV_OFFSET    8
#define RTCM_PPPB2B_TYPE_OFFSET   9
#define RTCM_PPPB2B_BODY_OFFSET  10

static uint16_t U2(uint8_t *p) {uint16_t u; memcpy(&u,p,2); return u;}
static uint32_t U4(uint8_t *p) {uint32_t u; memcpy(&u,p,4); return u;}

/* CRC24Q over arbitrary MSB-first bits. PPP-B2b protects the non-byte-aligned
 * MesTypeID(6)+data(456) block. */
static uint32_t crc24q_bits(const uint8_t *buff, int pos, int len)
{
    uint32_t crc=0,feedback;
    int i;

    for (i=0;i<len;i++) {
        feedback=((crc>>23)&1u)^getbitu(buff,pos+i,1);
        crc=(crc<<1)&0xFFFFFFu;
        if (feedback) crc^=0x864CFBu;
    }
    return crc;
}

B2bmask_t sinan_mask ={
    .MASK_BD = {0},
    .MASK_GPS = {0},
    .MASK_GALILEO = {0},
    .MASK_GLONASS = {0},
    .IOD_SSR = -1,
    .IODP = -1,
    .satnum = -1,
    .satno = {0}
};

/* A recorded RTCM stream can interleave pages broadcast by several PPP-B2b
 * GEO satellites.  Their IOD transitions are not necessarily simultaneous,
 * so sharing the receiver-specific sinan_mask would incorrectly reject valid
 * type 2/3/4 pages when another broadcaster changes its mask first. */
static B2bmask_t rtcm_b2b_mask[64];
static int rtcm_b2b_mask_initialized = 0;

static B2bmask_t *get_rtcm_b2b_mask(uint32_t prn)
{
    int i;

    if (prn < 1 || prn > 63) return NULL;
    if (!rtcm_b2b_mask_initialized) {
        memset(rtcm_b2b_mask, 0, sizeof(rtcm_b2b_mask));
        for (i = 0; i < 64; i++) {
            rtcm_b2b_mask[i].IOD_SSR = -1;
            rtcm_b2b_mask[i].IODP = -1;
        }
        rtcm_b2b_mask_initialized = 1;
    }
    return &rtcm_b2b_mask[prn];
}

/* sync header ---------------------------------------------------------------*/
static int sync_sinan(uint8_t *buff, uint8_t data)
{
    buff[0]=buff[1]; buff[1]=buff[2]; buff[2]=data;
    return buff[0]==SINAN_SYNC1&&buff[1]==SINAN_SYNC2&&buff[2]==SINAN_SYNC3;
}

static int process_message_type_1(raw_t *raw, int* pose)
{
    const uint8_t* buffer = raw->buff;
    int num,udi = 0;
    gtime_t temp_time = {0};

    temp_time = sinan_mask.time;
    
    memset(&sinan_mask.time, 0, sizeof(sinan_mask.time));
    memset(&sinan_mask.m_time, 0, sizeof(sinan_mask.m_time));
    memset(sinan_mask.MASK_BD, 0, sizeof(sinan_mask.MASK_BD));
    memset(sinan_mask.MASK_GPS, 0, sizeof(sinan_mask.MASK_GPS));
    memset(sinan_mask.MASK_GALILEO, 0, sizeof(sinan_mask.MASK_GALILEO));
    memset(sinan_mask.MASK_GLONASS, 0, sizeof(sinan_mask.MASK_GLONASS));
    sinan_mask.IOD_SSR = sinan_mask.IODP = 0;

    uint32_t sow = getbitu(buffer, *pose, 17); *pose += 17;
    sinan_mask.time = B2btod2time(raw->time,sow);
    sinan_mask.m_time = raw->time;
    
    if (temp_time.time == 0 && temp_time.sec == 0.0) {
        udi = 0;
    } else {
        udi = timediff(sinan_mask.time,temp_time);
    }

    uint32_t reserve2 = getbitu(buffer, *pose, 4); *pose += 4;
    sinan_mask.IOD_SSR = getbitu(buffer, *pose, 2); *pose += 2;
    sinan_mask.IODP = getbitu(buffer, *pose, 4); *pose += 4;

    for (int i = 0; i < 63; i++) {
        sinan_mask.MASK_BD[i] = getbitu(buffer, *pose, 1);
        *pose += 1;
    }

    for (int i = 0; i < 37; i++) {
        sinan_mask.MASK_GPS[i] = getbitu(buffer, *pose, 1);
        *pose += 1;
    }

    for (int i = 0; i < 37; i++) {
        sinan_mask.MASK_GALILEO[i] = getbitu(buffer, *pose, 1);
        *pose += 1;
    }

    for (int i = 0; i < 37; i++) {
        sinan_mask.MASK_GLONASS[i] = getbitu(buffer, *pose, 1);
        *pose += 1;
    }

    num = mask2satno(&sinan_mask);

    /* Reserved fields are longer than the 32-bit getbitu() return type. */
    *pose += 81;
    *pose += 174;

    output_B2bInfo1(raw, &sinan_mask,udi);
    return 20;

}

static int process_message_type_2(raw_t *raw, int* pose)
{
    const uint8_t* buffer = raw->buff;
    int satno,satslot;
    char satid[8];
    gtime_t m_time;
    double m_time_ep[8];
    

	uint32_t sow = getbitu(buffer, *pose, 17); *pose += 17;
    m_time = B2btod2time(raw->time,sow);
    time2epoch(m_time,m_time_ep);
    int verify_sow = (int)(m_time_ep[3]*60*60 + m_time_ep[4]*60 + m_time_ep[5]);

	uint32_t reserve2 = getbitu(buffer, *pose, 4); *pose += 4;

	uint32_t IOD_SSR = getbitu(buffer, *pose, 2); *pose += 2;

    if (IOD_SSR != sinan_mask.IOD_SSR) {
        // printf("error: msg2(eph), new iodssr=%d,msg->iodssr=%d\n", 
        // IOD_SSR, sinan_mask.IOD_SSR);
        trace(22,"error: msg2(eph), new iodssr=%d,msg->iodssr=%d\n", 
        IOD_SSR, sinan_mask.IOD_SSR);
        return 0;
    }

	for (int i = 0; i < 6; i++) {
		satslot = getbitu(buffer, *pose, 9); *pose += 9;
        satno = slot2satno(satslot);

        raw->nav.B2bssr[satno].iodssr[0] = IOD_SSR;
        raw->nav.B2bssr[satno].t0[0] = m_time;
        raw->nav.B2bssr[satno].sow = sow;
        raw->nav.B2bssr[satno].verify_sow = verify_sow;

        raw->nav.B2bssr[satno].iodn = getbitu(buffer, *pose, 10); *pose += 10;

        raw->nav.B2bssr[satno].iodcorr[0] = getbitu(buffer, *pose, 3); *pose += 3;
        
        int sRadial = getbits(buffer, *pose, 15); *pose += 15;
        
        int sInTrack = getbits(buffer, *pose, 13); *pose += 13;
        
        int sCross = getbits(buffer, *pose, 13); *pose += 13;
        
        if (abs(sRadial) >= 16383 ||abs(sInTrack) >= 4095 ||abs(sCross) >= 4095) {
                satno2id(satno,satid);
                trace(22,"error: bad-orb, sat=%s,deph=%8.3f,%8.3f,%8.3f\n", 
                satid, raw->nav.B2bssr[satno].deph[0], raw->nav.B2bssr[satno].deph[1], raw->nav.B2bssr[satno].deph[2]);
                // printf("error: bad-orb, sat=%s,deph=%8.3f,%8.3f,%8.3f\n", 
                // satid, raw->nav.B2bssr[satno].deph[0], raw->nav.B2bssr[satno].deph[1], raw->nav.B2bssr[satno].deph[2]);
                continue;
        }

        raw->nav.B2bssr[satno].deph[0] = sRadial * 0.0016; 
        raw->nav.B2bssr[satno].deph[1] = sInTrack* 0.0064;
        raw->nav.B2bssr[satno].deph[2] = sCross* 0.0064;

        raw->nav.B2bssr[satno].ura = getbitu(buffer, *pose, 6); *pose += 6;

        raw->nav.B2bssr[satno].update = 1;
	}
	uint32_t reserve_info2 = getbitu(buffer, *pose, 6); *pose += 6;
    output_B2bInfo2(raw,&raw->nav);
    return 20;
}

static int process_message_type_3(raw_t *raw, int* pose)
{
    const uint8_t* buffer = raw->buff;
    int satno,satslot,sys,type,*cods = NULL;
    double m_time_ep[8];
    char satid[8];
    gtime_t m_time;

    uint32_t sow = getbitu(buffer, *pose, 17); *pose += 17;
    m_time = B2btod2time(raw->time,sow);
    time2epoch(m_time,m_time_ep);
    int verify_sow = (int)(m_time_ep[3]*60*60 + m_time_ep[4]*60 + m_time_ep[5]);


    uint32_t reserve2 = getbitu(buffer, *pose, 4); *pose += 4;
    uint32_t IOD_SSR = getbitu(buffer, *pose, 2); *pose += 2;
    if (IOD_SSR != sinan_mask.IOD_SSR) {
        // printf("error: msg3(cbia), new iodssr=%d,msg->iodssr=%d\n", 
        // IOD_SSR, sinan_mask.IOD_SSR);
        trace(22,"error: msg3(cbia), new iodssr=%d,msg->iodssr=%d\n", 
        IOD_SSR, sinan_mask.IOD_SSR);
        return 0;
    }

    uint32_t SatNum = getbitu(buffer, *pose, 5); *pose += 5;

    for (int i = 0; i < SatNum; i++) {
        satslot = getbitu(buffer, *pose, 9); *pose += 9;
        satno = slot2satno(satslot);
        satno2id(satno,satid);
        sys = satsys(satno, NULL);
        if (sys == SYS_GPS) cods = b2b_gps_codebias_mode;
        else if (sys == SYS_GLO) cods = b2b_glo_codebias_mode;
        else if (sys == SYS_GAL) cods = b2b_gal_codebias_mode;
        else if (sys == SYS_CMP) cods = b2b_bds_codebias_mode;
        else continue;

        raw->nav.B2bssr[satno].iodssr[1] = IOD_SSR;
        raw->nav.B2bssr[satno].t0[1] = m_time;
        raw->nav.B2bssr[satno].sow = sow;
        raw->nav.B2bssr[satno].verify_sow = verify_sow;

        uint32_t SigNum = getbitu(buffer, *pose, 4); *pose += 4;


        for (int j = 0; j < SigNum; j++) {
            uint32_t mode = getbitu(buffer, *pose, 4); *pose += 4;
            int32_t DCB = getbits(buffer, *pose, 12); *pose += 12;

            if (abs(DCB) >= 2103) {
                continue;
            }
            if (mode < 0 || mode >= B2B_CodeBiasModeNum) {
                continue;
            }
            type = cods[mode];
            if (type == CODE_NONE) continue;
            raw->nav.B2bssr[satno].cbias[type] = DCB*0.017;
            raw->nav.B2bssr[satno].update = 1;
        }
    }

    uint32_t reserve3 = getbitu(buffer, *pose, 10); *pose += 10;
    output_B2bInfo3(raw,&raw->nav);
    return 20;
}



static int process_message_type_4(raw_t *raw, int* pose)
{
    const uint8_t* buffer = raw->buff;
    int satno,satslot,sys,type,flag, i, *cods = NULL;
    double m_time_ep[8];
    char satid[8],time_str[64];
    gtime_t m_time;

    uint32_t sow = getbitu(buffer, *pose, 17); *pose += 17;
    m_time = B2btod2time(raw->time,sow);

    uint32_t reserve2 = getbitu(buffer, *pose, 4); *pose += 4;
    uint32_t IOD_SSR = getbitu(buffer, *pose, 2); *pose += 2;
    uint32_t IODP = getbitu(buffer, *pose, 4); *pose += 4;
    uint32_t SubTypel = getbitu(buffer, *pose, 5); *pose += 5;

     if (IOD_SSR != sinan_mask.IOD_SSR || IODP != sinan_mask.IODP) {
        // printf("error: msg4(clk), new iodssr=%d,msg->iodssr=%d,iodp=%d,msg->iodp=%d\n", 
        // IOD_SSR, sinan_mask.IOD_SSR, IODP, sinan_mask.IODP);
        trace(22,"error: msg4(clk), new iodssr=%d,msg->iodssr=%d,iodp=%d,msg->iodp=%d\n", 
        IOD_SSR, sinan_mask.IOD_SSR, IODP, sinan_mask.IODP);
        return 0;
    }


    if (SubTypel > 31) {
        trace(4, "error: msg4(clk), bad subtype1=%d\n", SubTypel);
    }

    time2str(m_time,time_str,3);

    int begin = SubTypel * 23;

    for (i = 0; i < 23; i++){
        if(sinan_mask.satno[begin+i] != 0){
            satno = sinan_mask.satno[begin+i];
            raw->nav.B2bssr[satno].t0[2] = m_time;
            raw->nav.B2bssr[satno].sow = sow;
            time2epoch(m_time,m_time_ep);
            raw->nav.B2bssr[satno].verify_sow = (int)(m_time_ep[3]*60*60 + m_time_ep[4]*60 + m_time_ep[5]);
            raw->nav.B2bssr[satno].iodssr[2] = IOD_SSR;
            raw->nav.B2bssr[satno].iodp[0] = IODP;
            raw->nav.B2bssr[satno].iodcorr[1] = getbitu(buffer, *pose, 3); *pose += 3;

            int sC0 = getbits(buffer, *pose, 15); *pose += 15;

            if (abs(sC0) >= 16383 || raw->nav.B2bssr[satno].iodcorr[1] > 7) {
                satno2id(satno,satid);
					trace(22, "error: msg4(clk) bad-clk, sat=%s,iodcor=%d,c0=%d \n", satid, 
                    raw->nav.B2bssr[satno].iodcorr[1], sC0);
					continue;
				}
            
            raw->nav.B2bssr[satno].dclk[0] = 0.0016 * sC0;
            raw->nav.B2bssr[satno].update = 1;
        }
    }
    uint32_t reserve3 = getbitu(buffer, *pose, 10); *pose += 10;
    output_B2bInfo4(raw,&raw->nav);
    return 20;
}

/* decode one standard PPP-B2b information page -----------------------------
 * The 4047/64 body consists of sixteen little-endian uint32_t words. After
 * restoring each word, its first 486 bits exactly match the SIS message before
 * LDPC: MesTypeID(6)+data(456)+CRC24Q(24); the remaining 26 bits are padding.
 * The legacy receiver path passes only an already extracted 456-bit domain. */
static int decode_PPPB2b_body(raw_t *raw, B2bmask_t *mask, uint32_t prn_6,
                              uint32_t status_6, uint32_t mes_type, int pose,
                              int endbit, int word_little_endian)
{
    uint32_t crc_rx,crc_calc;
    int i,n,srcbyte,srcbit,dstbyte,dstbit,message_end;

    if (!raw||!mask||pose<0||
        pose+(word_little_endian?PPPB2B_PACKED_BITS:PPPB2B_DATA_BITS)>endbit) {
        trace(2,"PPP-B2b body length error: bitpos=%d endbit=%d\n",
              pose,endbit);
        return -1;
    }
    if (prn_6==0||prn_6>63||mes_type==0||mes_type>63) {
        trace(2,"PPP-B2b header error: prn=%u type=%u body=%d\n",
              prn_6,mes_type,pose);
        return -1;
    }

    /* Build PRN(6)+status(6)+the standard 486-bit message after the complete
     * source frame, leaving the two CRC-protected source layers untouched. */
    dstbyte=raw->len+3;
    if (dstbyte+(12+PPPB2B_MESSAGE_BITS+7)/8>MAXRAWLEN) return -1;
    memset(raw->buff+dstbyte,0,(12+PPPB2B_MESSAGE_BITS+7)/8);
    dstbit=dstbyte*8;
    setbitu(raw->buff,dstbit,6,prn_6);
    setbitu(raw->buff,dstbit+6,6,status_6&0x3Fu);

    if (word_little_endian) {
        if ((pose&7)||pose/8+PPPB2B_PACKED_BYTES>raw->len) return -1;
        for (i=0;i<PPPB2B_MESSAGE_BITS;i++) {
            n=i;
            srcbyte=pose/8+(n/32)*4+3-(n%32)/8;
            srcbit=n&7;
            setbitu(raw->buff,dstbit+12+i,1,
                    getbitu(raw->buff,srcbyte*8+srcbit,1));
        }
        if (getbitu(raw->buff,dstbit+12,6)!=mes_type) {
            trace(2,"PPP-B2b message type mismatch: outer=%u inner=%u\n",
                  mes_type,getbitu(raw->buff,dstbit+12,6));
            return -1;
        }
        crc_rx=getbitu(raw->buff,dstbit+12+PPPB2B_CRC_INPUT_BITS,
                       PPPB2B_CRC_BITS);
        crc_calc=crc24q_bits(raw->buff,dstbit+12,PPPB2B_CRC_INPUT_BITS);
        if (crc_rx!=crc_calc) {
            trace(2,"PPP-B2b CRC24Q error: prn=%u type=%u recv=%06X calc=%06X\n",
                  prn_6,mes_type,crc_rx,crc_calc);
            return -1;
        }
        for (i=PPPB2B_MESSAGE_BITS;i<PPPB2B_PACKED_BITS;i++) {
            n=i;
            srcbyte=pose/8+(n/32)*4+3-(n%32)/8;
            srcbit=n&7;
            if (getbitu(raw->buff,srcbyte*8+srcbit,1)) {
                trace(2,"PPP-B2b non-zero packing bit: prn=%u type=%u bit=%d\n",
                      prn_6,mes_type,i-PPPB2B_MESSAGE_BITS);
                break;
            }
        }
    }
    else {
        setbitu(raw->buff,dstbit+12,6,mes_type);
        for (i=0;i<PPPB2B_DATA_BITS;i+=n) {
            n=PPPB2B_DATA_BITS-i<24?PPPB2B_DATA_BITS-i:24;
            setbitu(raw->buff,dstbit+18+i,n,getbitu(raw->buff,pose+i,n));
        }
    }
    if (mes_type>4) return 0; /* reserved/null message; standard CRC checked */
    message_end=word_little_endian?dstbit+12+PPPB2B_MESSAGE_BITS:
                                   dstbit+PPPB2B_INFO_BITS;
    return decode_B2b(raw,mask,dstbit,message_end);
}

/* Decode the compact SIS page: PRN(6), reserved(6), type(6), body(456). */
static int decode_PPPB2b_page(raw_t *raw, B2bmask_t *mask,
                              int bitpos, int endbit)
{
    uint32_t prn_6,status_6,mes_type;
    int pose=bitpos;

    if (!raw||bitpos<0||bitpos+PPPB2B_INFO_BITS>endbit) {
        trace(2,"PPP-B2b page length error: bitpos=%d endbit=%d\n",
              bitpos,endbit);
        return -1;
    }
    prn_6=getbitu(raw->buff,pose,6); pose+=6;
    status_6=getbitu(raw->buff,pose,6); pose+=6;
    mes_type=getbitu(raw->buff,pose,6); pose+=6;
    if (!mask) mask=get_rtcm_b2b_mask(prn_6);
    return decode_PPPB2b_body(raw,mask,prn_6,status_6,mes_type,pose,endbit,0);
}

/* decode PPP-B2b page carried by a KXW D3/CRC24Q frame ---------------------
 * The complete outer frame must already be present in raw->buff and raw->len
 * must be the byte count from the D3 preamble through the payload (excluding
 * the three CRC bytes), as set by input_SSR(). */
extern int decode_PPPB2b(raw_t *raw)
{
    int i,bitpos=-1,start=KXW_PPPB2B_PAGE_OFFSET;
    uint32_t rtcm_type,subtype,prn_6,mes_type;

    if (!raw||raw->len<=start) return -1;

    if (raw->time.time==0) {
        raw->time=utc2gpst(timeget());
        trace(2,"PPP-B2b: receiver time unavailable, using system time\n");
    }

    /* RTKNAVI can record the correction stream in a byte-aligned private
     * RTCM wrapper. Its payload is:
     *
     *   DF002=4047 (12), subtype=64 (12), version(8), PRN(8),
     *   status/reserved(8), message-type(8), then the complete 486-bit
     *   pre-LDPC SIS message packed in sixteen little-endian uint32_t words.
     *
     * The PRN/reserved/type fields are byte aligned and therefore must not be
     * interpreted as the compact 6/6/6-bit SIS header. */
    rtcm_type=getbitu(raw->buff,24,12);
    subtype=getbitu(raw->buff,36,12);
    if (rtcm_type==RTCM_PRIVATE_TYPE&&subtype==RTCM_PPPB2B_SUBTYPE&&
        raw->len>=RTCM_PPPB2B_BODY_OFFSET+PPPB2B_PACKED_BYTES) {
        prn_6=raw->buff[RTCM_PPPB2B_PRN_OFFSET];
        mes_type=raw->buff[RTCM_PPPB2B_TYPE_OFFSET];
        trace(3,"decode_PPPB2b: RTCM 4047/64 len=%d prn=%u type=%u\n",
              raw->len,prn_6,mes_type);
        return decode_PPPB2b_body(raw,get_rtcm_b2b_mask(prn_6),prn_6,
                                  raw->buff[RTCM_PPPB2B_RSV_OFFSET],mes_type,
                                  RTCM_PPPB2B_BODY_OFFSET*8,raw->len*8,1);
    }

    /* Prefer the standard 16-bit B2b preamble when it is present. */
    for (i=start;i+1<raw->len;i++) {
        if ((((uint16_t)raw->buff[i]<<8)|raw->buff[i+1])==PPPB2B_PREAMBLE) {
            bitpos=(i+2)*8;
            break;
        }
    }
    /* Some receivers output the decoded page without its preamble. */
    if (bitpos<0) bitpos=start*8;

    trace(3,"decode_PPPB2b: len=%d bitpos=%d preamble=%s\n",
          raw->len,bitpos,bitpos!=start*8?"yes":"no");

    return decode_PPPB2b_page(raw,NULL,bitpos,raw->len*8);
}

/* URA value (m) to URA index ------------------------------------------------*/
static int uraindex(double value)
{
    static const double ura_eph[] = {
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0,0.0
    };
    int i;
    for (i = 0; i < 15; i++) if (ura_eph[i] >= value) break;
    return i;
}

static int decode_GPSEPHEM(raw_t *raw) {
    PACKED_SINO_GPSEPHEM data_GPSEPH = {0};
    memcpy(&data_GPSEPH, raw->buff + SINAN_LEN, sizeof(data_GPSEPH));

    eph_t eph={0};
    double tow;
    int sat, week;

    if(!data_GPSEPH.blFlag){
        trace(2, "unicore gpsephb satellite error: prn=%d\n", data_GPSEPH.ID);
        return 0;
    }

    if (!(sat = satno(SYS_GPS, data_GPSEPH.ID))) {
        trace(2, "unicore gpsephb satellite error: prn=%d\n", data_GPSEPH.ID);
        return 0;
    }

    eph.sva = data_GPSEPH.bHealth;
    eph.sat = sat;
    eph.iodc = data_GPSEPH.iodc;
    eph.sva  = uraindex(data_GPSEPH.accuracy);  /* URA index (m->index) */

    eph.week=data_GPSEPH.week + 2048;      /* GPS week */
    eph.iode=data_GPSEPH.iode;
    eph.toes = data_GPSEPH.toe;

    eph.toe = gpst2time( eph.week,  eph.toes);
    double tt = timediff(eph.toe, raw->time);
    if (tt < -302400.0) eph.week++;
    else if (tt > 302400.0) eph.week--;

    eph.toe = gpst2time( eph.week,  eph.toes);
    eph.toc = gpst2time( eph.week,  data_GPSEPH.toc);
    eph.ttr = raw->time;

    eph.f0  = data_GPSEPH.af0;
    eph.f1  = data_GPSEPH.af1;
    eph.f2  = data_GPSEPH.af2;
    eph.M0  = data_GPSEPH.Ms0;
    eph.deln= data_GPSEPH.deltan;
    eph.e   = data_GPSEPH.es;
    eph.A   = data_GPSEPH.roora * data_GPSEPH.roora;
    eph.OMG0= data_GPSEPH.Omega0;
    eph.i0  = data_GPSEPH.i0;
    eph.omg = data_GPSEPH.ws;
    eph.OMGd= data_GPSEPH.omegaot;
    eph.idot= data_GPSEPH.itoet;
    eph.cuc = data_GPSEPH.cuc;
    eph.cus = data_GPSEPH.cus;
    eph.crc = data_GPSEPH.crc;
    eph.crs = data_GPSEPH.crs;
    eph.cic = data_GPSEPH.cic;
    eph.cis = data_GPSEPH.cis;
    eph.tgd[0] = data_GPSEPH.tgd;

    if (!strstr(raw->opt, "-EPHALL")) {
        if (fabs(timediff(raw->nav.eph[sat - 1].toe, eph.toe)) < 1e-9 &&
            fabs(timediff(raw->nav.eph[sat - 1].toc, eph.toc)) < 1e-9) return 0;
    }

    raw->nav.eph[sat - 1] = eph;
    raw->ephsat = sat;
    raw->ephset = 0;
    return 2;
}


static int decode_BD3EPHEM(raw_t *raw) {
    PACKED_SINO_BD3EPHEM data_BDSEPH = {0};
    memcpy(&data_BDSEPH, raw->buff + SINAN_LEN, sizeof(data_BDSEPH));

    eph_t eph = { 0 };
    gtime_t gps_toe = {0};
    int sat, sat_type, gps_week, gps_tow, ref_A;

    if (!(sat = satno(SYS_CMP, data_BDSEPH.Prn))) {
        trace(22, "unicore gpsephb satellite error: prn=%d\n", data_BDSEPH.Prn);
        return 0;
    }
    eph.sat = sat;
    sat_type = data_BDSEPH.sattype;
    eph.svh = data_BDSEPH.health;
    eph.sva = data_BDSEPH.URAI;
    eph.iode= data_BDSEPH.IODE;
    eph.iodc= data_BDSEPH.IODC;

    int sif = data_BDSEPH.SIF;
    if (sif != 0){ // sif is narmal
        trace(22,"The signal is incomplete! SIF = %d \n",sif);
        return 0;
    }

    int aif = data_BDSEPH.AIF;
    if (aif != 0){ // sif is narmal
        trace(22,"The SISMAI is effective! AIF = %d \n",aif);
        return 0;
    }
    gps_tow  = time2gpst(raw->time,&gps_week); // week con by ttr
    gps_toe  = gpst2time(gps_week,data_BDSEPH.toe); // TOE in gpst

    double tt = timediff(gps_toe, raw->time);
    if (tt < -302400.0) gps_week++;
    else if (tt > 302400.0) gps_week--;

    eph.week = gps_week - 1356; // GPS WEEK --->BDS WEEK
    eph.toes = data_BDSEPH.toe;
    eph.toe = bdt2gpst(bdt2time(eph.week,eph.toes));
    eph.toc = bdt2gpst(bdt2time(eph.week,data_BDSEPH.toc));
    eph.ttr = raw->time;

    // ref_A + delt_A = sqrt(A) * sqrt(A)
    if (sat_type ==1) ref_A  = 42162200; //(GEO)
    else if (sat_type ==2) ref_A  = 42162200; //(IGSO)
    else if (sat_type ==3) ref_A  = 27906100; //(MEO)

    eph.A = ref_A + data_BDSEPH.Delt_A;
    eph.Adot = data_BDSEPH.Dot_A;
    eph.deln = data_BDSEPH.Delt_n0;
    eph.ndot = data_BDSEPH.Dot_n0;
    eph.M0   = data_BDSEPH.M0;
    eph.e    = data_BDSEPH.e;
    eph.omg  = data_BDSEPH.w;
    eph.OMG0 = data_BDSEPH.Omega0;
    eph.i0   = data_BDSEPH.i0;
    eph.OMGd = data_BDSEPH.Omega_dot;
    eph.idot = data_BDSEPH.i_dot;

    eph.cuc = data_BDSEPH.Cuc;
    eph.cus = data_BDSEPH.Cus;
    eph.crc = data_BDSEPH.Crc;
    eph.crs = data_BDSEPH.Crs;
    eph.cic = data_BDSEPH.Cic;
    eph.cis = data_BDSEPH.Cis;

    eph.f0  = data_BDSEPH.a0;
    eph.f1  = data_BDSEPH.a1;
    eph.f2  = data_BDSEPH.a2;

    eph.tgd[2] = data_BDSEPH.tgdB1Cp;
    eph.tgd[4] = data_BDSEPH.tgdB1Cd;

    if (!strstr(raw->opt, "-EPHALL")) {
        if (fabs(timediff(raw->nav.eph[sat - 1].toe, eph.toe)) < 1e-9 &&
            fabs(timediff(raw->nav.eph[sat - 1].toc, eph.toc)) < 1e-9) return 0;
    }

    raw->nav.eph[sat - 1] = eph;
    raw->ephsat = sat;
    raw->ephset = 0;

    PPP_Glo.BDS_CNV1_flag = 1;

    return 2;
}



/* decode SINAN K803 LITE/W/S message -----------------------------------------*/
/*H28*8---prn32---prn6---reserve6---mes_type6---data966---CRC32*/
static int decode_sino(raw_t *raw)
{
    double tow;
    char tstr[32];
    int msg,stat,week,type=U2(raw->buff+4);
    int ret = 0;
    
    trace(3,"decode_sino: type=%3d len=%d\n",type,raw->len);
    
    /* check crc32 */
    if (rtk_crc32(raw->buff,raw->len)!=U4(raw->buff+raw->len)) {
        trace(2,"SINO crc error: type=%3d len=%d\n",type,raw->len);
        return -1;
    }
    // msg =(U1(raw->buff+6)>>4)&0x3; /* message type: 0=binary,1=ascii */
    // stat=U1(raw->buff+13);
    week=U2(raw->buff+14);
    
    if (stat==20||week==0) {
        trace(3,"SINO time error: type=%3d msg=%d stat=%d week=%d\n",type,msg,
              stat,week);
        return 0;
    }
    // week=adjgpsweek(week); 
    tow =U4(raw->buff+16)*0.001;
    raw->time=gpst2time(week,tow);
    char time_str[64];
    time2str(raw->time,time_str,5);
    // if (msg!=0) return 0;
    
    if (raw->outtype) {
        time2str(gpst2time(week,tow),tstr,2);
        sprintf(raw->msgtype,"SINO %4d (%4d): %s",type,raw->len,tstr);
    }


    switch (type) {
        case B2BRAWNAVSUBFRAMEb:
            /* Sinan payload: 28-byte header + 32-bit receiver PRN tag, then
             * the standard PPP-B2b page beginning at its 6-bit PRN. */
            ret=decode_PPPB2b_page(raw,&sinan_mask,SINAN_LEN*8+32,
                                   raw->len*8);
            break;
        case GPSEPHEM:           ret = decode_GPSEPHEM(raw); break;
        case BD3EPHEM:           ret = decode_BD3EPHEM(raw); break;

    }

    if (ret >= 0) {
        if(type == GPSEPHEM){
            raw->raw_nmsg[8]++;
        }
        else if(type == BD3EPHEM){
            raw->raw_nmsg[9]++;
        }
        
    }
    return ret;
}

/* input KXW D3/CRC24Q SSR stream -------------------------------------------*/
extern int input_SSR(raw_t *raw, uint8_t data)
{
    int payload_len,ret;

    if (!raw) return -1;
    trace(5,"input_SSR: data=%02x\n",data);

    /* synchronize frame */
    if (raw->nbyte==0) {
        if (data!=SSR3PREAMB) return 0;
        raw->buff[raw->nbyte++]=data;
        return 0;
    }
    if (raw->nbyte>=MAXRAWLEN) {
        trace(2,"SSR frame overflow\n");
        raw->nbyte=raw->len=0;
        return -1;
    }
    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte==3) {
        payload_len=getbitu(raw->buff,14,10);
        raw->len=payload_len+3; /* bytes excluding the 3-byte CRC */
        if (raw->len+3>MAXRAWLEN||raw->len<6) {
            trace(2,"SSR length error: payload=%d total=%d\n",
                  payload_len,raw->len+3);
            raw->nbyte=raw->len=0;
            return -1;
        }
    }
    if (raw->nbyte<3||raw->nbyte<raw->len+3) return 0;

    raw->nbyte=0;

    /* check outer-frame parity */
    if (rtk_crc24q(raw->buff,raw->len)!=
        getbitu(raw->buff,raw->len*8,24)) {
        trace(2,"SSR CRC24Q error: len=%d\n",raw->len);
        return -1;
    }
    if (raw->buff[5]==KXW_PPPB2B_MSG_ID) {
        /* decode_B2b returns its SIS message type (1..4), while RTKLIB raw
         * decoder return values 1..4 mean observation/ephemeris/SBAS/etc.
         * Normalize every decoded B2b page to the server's B2b update code. */
        ret=decode_PPPB2b(raw);
        return ret>0?20:ret;
    }
    trace(3,"SSR unsupported message id: 0x%02x\n",raw->buff[5]);
    return 0;
}


/*H28*8---prn32---prn6---reserve6---mes_type6---data966---CRC32*/
extern int input_sino(raw_t *raw, uint8_t data)
{
    // trace(5,"input_sinan_B2b: data=%02x\n",data);

    /* RTKNAVI records decoded PPP-B2b pages as D3/CRC24Q private RTCM
     * frames (4047/64). Keep accepting the native Sinan receiver protocol,
     * but dispatch a D3 frame to its dedicated parser before looking for the
     * native receiver sync word. */
    if ((raw->nbyte==0&&data==SSR3PREAMB)||
        (raw->nbyte>0&&raw->buff[0]==SSR3PREAMB)) {
        return input_SSR(raw,data);
    }
    
    /* synchronize frame */
    if (raw->nbyte==0) {
        if (sync_sinan(raw->buff,data)) raw->nbyte=3;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;
    
    // Get Message length
    if (raw->nbyte==10&&(raw->len=U2(raw->buff+8)+SINAN_LEN)>MAX_SINAN_LEN-4) {
        trace(2,"sinan length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (raw->nbyte<10||raw->nbyte<raw->len+4) return 0;
    raw->nbyte=0;
    
    /* decode sinan message */
    return decode_sino(raw);
}


extern int input_sinof(raw_t *raw, FILE *fp)
{
    int i, data = 0, ret;

    if (raw == NULL) {
        fprintf(stderr, "Error: raw is NULL\n");
        return -1;
    }
    if (fp == NULL) {
        fprintf(stderr, "Error: File pointer is NULL\n");
        return -1;
    }

    // trace(5,"input_sinof");

    for (i = 0; i < 4096; i++) {
        data = fgetc(fp);
        if (data == EOF) return -2;
        ret = input_sino(raw, (uint8_t)data);
        if (ret) return ret;
    }

    return 0;
}
