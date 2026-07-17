/*------------------------------------------------------------------------------
* enum.c : Enumeration type conversion functions
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "enum.h"

/* Convert E_NavRecType to string ------------------------------------------*/
const char* NavRecToString(E_NavRecType type) {
    switch (type) {
        case E_NavRec_NONE: return "NONE";
        case EPH: return "EPH";
        case STO: return "STO";
        case EOP: return "EOP";
        case ION: return "ION";
        default: return "Unknown";
    }
}

/* Convert E_NavMsgType to string ------------------------------------------*/
const char* NavMsgToString(E_NavMsgType type) {
    switch (type) {
        case E_NavMsg_NONE: return "NONE";
        case LNAV: return "LNAV";
        case FDMA: return "FDMA";
        case FNAV: return "FNAV";
        case INAV: return "INAV";
        case IFNV: return "IFNV";
        case D1: return "D1";
        case D2: return "D2";
        case D1D2: return "D1D2";
        case SBAS: return "SBAS";
        case CNAV: return "CNAV";
        case CNV1: return "CNV1";
        case CNV2: return "CNV2";
        case CNV3: return "CNV3";
        case CNVX: return "CNVX";
        default: return "Unknown";
    }
}

/* Convert string to E_NavRecType ------------------------------------------*/
E_NavRecType StringToNavRec(const char* str) {
    if (strcmp(str, "NONE") == 0) return E_NavRec_NONE;
    if (strcmp(str, "EPH") == 0) return EPH;
    if (strcmp(str, "STO") == 0) return STO;
    if (strcmp(str, "EOP") == 0) return EOP;
    if (strcmp(str, "ION") == 0) return ION;
    return E_NavRec_NONE;  // Default
}

/* Convert string to E_NavMsgType ------------------------------------------*/
E_NavMsgType StringToNavMsg(const char* str) {
    if (strcmp(str, "NONE") == 0) return E_NavMsg_NONE;
    if (strcmp(str, "LNAV") == 0) return LNAV;
    if (strcmp(str, "FDMA") == 0) return FDMA;
    if (strcmp(str, "FNAV") == 0) return FNAV;
    if (strcmp(str, "INAV") == 0) return INAV;
    if (strcmp(str, "IFNV") == 0) return IFNV;
    if (strcmp(str, "D1  ") == 0) return D1;
    if (strcmp(str, "D2  ") == 0) return D2;
    if (strcmp(str, "D1D2") == 0) return D1D2;
    if (strcmp(str, "SBAS") == 0) return SBAS;
    if (strcmp(str, "CNAV") == 0) return CNAV;
    if (strcmp(str, "CNV1") == 0) return CNV1;
    if (strcmp(str, "CNV2") == 0) return CNV2;
    if (strcmp(str, "CNV3") == 0) return CNV3;
    if (strcmp(str, "CNVX") == 0) return CNVX;
    return E_NavMsg_NONE;  // Default
}

/* Default navigation message type for RINEX 3 and 2 ----------------------*/
const E_NavMsgType defNavMsgType[] = {
    [SYS_GPS] = LNAV,
    [SYS_GLO] = FDMA,
    [SYS_GAL] = IFNV,
    [SYS_CMP] = D1D2,
    [SYS_QZS] = LNAV,
    [SYS_IRN] = LNAV,
    [SYS_SBS] = SBAS
};
