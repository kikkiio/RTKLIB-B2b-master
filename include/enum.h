// 定义枚举类型 E_NavRecType
#ifndef ENUM_H
#define ENUM_H

#ifdef __cplusplus
extern "C" {
#endif

//rinex4 包含的数据类型
typedef enum {
    E_NavRec_NONE,			///< NONE for unknown */
    EPH,			///< Ephemerides data including orbit, clock, biases, accuracy and status parameters */
    STO,			///< System Time and UTC proxy offset parameters */
    EOP,			///< Earth Orientation Parameters */
    ION			///< Global/Regional ionospheric model parameters */
} E_NavRecType;

//星历类型
typedef enum {
    E_NavMsg_NONE,			///< NONE for unknown
    LNAV,			///< GPS/QZSS/NavIC Legacy Navigation Messages
    FDMA,			///< GLONASS Legacy FDMA Navigation Message
    FNAV,			///< Galileo Free Navigation Message
    INAV,			///< Galileo Integrity Navigation Message
    IFNV,			///< Galileo INAV or FNAV Navigation Message
    D1,				///< BeiDou-2/3 MEO/IGSO Navigation Message
    D2,				///< BeiDou-2/3 GEO Navigation Message
    D1D2,			///< BeiDou-2/3 MEO/IGSO and GEO Navigation Message
    SBAS,			///< SBAS Navigation Message
    CNAV,			///< GPS/QZSS CNAV Navigation Message
    CNV1,			///< BeiDou-3 CNAV-1 Navigation Message
    CNV2,			///< GPS/QZSS CNAV-2 Navigation Message	  BeiDou-3 CNAV-2 Navigation Message
    CNV3,			///< BeiDou-3 CNAV-3 Navigation Message
    CNVX			///< GPS/QZSS CNAV or CNAV-2 Navigation Message  BeiDou-3 CNAV-1, CNAV-2 or CNAV-3 Navigation Message
} E_NavMsgType;

//星历种类
typedef enum {
    E_Eph_NONE,			///< NONE for unknown
    Eph_EPH,			///< GPS/QZS LNAV, GAL IFNV, BDS D1D2 Ephemeris
    Eph_GEPH,			///< GLO Ephemeris
    Eph_SEPH,			///< SBAS Ephemeris
    Eph_CEPH,			///< GPS/QZS/BDS CNVX Ephemeris
    // Eph_STO,			///< STO message
    // Eph_EOP,			///< EOP message
    // Eph_ION			///< ION message
}E_EphType;

// 函数声明
extern const char* NavRecToString(E_NavRecType type);
extern const char* NavMsgToString(E_NavMsgType type);
extern E_NavRecType StringToNavRec(const char* str);
extern E_NavMsgType StringToNavMsg(const char* str);

// 声明全局数组
extern const E_NavMsgType defNavMsgType[];

#ifdef __cplusplus
}
#endif

#endif /* ENUM_H */