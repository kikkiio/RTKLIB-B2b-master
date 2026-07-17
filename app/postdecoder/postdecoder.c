/*------------------------------------------------------------------------------
* postdecoder.c : test decoding for various receiver-specific custom data formats
*
* Built by Liu@APM(SDUST), Jiang@SEU(SDUST), Li@SEU(SDUST), Zhu@SEU(SDUST) for RTCM3 decoder
*
* version : $Revision:$ $Date:$
* history : 14 Sep 2023  1.0 created file for initial decoding tests by Liu@APM, Jiang@SEU, Li@SEU, Zhu@SEU
*           05 Oct 2023  1.1 added RTCM3 SSR format decoding support by Liu@APM
*           19 Nov 2024  1.2 added Sinan/Unicore B2b custom format decoding by Liu@APM,Xu@CUGB
*
* description : This program tests decoding for receiver manufacturers' custom data formats.
*               Input design simulates the RTKLIB framework for enhanced compatibility.
*-----------------------------------------------------------------------------*/
#include "postdecoder.h"
#include "rtklib.h"
#include <ctype.h>

#define TRACEFILE   "./B2bPPP_%Y_%m_%d.trace"    /* debug trace file */
#define B2bTRACEFILE "./B2bPPP_%Y_%m_%d.B2bssr"   /* B2b SSR debug file */

// static nav_t navs={0};          /* navigation data structure */


PPPGlobal_t PPP_Glo = {0};        /* PPP global parameters */

/* detect whitespace-separated hexadecimal text without consuming input -----*/
static int is_hex_text_file(FILE *fp)
{
    long pos;
    int c,n=0,digits=0,is_hex=1;

    if (!fp||(pos=ftell(fp))<0) return 0;
    while (n++<512&&(c=fgetc(fp))!=EOF) {
        if (isspace((unsigned char)c)) continue;
        if (!isxdigit((unsigned char)c)) {is_hex=0; break;}
        digits++;
    }
    fseek(fp,pos,SEEK_SET);
    return is_hex&&digits>=2;
}

static int hex_value(int c)
{
    if (c>='0'&&c<='9') return c-'0';
    if (c>='a'&&c<='f') return c-'a'+10;
    if (c>='A'&&c<='F') return c-'A'+10;
    return -1;
}

/* input binary or whitespace-separated hexadecimal PPP-B2b stream ---------*/
static int input_B2bf(raw_t *raw, FILE *fp, int hex_text)
{
    int i=0,c,value,high=-1,ret;

    if (!raw || !fp) return -1;

    while (i<4096) {
        if ((c=fgetc(fp))==EOF) {
            if (high>=0) {
                printf("Error: odd number of hexadecimal digits in input.\n");
                return -3;
            }
            return -2;
        }
        if (hex_text) {
            if (isspace((unsigned char)c)) continue;
            if ((value=hex_value(c))<0) {
                printf("Error: invalid hexadecimal character 0x%02X.\n",
                       (unsigned char)c);
                return -3;
            }
            if (high<0) {high=value; continue;}
            c=(high<<4)|value;
            high=-1;
        }
        i++;
        if ((ret=input_SSR(raw,(uint8_t)c))!=0) return ret;
    }
    return 0;
}

/* set a BDT calendar-day reference for pages carrying seconds-of-day only --*/
static int set_B2b_date(raw_t *raw, const char *date)
{
    double ep[6] = {0};
    int year, month, day;

    if (!raw || !date ||
        sscanf(date, "%d/%d/%d", &year, &month, &day) != 3 ||
        year < 1980 || month < 1 || month > 12 || day < 1 || day > 31) {
        return 0;
    }
    ep[0] = year;
    ep[1] = month;
    ep[2] = day;
    ep[3] = 12.0; /* noon keeps every page TOD within the +/-12 h adjustment */

    raw->time = bdt2gpst(epoch2time(ep));
    return 1;
}

static void print_usage(void)
{
    printf("Usage: postdecoder [-U|-S|-B] -in B2b_filepath "
           "[-out B2bSSR_resultpath] [-date YYYY/MM/DD]\n");
    printf("  -U  Unicore receiver binary format\n");
    printf("  -S  SinoGNSS receiver binary format\n");
    printf("  -B  D3/CRC24Q PPP-B2b stream (compact or RTKNAVI RTCM 4047/64)\n");
    printf("  -date  BDT calendar date for -B pages carrying seconds-of-day only\n");
}


int main(int argc, char **argv)
{
    int i, ret;
    char output_path[1024];        // Buffer for output file path
    const char *B2bfilepath = NULL; // Input file path
    const char *format_str = NULL;  // Format type identifier
    const char *output_specified = NULL; // User-specified output path
    const char *date_specified = NULL; /* BDT date for standard raw pages */
    int standard_B2b = 0;
    int hex_text = 0;

    // Parse command-line arguments
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-U") == 0) {
            format_str = "U";
        } else if (strcmp(argv[i], "-S") == 0) {
            format_str = "S";
        } else if (strcmp(argv[i], "-B") == 0) {
            format_str = "B";
        } else if (strcmp(argv[i], "-in") == 0 && i + 1 < argc) {
            B2bfilepath = argv[++i];
        } else if (strcmp(argv[i], "-out") == 0 && i + 1 < argc) {
            output_specified = argv[++i];
        } else if (strcmp(argv[i], "-date") == 0 && i + 1 < argc) {
            date_specified = argv[++i];
        } else {
            printf("Invalid option: %s\n", argv[i]);
            print_usage();
            return -1;
        }
    }

    // Validate required parameters
    if (format_str == NULL || B2bfilepath == NULL) {
        printf("Missing required options.\n");
        print_usage();
        return -1;
    }

    // Determine data format
    int format;
    if (format_str[0] == 'S' || format_str[0] == 's') {
        format = STRFMT_SINO;
    } else if (format_str[0] == 'U' || format_str[0] == 'u') {
        format = STRFMT_UNICORE;
    } else if (format_str[0] == 'B' || format_str[0] == 'b') {
        /* input_SSR() is receiver-independent; use an initialized raw_t. */
        format = STRFMT_SINO;
        standard_B2b = 1;
    } else {
        printf("Invalid format type.\n");
        print_usage();
        return -1;
    }

    // Configure output path
    if (output_specified != NULL) {
        if (snprintf(output_path, sizeof(output_path), "%s",
                     output_specified) >= (int)sizeof(output_path)) {
            printf("Output path is too long.\n");
            return -1;
        }
    } else if (snprintf(output_path, sizeof(output_path), "%s.B2bSSR",
                        B2bfilepath) >= (int)sizeof(output_path)) {
        printf("Output path is too long.\n");
        return -1;
    }

    /* These structures are too large for the default Windows thread stack. */
    static raw_t raw0;           // Raw data container
    static nav_t navs;           // Navigation data container

    if (init_raw(&raw0, format) != 1) {
        printf("Error initializing raw data handler!\n");
        return -1;
    }

    if (standard_B2b && date_specified &&
        !set_B2b_date(&raw0, date_specified)) {
        printf("Invalid -date value: %s (expected YYYY/MM/DD)\n",
               date_specified);
        free_raw(&raw0);
        return -1;
    }

    FILE *B2bfp = fopen(B2bfilepath, "rb");
    if (B2bfp == NULL) {
        printf("Error: Failed to open: %s\n", B2bfilepath);
        free_raw(&raw0);
        return -1;
    }
    if (standard_B2b) hex_text=is_hex_text_file(B2bfp);

    // tracelevel(100);             // Enable full tracing
    B2b_tracelevel(22);
    B2b_traceopen(output_path);  // Initialize tracing with output path

    printf("Decoding started (%s)...\n",
           standard_B2b ? (hex_text ? "PPP-B2b hexadecimal text" :
                                      "PPP-B2b D3/RTKNAVI stream") :
           format == STRFMT_UNICORE ? "Unicore" : "SinoGNSS");
    while (1) {
        ret = standard_B2b ? input_B2bf(&raw0,B2bfp,hex_text) :
                             input_rawf(&raw0, format, B2bfp);
        if (ret < -1) break;
        for (i = 0; i < MAXSAT; i++) {
            if (!raw0.nav.B2bssr[i].update) continue;
            if (!standard_B2b && raw0.num_PPPB2BINF02 != 0) {
                raw0.nav.B2bssr[i].udi[0] = timediff(raw0.nav.B2bssr[i].t0[0], navs.B2bssr[i].t0[0]);
                if (raw0.nav.B2bssr[i].udi[0] > 86400) raw0.nav.B2bssr[i].udi[0] = 0;
            }
            if (!standard_B2b && raw0.num_PPPB2BINF03 != 0) {
                raw0.nav.B2bssr[i].udi[1] = timediff(raw0.nav.B2bssr[i].t0[1], navs.B2bssr[i].t0[1]);
                if (raw0.nav.B2bssr[i].udi[1] > 86400) raw0.nav.B2bssr[i].udi[1] = 0;
            }
            if (!standard_B2b && raw0.num_PPPB2BINF04 != 0) {
                raw0.nav.B2bssr[i].udi[2] = timediff(raw0.nav.B2bssr[i].t0[2], navs.B2bssr[i].t0[2]);
                if (raw0.nav.B2bssr[i].udi[2] > 86400) raw0.nav.B2bssr[i].udi[2] = 0;
            }
            navs.B2bssr[i] = raw0.nav.B2bssr[i];
            raw0.nav.B2bssr[i].update = 0;
        }
        raw0.num_PPPB2BINF01 = 0;
        raw0.num_PPPB2BINF02 = 0;
        raw0.num_PPPB2BINF03 = 0;
        raw0.num_PPPB2BINF04 = 0;
    }

    traceclose();                // Close trace files
    B2b_traceclose();
    fclose(B2bfp);

    printf("Decoding completed: type1=%u type2=%u type3=%u type4=%u\n",
           raw0.raw_nmsg[0], raw0.raw_nmsg[1], raw0.raw_nmsg[2],
           raw0.raw_nmsg[3]);
    free_raw(&raw0);

    return 0;
}
