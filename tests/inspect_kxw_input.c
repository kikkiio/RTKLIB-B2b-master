/* Read-only audit of a tagged KXW correction file using the production decoder. */
#include "rtklib.h"
PPPGlobal_t PPP_Glo={0};
int main(int argc,char **argv)
{
    raw_t *r;stream_t stream;FILE *fp;char path[MAXSTRPATH];
    int ret,errors=0,nav=0,corrections=0,i;
    if (argc!=2) {fprintf(stderr,"usage: inspect_kxw_input PPP_file (requires .tag)\n");return 1;}
    r=calloc(1,sizeof(*r));if(!r||!init_raw(r,STRFMT_KXW))return 1;
    snprintf(path,sizeof(path),"%s::T",argv[1]);strinit(&stream);
    if(!stropen(&stream,STR_FILE,STR_MODE_R,path))return 1;
    r->time=strgettime(&stream);strclose(&stream);
    if(!(fp=fopen(argv[1],"rb")))return 1;
    while((ret=input_rawf(r,STRFMT_KXW,fp))!=-2) {
        if(ret==-1)errors++;
        if(ret==2)nav++;
        if(ret==20)corrections++;
        for(i=0;i<MAXSAT;i++)r->nav.B2bssr[i].update=0;
    }
    printf("{\"decode_errors\":%d,\"accepted_cnav1\":%d,\"correction_events\":%d,",errors,nav,corrections);
    printf("\"message_counts\":{\"mask\":%d,\"orbit\":%d,\"dcb\":%d,\"clock\":%d,\"other\":%d,\"cnav1\":%d},",
           r->raw_nmsg[0],r->raw_nmsg[1],r->raw_nmsg[2],r->raw_nmsg[3],r->raw_nmsg[7],r->raw_nmsg[9]);
    printf("\"last_correction_time_gpst\":\"%s\"}\n",time_str(r->time,3));
    fclose(fp);free_raw(r);free(r);return errors?1:0;
}
