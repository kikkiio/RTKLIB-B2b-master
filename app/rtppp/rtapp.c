#include "rtapp.h"
#include "rtklib.h"
#include <signal.h>

PPPGlobal_t PPP_Glo = {0};

/* Bounded file-replay runs can shut down through app_rtkrcv's normal SIGTERM
 * handler, flushing all output streams instead of killing the process. */
#ifdef _WIN32
static DWORD WINAPI replay_timer(void *arg)
#else
static void *replay_timer(void *arg)
#endif
{
    sleepms((int)(intptr_t)arg*1000);
    raise(SIGTERM);
    return 0;
}

int main(int argc, char **argv)
{
    int i,j,seconds=0;
    for (i=1;i<argc;i++) {
        if (strcmp(argv[i],"--run-seconds")) continue;
        if (i+1>=argc||(seconds=atoi(argv[i+1]))<30||seconds>86400) {
            fprintf(stderr,"--run-seconds requires 30..86400 wall-clock seconds\n");
            return EXIT_FAILURE;
        }
        for (j=i;j+2<argc;j++) argv[j]=argv[j+2];
        argc-=2;argv[argc]=NULL;i--;
    }
    if (seconds) {
#ifdef _WIN32
        HANDLE timer=CreateThread(NULL,0,replay_timer,(void *)(intptr_t)seconds,0,NULL);
        if (!timer) return EXIT_FAILURE;
        CloseHandle(timer);
#else
        pthread_t timer;
        if (pthread_create(&timer,NULL,replay_timer,(void *)(intptr_t)seconds)) return EXIT_FAILURE;
        pthread_detach(timer);
#endif
    }
    //app_convbin(argc, argv);

    return app_rtkrcv(argc, argv);
}
