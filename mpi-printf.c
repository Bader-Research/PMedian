#include "mpi-printf.h"
#ifdef SUN
#include <stdarg.h>
#endif

void MPI_fprintf(FILE *fp, char *fmt, ...)
{
#define BUFSIZE 128
    int rank;
    int tag=10, p,np;
    va_list args;
    char str[BUFSIZE];
    MPI_Status status;
#ifdef OLD_SUN
    va_start(args);
#else
    va_start(args,fmt);
#endif

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (rank==0) {
        fprintf(fp,"(PE%3d): ",rank);
        vfprintf(fp,fmt,args);
        for (p=1; p<np; p++) {
            MPI_Recv(str,BUFSIZE,MPI_CHAR,p,tag, MPI_COMM_WORLD, &status);
            fprintf(fp,"(PE%3d): ",p);
            fprintf(fp,"%s",str);
        }
    } else {
        vsprintf(str,fmt,args);
        MPI_Send(str,BUFSIZE,MPI_CHAR,0,tag,MPI_COMM_WORLD);
    }
    va_end(args);
}

void MPI_printf(char *fmt, ...)
{
#if 0
  
  va_list args;
#ifdef OLD_SUN
  va_start(args);
#else
  va_start(args,fmt);
#endif
  MPI_fprintf(stdout,fmt,args);
  va_end(args);

#else

#define BUFSIZE 128
    int rank;
    int tag=10, p,np;
    va_list args;
    char str[BUFSIZE];
    MPI_Status status;
    FILE *fp;
#ifdef OLD_SUN
    va_start(args);
#else
    va_start(args,fmt);
#endif

    fp = stdout;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (rank==0) {
        fprintf(fp,"(PE%3d): ",rank);
        vfprintf(fp,fmt,args);
        for (p=1; p<np; p++) {
            MPI_Recv(str,BUFSIZE,MPI_CHAR,p,tag, MPI_COMM_WORLD, &status);
            fprintf(fp,"(PE%3d): ",p);
            fprintf(fp,"%s",str);
        }
    } else {
        vsprintf(str,fmt,args);
        MPI_Send(str,BUFSIZE,MPI_CHAR,0,tag,MPI_COMM_WORLD);
    }
    va_end(args);
#endif

}


