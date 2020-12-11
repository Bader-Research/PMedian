#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "main.h"

int MYPROC, PROCS;
int PROCSLOG;

main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &PROCS);
  MPI_Comm_rank(MPI_COMM_WORLD, &MYPROC);

  outfile = stdout;
  
  fprintf(outfile,"PE%3d: PROCS: %d\n",  MYPROC, PROCS);
  fflush(outfile);

  if (PROCS > MAXPROCS) {
    if (MYPROC==0) fprintf(stderr,"PE%3d: ERROR. MAXPROCS: %d PROCS: %d\n",
		   MYPROC, MAXPROCS, PROCS);
    sleep(1);
    exit(-1);
  }

  PROCSLOG = (int)rint(log2((double)PROCS));
  MPI_Barrier(MPI_COMM_WORLD);

  srrandom(1001*MYPROC + 21);
  
  all_test_select();
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
