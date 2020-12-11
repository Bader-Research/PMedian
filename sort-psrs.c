#include "sort-psrs.h"
#include "sorting.h"

#define OVERSMP 1

#define DEBUG   0

void all_sort_psrs_check_i(int *A, int A_size) {
  int i;

  for (i=1 ; i<A_size ; i++)
    if (A[i] < A[i-1])
      fprintf(stderr,"(PE%3d)ERROR: A[%8d] < A[%8d]  (%12d %12d)\n",
	      MYPROC, i, i-1, A[i], A[i-1]);

  MPI_fprintf(stderr,"HERE CHECKED\n");
  
}

int *findpos (int *base, int size, int key) {
/*  find position of key (or next larger) in sorted array */

  int done = 0, *l, *u, *m, *end;

  l = base;
  u = end = base + size - 1;

  while (! done && size > 0) {
    m = l + size/2;
    if (key == *m)
      done = 1;
    else {
      size /= 2;
      if (key > *m)
	l = m+1;
      else
	u = m-1;
    }
  }
  
  while (m < end && key > *m)
    ++m;

  return (m);
}


int all_sort_psrs_i(int *A, int A_size, int *B) {


  int i, B_size;

  int 
      *recv_cnt,     /* number of keys to receive from each PN */
      *send_cnt,     /* number of keys to send to PNs */
      *recv_off,
      *send_off;

  int no_samples,    /* number of samples to take from sorted data */
      *pivbuffer,    /* array of pivots */
      pivdist,       /* distance between pivots in set of samples */
      possmp,        /* position of pivot in set of all samples */
      **prtbound,    /* boundaries for partitioning local data */
      *recv_buf,     /* incoming data to merge */
      *smpbuffer,    /* array of local samples */
      *smpallbuf,    /* array of global samples */
      smpdist,       /* distance between consecutive samples */
      trecv;

  MPI_Status stat;

  prtbound = (int **) malloc((PROCS+1)*sizeof(int *));
  assert_malloc(prtbound);
  
  send_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(send_cnt);
  
  send_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(send_off);
  
  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);
  
  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  pivbuffer = (int *)malloc((PROCS-1)*sizeof(int));
  assert_malloc(pivbuffer);
  
  /*------------------------------------------------------------
   -             Phase 2: Sort
   *------------------------------------------------------------*/

  fastsort(A, A_size);

  /*------------------------------------------------------------
   -             Phase 3: Find pivots
   *------------------------------------------------------------*/
  /*
   -  Sample the sorted array
   */
  no_samples = OVERSMP*PROCS - 1;
  smpdist = A_size / (no_samples + 1);

  smpbuffer = (int *)malloc(no_samples*sizeof(int));
  assert_malloc(smpbuffer);

  for (i=0; i<no_samples; i++)
    smpbuffer[i] = A[(i+1)*smpdist];

#if DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0 ; i<no_samples ; i++)
    MPI_fprintf(stderr,"(s %5d) smp[%5d]: %12d\n",
		no_samples, i, smpbuffer[i]);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  /*
   -  Concatenate samples, sort and select pivots
   */
  smpallbuf = (int *)malloc(PROCS*no_samples*sizeof(int));
  assert_malloc(smpallbuf);

  MPI_Gather(smpbuffer, no_samples, MPI_INT,
	     smpallbuf, no_samples, MPI_INT,
	     0, MPI_COMM_WORLD);


#if 0
  MPI_Barrier(MPI_COMM_WORLD);
  if (MYPROC==0)
    for (i=0 ; i<no_samples*PROCS ; i++)
      fprintf(stderr,"(PE%3d) sample[%5d]: %12d\n",
	      MYPROC, i, smpallbuf[i]);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  fastsort(smpallbuf, PROCS*no_samples);

#if 0
  MPI_Barrier(MPI_COMM_WORLD);
  if (MYPROC==0)
    for (i=0 ; i<no_samples*PROCS ; i++)
      fprintf(stderr,"(PE%3d) SAMPLE[%5d]: %12d\n",
	      MYPROC, i, smpallbuf[i]);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (MYPROC==0) {
    pivdist = PROCS*no_samples/(PROCS-1);
    possmp = pivdist/2;
    for (i=0; i<PROCS-1; i++) {
      pivbuffer[i] = smpallbuf[possmp];
      possmp += pivdist;
    }
  }

  MPI_Bcast(pivbuffer, PROCS-1, MPI_INT, 0, MPI_COMM_WORLD);
  
#if DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  if (MYPROC==0)
    fprintf(stderr,"pivdist: %12d\n",pivdist);
  for (i=0 ; i<PROCS-1 ; i++)
    MPI_fprintf(stderr,"pivbuffer[%5d]: %12d\n", i, pivbuffer[i]);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  /*
   -  Search for pivots in local buffer
   */
  prtbound[0] = A;
  for (i=1; i<PROCS; i++)
    prtbound[i] = findpos (A, A_size, pivbuffer[i-1]);

  prtbound[PROCS] = A + A_size;

#if DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0 ; i<=PROCS ; i++)
    MPI_fprintf(stderr,"prtbound[%3d]: %5d\n",i,prtbound[i]-A);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  free (smpbuffer);
  free (smpallbuf);
  free (pivbuffer);

  /*------------------------------------------------------------
   -            Phase 4: Redistribution
   *------------------------------------------------------------*/

  /*
   -  Find out how many keys have to be sent to each PN
   */

  for (i=0; i<PROCS; i++)
    send_cnt[i] = max(0, prtbound[i+1] - 1 - prtbound[i]);

  /*
   -  Find out how many keys have to be received from each PN
   -
   -	 communication pattern: at iteration i,
   -  	     send send_cnt to PN distant i ahead and 
   -   	     receive recv_cnt from PN distant i behind
   -
   */

#if DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0 ; i<PROCS ; i++)
    MPI_fprintf(stderr,"send_cnt[%3d]: %5d\n",i,send_cnt[i]);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  MPI_Alltoall(send_cnt, 1, MPI_INT,
	       recv_cnt, 1, MPI_INT,
	       MPI_COMM_WORLD);
   	

#if DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0 ; i<PROCS ; i++)
    MPI_fprintf(stderr,"recv_cnt[%3d]: %5d\n",i,recv_cnt[i]);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  /*
   -  allocate space for incoming data 
   */

  send_off[0] = 0;
  recv_off[0] = 0;
  for (i=1 ; i<PROCS ; i++) {
    send_off[i] = send_off[i-1] + send_cnt[i-1];
    recv_off[i] = recv_off[i-1] + recv_cnt[i-1];
  }

  trecv = recv_off[PROCS-1] + recv_cnt[PROCS-1];

  if (trecv>0) {
    recv_buf = (int *)malloc(trecv*sizeof(int));
    assert_malloc(recv_buf);
  }
  else {
    recv_buf = NULL;
  }

  /* 
   -  send and receive keys
   */
  MPI_Alltoallv(A,        send_cnt, send_off, MPI_INT,
		recv_buf, recv_cnt, recv_off, MPI_INT,
		MPI_COMM_WORLD);
  

  /*-----------------------------------------------------------
   *                  Phase 5: Merge received keys
   *-----------------------------------------------------------*/

  /*
   -  in this version, simply sort again
   -  (to be substituted by an p-way merge function)
   */
  B_size = trecv;
  
  bcopy(recv_buf, B, B_size*sizeof(int));

  fastsort(B, B_size);

  MPI_Barrier(MPI_COMM_WORLD);


  free(recv_buf);
  free(recv_off);
  free(send_off);
  free(send_cnt);
  free(recv_cnt);
  free(prtbound);
  return (B_size);
}
