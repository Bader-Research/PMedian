/*                                                                    tab:8
 *
 * select_par.c - Parallel Selection Algorithm
 *
 * 
 * "Copyright (c) 1996 The Regents of the University of Maryland.
 * All rights reserved.
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose, without fee, and without written agreement is
 * hereby granted, provided that the above copyright notice and the following
 * two paragraphs appear in all copies of this software.
 * 
 * IN NO EVENT SHALL THE UNIVERSITY OF MARYLAND BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
 * OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
 * MARYLAND HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE UNIVERSITY OF MARYLAND SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
 * ON AN "AS IS" BASIS, AND THE UNIVERSITY OF MARYLAND HAS NO OBLIGATION TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS."
 *
 * Authors:             David A. Bader   <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:             1.0
 * Creation Date:       February 6, 1996
 * Filename:            select_par.c
 * History:
 */

#include "select_par.h"
#include "select_seq.h"

#define DEBUG         0

#define BUNDLE        4

#if USE_MY_ALLTOALLV
#define MPIR_ALLTOALLV_TAG 10
int init_Alltoallv_param = 0;
MPI_Request *reqarray;
MPI_Status  *starray;

void init_Alltoallv() {
  reqarray = (MPI_Request *)malloc(PROCS*sizeof(MPI_Request));
  assert_malloc(reqarray);
  starray  = (MPI_Status  *)malloc(PROCS*sizeof(MPI_Status));
  assert_malloc(starray);
  init_Alltoallv_param = 1;
}

void destroy_Alltoallv() {
  free(starray);
  free(reqarray);
  init_Alltoallv_param = 0;
}
#endif

void all_randomize_i(int *A, int *B, int M, int bin_size, int bins,
		     int logbins, int *bin_counts) {

  int *A_ptr;
  int *B_ptr;
  int *ptr;
  int *record_ptr;
  int **Loc;
  
  int
    v, k,
    i, j,
    bins_m1,
    times,
    shift,
    result,
    rem,
    maxsz;

  Loc = (int **)malloc(bins*sizeof(int *));
  assert_malloc(Loc);
  
  srandom(1001*MYPROC + 21);

  bins_m1 = bins - 1;
  
  A_ptr = A;
  B_ptr = B - bin_size - BUNDLE;
  
  
  for (i=0 ; i<bins ; i++)
    Loc[i] = (B_ptr += bin_size);

#if BUNDLE == 1

  times  = 31/logbins;
  shift  = 31 - times*logbins;
  result = M/times;
  rem    = M - times*result;
  for (i=0 ; i<result ; i++) {
    v = rrandom() >> shift;
    ptr = Loc[v & bins_m1]++;
    *ptr = *(A_ptr++);
    for (j=1 ; j<times ; j++) {
      v >>= logbins;
      ptr = Loc[v & bins_m1]++;
      *ptr = *(A_ptr++);
    }
  }

  v = rrandom();
  for (i=0 ; i<rem ; i++) {
    v >>= logbins;
    ptr = Loc[v & bins_m1]++;
    *ptr = *(A_ptr++);
  }

#endif

#if BUNDLE > 1

  times  = 31/logbins;
  shift  = 31 - times*logbins;
  result = M/(times*BUNDLE);
  rem    = (M - BUNDLE*times*result)/BUNDLE;

  for (i=0 ; i<result ; i++) {
    v = rrandom() >> shift;
    record_ptr = (Loc[v & bins_m1] += BUNDLE); 
    for (k=0 ; k<BUNDLE ; k++)
     record_ptr[k] = *(A_ptr++);
    for (j=1 ; j<times ; j++) {
      v >>= logbins;
      record_ptr = (Loc[v & bins_m1] += BUNDLE);
      for (k=0 ; k<BUNDLE ; k++)
        record_ptr[k] = *(A_ptr++);
    }
  }

  v = rrandom();
  for (i=0 ; i<rem ; i++) {
    v >>= logbins;
    record_ptr = (Loc[v & bins_m1] += BUNDLE);
    for (k=0 ; k<BUNDLE ; k++)
      record_ptr[k] = *(A_ptr++);
  }
  
#endif

  B_ptr = B;
  maxsz = *bin_counts = *Loc + BUNDLE - B_ptr; 
  for (i=1 ; i<bins ; i++) {
    bin_counts[i] = Loc[i] + BUNDLE - (B_ptr+=bin_size);
    maxsz = max(bin_counts[i], maxsz);
  }

  if ((maxsz-2) > bin_size) {
    fprintf(stderr,"ERROR: bins too small\n");
    exit(1);
  }

  free(Loc);
}


void all_randomize_d(double *A, double *B, int M, int bin_size, int bins,
		     int logbins, int *bin_counts) {

  double *A_ptr;
  double *B_ptr;
  double *ptr;
  double *record_ptr;
  double **Loc;
  
  int
    v, k,
    i, j,
    bins_m1,
    times,
    shift,
    result,
    rem,
    maxsz;

  Loc = (double **)malloc(bins*sizeof(double *));
  assert_malloc(Loc);
  
  srandom(1001*MYPROC + 21);

  bins_m1 = bins - 1;
  
  A_ptr = A;
  B_ptr = B - bin_size - BUNDLE;
  
  
  for (i=0 ; i<bins ; i++)
    Loc[i] = (B_ptr += bin_size);

#if BUNDLE == 1

  times  = 31/logbins;
  shift  = 31 - times*logbins;
  result = M/times;
  rem    = M - times*result;
  for (i=0 ; i<result ; i++) {
    v = rrandom() >> shift;
    ptr = Loc[v & bins_m1]++;
    *ptr = *(A_ptr++);
    for (j=1 ; j<times ; j++) {
      v >>= logbins;
      ptr = Loc[v & bins_m1]++;
      *ptr = *(A_ptr++);
    }
  }

  v = rrandom();
  for (i=0 ; i<rem ; i++) {
    v >>= logbins;
    ptr = Loc[v & bins_m1]++;
    *ptr = *(A_ptr++);
  }

#endif

#if BUNDLE > 1

  times  = 31/logbins;
  shift  = 31 - times*logbins;
  result = M/(times*BUNDLE);
  rem    = (M - BUNDLE*times*result)/BUNDLE;

  for (i=0 ; i<result ; i++) {
    v = rrandom() >> shift;
    record_ptr = (Loc[v & bins_m1] += BUNDLE); 
    for (k=0 ; k<BUNDLE ; k++)
     record_ptr[k] = *(A_ptr++);
    for (j=1 ; j<times ; j++) {
      v >>= logbins;
      record_ptr = (Loc[v & bins_m1] += BUNDLE);
      for (k=0 ; k<BUNDLE ; k++)
        record_ptr[k] = *(A_ptr++);
    }
  }

  v = rrandom();
  for (i=0 ; i<rem ; i++) {
    v >>= logbins;
    record_ptr = (Loc[v & bins_m1] += BUNDLE);
    for (k=0 ; k<BUNDLE ; k++)
      record_ptr[k] = *(A_ptr++);
  }
  
#endif

  B_ptr = B;
  maxsz = *bin_counts = *Loc + BUNDLE - B_ptr; 
  for (i=1 ; i<bins ; i++) {
    bin_counts[i] = Loc[i] + BUNDLE - (B_ptr+=bin_size);
    maxsz = max(bin_counts[i], maxsz);
  }

  if ((maxsz-2) > bin_size) {
    fprintf(stderr,"ERROR: bins too small\n");
    exit(1);
  }

  free(Loc);
}


#if USE_MY_ALLTOALLV
int my_Alltoallv ( void *sendbuf, int *sendcnts,
		   int *sdispls, MPI_Datatype sendtype, 
		   void *recvbuf, int *recvcnts, 
		   int *rdispls, MPI_Datatype recvtype,
		   MPI_Comm comm )
{
  int        i, k;
  MPI_Aint   send_extent, recv_extent;
  
  MPI_Type_extent(sendtype, &send_extent);
  MPI_Type_extent(recvtype, &recv_extent);

  if (!init_Alltoallv_param)
    init_Alltoallv();
       
  for ( k=0; k<PROCS; k++ ) {
    i = k ^ MYPROC;
    MPI_Irecv((void *)((char *)recvbuf+rdispls[i]*recv_extent), 
	      recvcnts[i], recvtype, i, MPIR_ALLTOALLV_TAG,
	      comm, reqarray+i);
    MPI_Send((void *)((char *)sendbuf+sdispls[i]*send_extent), 
	     sendcnts[i], sendtype, i, MPIR_ALLTOALLV_TAG,
	     comm);
  }
  
  MPI_Waitall(PROCS,reqarray,starray);
  
  return (0);
}
#endif

int all_transpose_i(int *A, int *B, int bin_size, int *bin_counts) {
  
  register int i;

  int
    j, t,
    *send_off,
    *send_cnt,
    *recv_off,
    *recv_cnt;
    
  send_cnt = bin_counts;

  send_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(send_off);

  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);


  MPI_Alltoall(send_cnt, 1, MPI_INT,
	       recv_cnt, 1, MPI_INT,
	       MPI_COMM_WORLD);

  send_off[0] = 0;
  recv_off[0] = 0;
  for (i=1 ; i<PROCS ; i++) {
    send_off[i] = send_off[i-1] + bin_size;
    recv_off[i] = recv_off[i-1] + recv_cnt[i-1];
  }

  t = recv_off[PROCS-1] + recv_cnt[PROCS-1];
  
#if USE_MY_ALLTOALLV
  my_Alltoallv(A, send_cnt, send_off, MPI_INT,
	       B, recv_cnt, recv_off, MPI_INT,
	       MPI_COMM_WORLD);
#else
  MPI_Alltoallv(A, send_cnt, send_off, MPI_INT,
		B, recv_cnt, recv_off, MPI_INT,
		MPI_COMM_WORLD);
#endif
		  
  free(recv_off);
  free(recv_cnt);
  free(send_off);
  return(t);
}


int all_transpose_d(double *A, double *B, int bin_size, int *bin_counts) {
  
  register int i;

  int
    j, t,
    *send_off,
    *send_cnt,
    *recv_off,
    *recv_cnt;
    
  send_cnt = bin_counts;

  send_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(send_off);

  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);


  MPI_Alltoall(send_cnt, 1, MPI_INT,
	       recv_cnt, 1, MPI_INT,
	       MPI_COMM_WORLD);

  send_off[0] = 0;
  recv_off[0] = 0;
  for (i=1 ; i<PROCS ; i++) {
    send_off[i] = send_off[i-1] + bin_size;
    recv_off[i] = recv_off[i-1] + recv_cnt[i-1];
  }

  t = recv_off[PROCS-1] + recv_cnt[PROCS-1];
  
#if USE_MY_ALLTOALLV
  my_Alltoallv(A, send_cnt, send_off, MPI_DOUBLE,
	       B, recv_cnt, recv_off, MPI_DOUBLE,
	       MPI_COMM_WORLD);
#else
  MPI_Alltoallv(A, send_cnt, send_off, MPI_DOUBLE,
		B, recv_cnt, recv_off, MPI_DOUBLE,
		MPI_COMM_WORLD);
#endif

  free(recv_off);
  free(recv_cnt);
  free(send_off);
  return(t);
}


int all_select_gatherv_i(int *A, int ni, int total_n, int i) {

  int j, t;
  int *recv_cnt, *recv_off;

  int
    *seqval,
    result;
  
  MPI_Status stat;

  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);

  MPI_Gather(&ni,  1, MPI_INT,
	     recv_cnt, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (MYPROC==0) {
    t = 0;
    for (j=0 ; j<PROCS ; j++) {
      recv_off[j] = t;
      t          += recv_cnt[j];
    }

    seqval = (int *)malloc(t*sizeof(int));
    assert_malloc(seqval);

#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n  xxx\n");
#endif
  }
  

  MPI_Gatherv(A, ni,  MPI_INT,
	      seqval, recv_cnt, recv_off, MPI_INT,
	      0, MPI_COMM_WORLD);

  if (MYPROC==0)
    result = select_mom_i(seqval, t, i);

  free(recv_off);
  free(recv_cnt);
  return(result);
}


int all_select_gatherv_alloc_i(int *A, int ni, int total_n, int i, int *B) {

  int j, t;
  int *recv_cnt, *recv_off;

  int
    *seqval,
    result;
  
  MPI_Status stat;

  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);

  MPI_Gather(&ni,  1, MPI_INT,
	     recv_cnt, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (MYPROC==0) {
    t = 0;
    for (j=0 ; j<PROCS ; j++) {
      recv_off[j] = t;
      t          += recv_cnt[j];
    }

    seqval = (int *)malloc(t*sizeof(int));
    assert_malloc(seqval);

#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n  xxx\n");
#endif
  }
  

  MPI_Gatherv(A, ni,  MPI_INT,
	      seqval, recv_cnt, recv_off, MPI_INT,
	      0, MPI_COMM_WORLD);

  if (MYPROC==0)
    result = select_mom_alloc_i(seqval, t, i, B);

  free(recv_off);
  free(recv_cnt);
  return(result);
}


double all_select_gatherv_d(double *A, int ni, int total_n, int i) {

  int j, t;
  int *recv_cnt, *recv_off;

  double
    *seqval,
    result;
  
  MPI_Status stat;

  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);

  MPI_Gather(&ni,  1, MPI_INT,
	     recv_cnt, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (MYPROC==0) {
    t = 0;
    for (j=0 ; j<PROCS ; j++) {
      recv_off[j] = t;
      t          += recv_cnt[j];
    }

    seqval = (double *)malloc(t*sizeof(double));
    assert_malloc(seqval);

#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n  xxx\n");
#endif
  }
  

  MPI_Gatherv(A, ni,  MPI_DOUBLE,
	      seqval, recv_cnt, recv_off, MPI_INT,
	      0, MPI_COMM_WORLD);

  if (MYPROC==0)
    result = select_mom_d(seqval, t, i);

  free(recv_off);
  free(recv_cnt);
  return(result);
}

double all_select_gatherv_alloc_d(double *A, int ni, int total_n, int i,
				  double *B) {

  int j, t;
  int *recv_cnt, *recv_off;

  double
    *seqval,
    result;
  
  MPI_Status stat;

  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);

  MPI_Gather(&ni,  1, MPI_INT,
	     recv_cnt, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (MYPROC==0) {
    t = 0;
    for (j=0 ; j<PROCS ; j++) {
      recv_off[j] = t;
      t          += recv_cnt[j];
    }

    seqval = (double *)malloc(t*sizeof(double));
    assert_malloc(seqval);

#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n  xxx\n");
#endif
  }
  

  MPI_Gatherv(A, ni,  MPI_DOUBLE,
	      seqval, recv_cnt, recv_off, MPI_INT,
	      0, MPI_COMM_WORLD);

  if (MYPROC==0)
    result = select_mom_alloc_d(seqval, t, i, B);

  free(recv_off);
  free(recv_cnt);
  return(result);
}


void all_gather_picks_i(int *A, int A_size, int **B, int *B_size) {
  /* gather A_size elements from each processor into array B on
     processor 0, and save B_size = Sum(A_size) on processor 0 */

  int
    i,
    *recv_cnt,
    *recv_off;
  
  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);


  MPI_Gather(&A_size,  1, MPI_INT,
	     recv_cnt, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (MYPROC==0) {
    *B_size     = 0;
    for (i=0 ; i<PROCS ; i++) {
      recv_off[i] = *B_size;
      *B_size    += recv_cnt[i];
    }

    *B = (int *)malloc(*B_size*sizeof(int));   /* Free later! */
    assert_malloc(*B);
  }
  
  MPI_Gatherv( A, A_size, MPI_INT,
	      *B, recv_cnt, recv_off, MPI_INT,
	       0, MPI_COMM_WORLD);

  free(recv_off);
  free(recv_cnt);

}

void all_gather_picks_d(double *A, int A_size, double **B, int *B_size) {
  /* gather A_size elements from each processor into array B on
     processor 0, and save B_size = Sum(A_size) on processor 0 */

  int
    i,
    *recv_cnt,
    *recv_off;
  
  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);


  MPI_Gather(&A_size,  1, MPI_INT,
	     recv_cnt, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (MYPROC==0) {
    *B_size     = 0;
    for (i=0 ; i<PROCS ; i++) {
      recv_off[i] = *B_size;
      *B_size    += recv_cnt[i];
    }

    *B = (double *)malloc(*B_size*sizeof(double));   /* Free later! */
    assert_malloc(*B);
  }
  
  MPI_Gatherv( A, A_size, MPI_DOUBLE,
	      *B, recv_cnt, recv_off, MPI_DOUBLE,
	       0, MPI_COMM_WORLD);

  free(recv_off);
  free(recv_cnt);

}

void partition_with_two_i(int *A, int *B, int M, int s0, int s1,
			  int *a0, int *a1) {
  int i, j;
  int *Aptr, *Bptr, *eptr;

  *a0 = 0;
  Aptr = A;
  eptr = Aptr + M;
  Bptr = B;
  while (Aptr < eptr) {
    if (*Aptr<=s0)
      *a0 += 1;
    else
      if (*Aptr<s1)
	*Bptr++ = *Aptr;
    Aptr++;
  }

  *a1 = Bptr-B;
}

void partition_with_two_d(double *A, double *B, int M, double s0, double s1,
			  int *a0, int *a1) {
  int i, j;
  double *Aptr, *Bptr, *eptr;

  *a0 = 0;
  Aptr = A;
  eptr = Aptr + M;
  Bptr = B;
  while (Aptr < eptr) {
    if (*Aptr<=s0)
      *a0 += 1;
    else
      if (*Aptr<s1)
	*Bptr++ = *Aptr;
    Aptr++;
  }

  *a1 = Bptr-B;
}


#if 0
int all_select_min_set_i(int ni, int *A, int total_n, int i,
			   int *sublist, int X) {
#define DATA_TYPE int

    DATA_TYPE
       *Aptr,
       *Eptr,
	result,
        r;

    Aptr = A;
    Eptr = A+ni;

    if (i==1) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = min(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    }
    else {
      /* Need to get the min X elements */
      sort_min_i(A, ni, X); 

      MPI_Gather(A,       X, MPI_INT,
		 sublist, X, MPI_INT,
		 0, MPI_COMM_WORLD);
      
      if (MYPROC==0) {
	result = select_merge_i(sublist, i, X, PROCS);
      }
    }

    return (result);
}
#else
int all_select_min_set_i(int ni, int *A, int total_n, int i,
			   int *sublist, int X) {
#define DATA_TYPE int

    DATA_TYPE
       *Aptr,
       *Eptr,
	result,
        r;

    int shortbuf,
      j,
      newX,
      t,
      *recv_cnt,
      *recv_off;

    Aptr = A;
    Eptr = A+ni;

    if (i==1) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = min(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    }
    else {
      recv_cnt = (int *)malloc(PROCS*sizeof(int));
      assert_malloc(recv_cnt);

      newX = min(ni, X);
      MPI_Allgather(&newX, 1, MPI_INT, recv_cnt, 1, MPI_INT, MPI_COMM_WORLD);

      shortbuf = 0;
      for (j=0 ; j<PROCS ; j++)
	shortbuf += ((recv_cnt[j] < X) ?1:0);

      if (shortbuf==0) {
	/* Need to get the min X elements */
	sort_min_i(A, ni, X); 

	MPI_Gather(A,       X, MPI_INT,
		   sublist, X, MPI_INT,
		   0, MPI_COMM_WORLD);
      
	if (MYPROC==0) {
	  result = select_merge_i(sublist, i, X, PROCS);
	}
      }
      else {
	recv_off = (int *)malloc(PROCS*sizeof(int));
	assert_malloc(recv_off);
	
	/* Need to get the min X elements */
	sort_min_i(A, ni, newX); 

	if (MYPROC==0) {
	  t=0;
	  for (j=0 ; j<PROCS ; j++) {
	    recv_off[j]  = t;
	    t           += recv_cnt[j];
	  }
	}
	
	MPI_Gatherv(A,    newX, MPI_INT,
		    sublist, recv_cnt, recv_off, MPI_INT,
		    0, MPI_COMM_WORLD);
      
	if (MYPROC==0) {
#if 0
	  result = select_merge_i(sublist, i, X, PROCS);
#else
	  result = select_mom_i(sublist, t, i);
#endif
	}
	
	free(recv_off);
      }
      free(recv_cnt);
    }

    return (result);
}
#endif

#if 0
int all_select_max_set_i(int ni, int *A, int total_n, int i,
			 int *sublist, int X) {
#define DATA_TYPE int

    int j;

    DATA_TYPE
        *Aptr,
        *Eptr,
        r,
	result;
    
    Aptr = A;
    Eptr = A+ni;

    if (i==total_n) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = max(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    }
    else {
      /* Need to get the max X elements... */
      sort_max_i(A, ni, X);

      MPI_Gather(A+(ni-X), X, MPI_INT,
		 sublist,  X, MPI_INT,
		 0, MPI_COMM_WORLD);

      if (MYPROC==0) {
	j = i - (total_n - X*PROCS);
	result = select_merge_i(sublist, j, X, PROCS);
      }
    }
    return (result);

#undef DATA_TYPE
}
#else
int all_select_max_set_i(int ni, int *A, int total_n, int i,
			 int *sublist, int X) {
#define DATA_TYPE int

    int j;

    DATA_TYPE
        *Aptr,
        *Eptr,
        r,
	result;
    
    int shortbuf,
      newX,
      t,
      *recv_cnt,
      *recv_off;

    Aptr = A;
    Eptr = A+ni;

    if (i==total_n) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = max(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    }
    else {
      recv_cnt = (int *)malloc(PROCS*sizeof(int));
      assert_malloc(recv_cnt);
      
      newX = min(ni, X);
      MPI_Allgather(&newX, 1, MPI_INT, recv_cnt, 1, MPI_INT, MPI_COMM_WORLD);

      shortbuf = 0;
      for (j=0 ; j<PROCS ; j++)
	shortbuf += ((recv_cnt[j] < X) ?1:0);
      
      if (shortbuf==0) {
	/* Need to get the max X elements... */
	sort_max_i(A, ni, X);
#if DEBUG
	MPI_fprintf(stdout,"(select_max_i) Shortbuf (%d) false. newX: %d\n",
		    shortbuf,newX);
#endif
	MPI_Gather(A+(ni-X), X, MPI_INT,
		   sublist,  X, MPI_INT,
		   0, MPI_COMM_WORLD);

	if (MYPROC==0) {
	  j = i - (total_n - X*PROCS);
	  result = select_merge_i(sublist, j, X, PROCS);
#if DEBUG
	  fprintf(stdout,"PE%3d: (select_max_i) i: %d j: %d\n",MYPROC,i,j);
#endif
	}
	  
      }
      else {
#if DEBUG
	MPI_fprintf(stdout,"(select_max_i) Shortbuf (%d) true. newX: %d\n",
		    shortbuf,newX);
#endif
	recv_off = (int *)malloc(PROCS*sizeof(int));
	assert_malloc(recv_off);
	
	/* Need to get the max X elements */
	sort_max_i(A, ni, newX); 

	if (MYPROC == 0) {
	  t=0;
	  for (j=0 ; j<PROCS ; j++) {
	    recv_off[j]  = t;
	    t           += recv_cnt[j];
	  }
	}
	
	MPI_Gatherv(A+(ni-newX), newX, MPI_INT,
		    sublist, recv_cnt, recv_off, MPI_INT,
		    0, MPI_COMM_WORLD);
      
	if (MYPROC==0) {
#if 0
	  j = i - (total_n - X*PROCS);
	  result = select_merge_i(sublist, j, X, PROCS);
#else
	  j = i - (total_n - t);
	  result = select_mom_i(sublist, t, j);
#endif
#if DEBUG
	  fprintf(stdout,"PE%3d: (select_max_i) i: %d j: %d\n",MYPROC,i,j);
#endif
	}
	
	free(recv_off);
      }
      free(recv_cnt);
    }
    return (result);

#undef DATA_TYPE
}
#endif


#if 0
double all_select_min_set_d(int ni, double *A, int total_n, int i,
			    double *sublist, int X) {
#define DATA_TYPE double

    DATA_TYPE
        *Aptr,
        *Eptr,
	result,
        r;

    Aptr = A;
    Eptr = A+ni;

    if (i==1) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = min(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    }
    else{
      /* Need to get the min X elements */
      sort_min_d(A, ni, X); 

      MPI_Gather(A,        X, MPI_DOUBLE,
		 sublist,  X, MPI_DOUBLE,
		 0, MPI_COMM_WORLD);
	
      if (MYPROC==0) {
	result = select_merge_d(sublist, i, X, PROCS);
      }
    }

    return (result);
}
#else
double all_select_min_set_d(int ni, double *A, int total_n, int i,
			    double *sublist, int X) {
#define DATA_TYPE double

    DATA_TYPE
       *Aptr,
       *Eptr,
	result,
        r;

    int shortbuf,
      j,
      newX,
      t,
      *recv_cnt,
      *recv_off;

    Aptr = A;
    Eptr = A+ni;

    if (i==1) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = min(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    }
    else {
      recv_cnt = (int *)malloc(PROCS*sizeof(int));
      assert_malloc(recv_cnt);
      
      newX = min(ni, X);
      MPI_Allgather(&newX, 1, MPI_INT, recv_cnt, 1, MPI_INT, MPI_COMM_WORLD);

      shortbuf = 0;
      for (j=0 ; j<PROCS ; j++)
	shortbuf += ((recv_cnt[j] < X) ?1:0);

      if (shortbuf==0) {
	/* Need to get the min X elements */
	sort_min_d(A, ni, X); 

	MPI_Gather(A,       X, MPI_DOUBLE,
		   sublist, X, MPI_DOUBLE,
		   0, MPI_COMM_WORLD);
      
	if (MYPROC==0) {
	  result = select_merge_d(sublist, i, X, PROCS);
	}
      }
      else {
	recv_off = (int *)malloc(PROCS*sizeof(int));
	assert_malloc(recv_off);
	
	/* Need to get the min X elements */
	sort_min_d(A, ni, newX); 

	if (MYPROC == 0) {
	  t=0;
	  for (j=0 ; j<PROCS ; j++) {
	    recv_off[j]  = t;
	    t           += recv_cnt[j];
	  }
	}
	
	MPI_Gatherv(A,    newX, MPI_DOUBLE,
		    sublist, recv_cnt, recv_off, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
      
	if (MYPROC==0) {
#if 0
	  result = select_merge_d(sublist, i, X, PROCS);
#else
	  result = select_mom_d(sublist, t, i);
#endif
	}
	
	free(recv_off);
      }
      free(recv_cnt);
    }

    return (result);
}
#endif

#if 0
double all_select_max_set_d(int ni, double *A, int total_n, int i,
			    double *sublist, int X) {
#define DATA_TYPE double

    int j;

    DATA_TYPE
        *Aptr,
        *Eptr,
        r,
	result;
    
    Aptr = A;
    Eptr = A+ni;

    if (i==total_n) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = max(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }
    else {
      /* Need to get the max X elements... */
      sort_max_d(A, ni, X);

      MPI_Gather(A+(ni-X), X, MPI_DOUBLE,
		 sublist,  X, MPI_DOUBLE,
		 0, MPI_COMM_WORLD);

      if (MYPROC==0) {
	j = i - (total_n - X*PROCS);
	result = select_merge_d(sublist, j, X, PROCS);
      }
    }
    return (result);

#undef DATA_TYPE
}
#else
double all_select_max_set_d(int ni, double *A, int total_n, int i,
			    double *sublist, int X) {
#define DATA_TYPE double

    int j;

    DATA_TYPE
        *Aptr,
        *Eptr,
        r,
	result;
    
    int shortbuf,
      newX,
      t,
      *recv_cnt,
      *recv_off;

    Aptr = A;
    Eptr = A+ni;

    if (i==total_n) {
      r = *Aptr++;
      while (Aptr < Eptr)
	r = max(r, *Aptr++);
      MPI_Reduce(&r, &result, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }
    else {
      recv_cnt = (int *)malloc(PROCS*sizeof(int));
      assert_malloc(recv_cnt);
      
      newX = min(ni, X);
      MPI_Allgather(&newX, 1, MPI_INT, recv_cnt, 1, MPI_INT, MPI_COMM_WORLD);

      shortbuf = 0;
      for (j=0 ; j<PROCS ; j++)
	shortbuf += ((recv_cnt[j] < X) ?1:0);
      
      if (shortbuf==0) {
	/* Need to get the max X elements... */
	sort_max_d(A, ni, X);

	MPI_Gather(A+(ni-X), X, MPI_DOUBLE,
		   sublist,  X, MPI_DOUBLE,
		   0, MPI_COMM_WORLD);

	if (MYPROC==0) {
	  j = i - (total_n - X*PROCS);
	  result = select_merge_d(sublist, j, X, PROCS);
	}
	  
      }
      else {
	recv_off = (int *)malloc(PROCS*sizeof(int));
	assert_malloc(recv_off);
	
	/* Need to get the max X elements */
	sort_max_d(A, ni, newX); 

	if (MYPROC == 0) {
	  t=0;
	  for (j=0 ; j<PROCS ; j++) {
	    recv_off[j]  = t;
	    t           += recv_cnt[j];
	  }
	}
		
	MPI_Gatherv(A+(ni-newX), newX, MPI_DOUBLE,
		    sublist, recv_cnt, recv_off, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
      
	if (MYPROC==0) {
#if 0
	  j = i - (total_n - X*PROCS);
	  result = select_merge_d(sublist, j, X, PROCS);
#else
	  j = i - (total_n - t);
	  result = select_mom_d(sublist, t, j);
#endif
	}
	
	free(recv_off);
      }
      free(recv_cnt);
    }
    return (result);

#undef DATA_TYPE
}
#endif


