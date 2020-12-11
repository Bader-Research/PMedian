/*                                                                    tab:8
 *
 * select_a.c - Parallel Selection Algorithm
 *
 * 
 * "Copyright (c) 1995 The Regents of the University of Maryland.
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
 * Creation Date:       April 17, 1995
 * Filename:            select_a.c
 * History:
 */

#include "select_a.h"
#include "select_seq.h"

#define THRESH_SEL  (max((PROCS*PROCS),4096))
#define PROFILE_DETAILED 0
#define DEBUG            0

static int binsearch_gt(int * list, int min_idx, int max_idx, int target) 
/* return the highest index of the array list with element <= target */
/* return min_idx-1 if all elements of list > target */
/* "THE PRICE IS RIGHT" */    
{
    register int mid_idx;

#if PROFILE_DETAILED
    fprintf(outfile,"PE %2d: in binsearch_gt (%5d, %5d): %5d\n",
	    MYPROC, min_idx, max_idx, target);
#endif
    
    if (max_idx < min_idx)
	return (min_idx-1);
    if (min_idx == max_idx)
	if (list[min_idx] <= target)
	    return (min_idx);
	else
	    return (min_idx-1);
    mid_idx = min_idx + (max_idx - min_idx) / 2;
    if (list[mid_idx] <= target)                 /* This person bid under */
	return binsearch_gt(list, mid_idx+1, max_idx,   target);
    else                                         /* This person bid over  */
	return binsearch_gt(list, min_idx,   mid_idx-1, target);
}


static int    sendcounts[MAXPROCS];
static int    recvcounts[MAXPROCS];
static int    sdispls[MAXPROCS];
static int    rdispls[MAXPROCS];

static int    ps[MAXPROCS];

/*************************************************************/
void all_LoadBalance_a(int M, int *A, int *LB,
		     int *N, int total_n)
/*************************************************************/
{
    register int
	j, t;
    
    int 
      q, qB,
      f_rank, f_proc,
      last_rank,
      myN,
      keys_send,
      keys_recv,
      curdis;

#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_LoadBalance: %d\n",
	    MYPROC, total_n);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    ps[0] = N[0];

    for (j=1 ; j<PROCS ; j++)
	ps[j] = ps[j-1] + N[j];

    qB  = DIVPROCS(total_n);

    if (qB*PROCS < total_n)
	qB++;

    q = qB;

    if (MYPROC==PROCS-1) {
	t  = (q*PROCS) - total_n;
	q -= t;     /* since total_n > P^2, this should be +ve */
    }

    /* q is the number of keys I will have after LB */
    
    myN = N[MYPROC];  /* The number of keys I have initially */
    
    /* Calculate SEND */
    if (MYPROC==0)           /* f_rank is the rank of my first key (1..n) */
      f_rank = 1;
    else
      f_rank = ps[MYPROC-1] + 1;

    f_proc    = (f_rank-1) / qB;/* the processor which gets my first key */
    keys_send = 0;              /* number of keys so far assigned to be sent */
    curdis    = 0;             
    
    for (j=0 ; j<PROCS ; j++) {
      if ((j < f_proc)||(keys_send == myN)) {
	sendcounts[j] = 0;
	sdispls[j]    = 0;
      }
      else {
	/* myN - keys_send:       number of keys I have left
	   qB:                    largest number of keys per node
	   last_rank-f_rank+1 :   number of keys the dest needs    */
	last_rank      = min((j+1)*qB, total_n);
	sendcounts[j]  = t = min(min((myN-keys_send), qB), last_rank-f_rank+1);
	sdispls[j]     = curdis;
	curdis        += t;
	keys_send     += t;
      }
    }

    /* Calculate RECV */
    f_rank = MYPROC*qB + 1;  /* f_rank is the rank of the first key I want */ 
    keys_recv = 0;           /* number of keys so far assigned to be recv'd */
    curdis    = 0;             
    
    for (j=0 ; j<PROCS ; j++) {
      if ((ps[j] < f_rank)||(keys_recv == q)) {
	recvcounts[j] = 0;
	rdispls[j]    = 0;
      }
      else {
	/* q - keys_recv:         I end on this proc
	   N[j]:                  I get all from this proc
	   ps[j]-f_rank+1 :       I start on this proc  */
	recvcounts[j]  = t = min(min((q-keys_recv), N[j]), ps[j]-f_rank+1);
	rdispls[j]     = curdis;
	curdis        += t;
	keys_recv     += t;
      }
    }

#if DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    for (j=0 ; j< PROCS ; j++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==j) {
	int i;
	fprintf(outfile,"PE%3d: myN:%6d qB: %6d q: %6d\n",
		  MYPROC, myN, qB, q);
	for (i=0 ; i<PROCS ; i++) {
	  fprintf(outfile,
		  "PE%3d: i: %2d sc: %4d sd: %4d rc: %4d rd: %4d\n",
		  MYPROC, i,
		  sendcounts[i], sdispls[i],
		  recvcounts[i], rdispls[i]);
	  fflush(outfile);
	}
      }
      MPI_Barrier(MPI_COMM_WORLD);
      sleep(1);
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    MPI_Alltoallv(A,  sendcounts, sdispls, MPI_INT,
		  LB, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&q, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);

#if DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
      int i;
      for (i=0 ; i< min(min(100,myN),q) ; i++)
	fprintf(outfile,"PE%3d: myN:%6d qB: %6d q: %6d A[%3d]: %12d  LB[%3d]: %12d\n",
		MYPROC, myN, qB, q, i, A[i], i, LB[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

}

int all_select_gather_i(int M, int *A, int *N, int total_n, int i) {

  int j, t;

  int
    *seqval,
    result;
  
  MPI_Status stat;

  if (MYPROC==0) {
    seqval = (int*)malloc(total_n*sizeof(int));
    assert_malloc(seqval);

    bcopy(A, seqval, N[0]*sizeof(int));
    t = N[0];
    for (j=1; j<PROCS ; j++) {
      MPI_Recv(seqval+t, N[j], MPI_INT, j, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &stat);
      t += N[j];
    }
#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n\n");
#endif
    result = select_mom_i(seqval, t, i);
    free(seqval);
  }
  else {
    MPI_Send(A, N[MYPROC], MPI_INT, 0, 0, MPI_COMM_WORLD);
  }

  return(result);
}

double all_select_gather_d(int M, double *A, int *N, int total_n, int i) {

  int j, t;

  double
    *seqval,
    result;
  
  MPI_Status stat;

  if (MYPROC==0) {
    seqval = (double *)malloc(total_n*sizeof(double));
    assert_malloc(seqval);

    bcopy(A, seqval, N[0]*sizeof(double));
    t = N[0];
    for (j=1; j<PROCS ; j++) {
      MPI_Recv(seqval+t, N[j], MPI_DOUBLE, j, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &stat);
      t += N[j];
    }
#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n\n");
#endif	
    result = select_mom_d(seqval, t, i);
    free(seqval);
  }
  else {
    MPI_Send(A, N[MYPROC], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  return(result);
}

int all_select_gather_alloc_i(int M, int *A, int *N, int total_n, int i,
			      int *sublist) {

  int j, t;

  int
    *seqval,
    result;
  
  MPI_Status stat;

  if (MYPROC==0) {
    seqval = (int*)malloc(total_n*sizeof(int));
    assert_malloc(seqval);

    bcopy(A, seqval, N[0]*sizeof(int));
    t = N[0];
    for (j=1; j<PROCS ; j++) {
      MPI_Recv(seqval+t, N[j], MPI_INT, j, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &stat);
      t += N[j];
    }
#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n\n");
#endif
    result = select_mom_alloc_i(seqval, t, i, sublist);
    free(seqval);
  }
  else {
    MPI_Send(A, N[MYPROC], MPI_INT, 0, 0, MPI_COMM_WORLD);
  }

  return(result);
}

double all_select_gather_alloc_d(int M, double *A, int *N, int total_n, int i,
				 double *sublist) {

  int j, t;

  double
    *seqval,
    result;
  
  MPI_Status stat;

  if (MYPROC==0) {
    seqval = (double *)malloc(total_n*sizeof(double));
    assert_malloc(seqval);

    bcopy(A, seqval, N[0]*sizeof(double));
    t = N[0];
    for (j=1; j<PROCS ; j++) {
      MPI_Recv(seqval+t, N[j], MPI_DOUBLE, j, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &stat);
      t += N[j];
    }
#if 1
    if (t != total_n)
      fprintf(stderr,"ERROR: t != total_n\n");
#endif	
    result = select_mom_alloc_d(seqval, t, i, sublist);
    free(seqval);
  }
  else {
    MPI_Send(A, N[MYPROC], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  return(result);
}

static int meds[MAXPROCS];

/*************************************************************/
int all_select_a(int M, int *A, int *LB,
		 int *N, int total_n, int i)
/*************************************************************/
{

    int
	mom,
	j, k, l,
	t, v,
	result,
        myN;


#if PROFILE_DETAILED
    double secs;
#endif

#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select:     total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if (total_n < THRESH_SEL) {
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif

	result = all_select_gather_i(M, A, N, total_n, i+1);

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"Seq1  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq1:  %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif	
	MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
	all_LoadBalance_a(M, A, LB, N, total_n);
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"LB1   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB1:   %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    

	myN = N[MYPROC];
	fastsort(LB, myN); 
	MPI_Gather(&LB[(myN + 1) >> 1], 1, MPI_INT,
		   meds,  1, MPI_INT,
		   0, MPI_COMM_WORLD);
	if (MYPROC==0) fastsort(meds, PROCS);

	mom = meds[PROCS>>1];
	MPI_Bcast(&mom, 1, MPI_INT, 0, MPI_COMM_WORLD);

#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);
	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(outfile,"PE%3d: meds[PROCS>>1]: %12d mom: %12d \n",MYPROC,
		meds[PROCS>>1], mom);
	fflush(outfile);
	sleep(1);
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
	k       = binsearch_gt(LB, 0, myN - 1, mom);
	l       = k+1;
	MPI_Allreduce(&l, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	v       = (i <= t);
	j       = (v ? t : (total_n - t) );

	total_n = j;

	if (v)
	    myN = l ;
	else {
	    myN -= l;
#if 1
	    bcopy(LB + l, LB, myN*sizeof(int));
#else
	    for (j=0 ; j< myN ; j++)
		LB[j] = LB[j+l];
#endif	    
	}

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"Sel1  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Sel1:  %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
	result = all_select_a(M, LB, A, N, total_n, (v ? i : (i-t)));

    }

    return (result);
}

int all_select_median_a(int M, int *A)
{
    register int i;
    int *N;
    int *LB; 
    int t;
    int r;

    N = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(N);

    for (i=0 ; i<PROCS ; i++)
      N[i] = M;

    t = PROCS*M;

    LB = (int *)malloc(M*sizeof(int));
    assert_malloc(LB);

    r = all_select_a(M, A, LB, N, t, (t+1)>>1);

    free(LB);
    free(N);
    return (r);
}

int all_select_median_unbalanced_a(int M, int *A,
				   int *N, int total_n)
{
    int *LB;
    int r;

    LB = (int *)malloc(M*sizeof(int));
    assert_malloc(LB);

    r = all_select_a(M, A, LB, N, total_n, (total_n+1)>>1);
    free(LB);
    return (r);
}

int all_select_median_unbalanced_make_data(int M, int *A,
					   int *N, int total_n)
{

    FILE *datafile;
    int i,j,n;
    MPI_Status stat;
    
    int *arr;
    arr = (int *)malloc(M*sizeof(int));
    assert_malloc(arr);

    MPI_Barrier(MPI_COMM_WORLD);

    if (MYPROC==0) {

	datafile = fopen("dfile.dat","w+");

	fprintf(datafile,"%d\n",M);
	fprintf(datafile,"%d\n",PROCS);
	n = N[0];
	fprintf(datafile,"%d\n",n);
	for (j=0 ; j<n ; j++)
	  fprintf(datafile,"%d ",A[j]);
	fprintf(datafile,"\n");

	for (i=1 ; i<PROCS ; i++) {
	    n = N[i];
	    fprintf(datafile,"%d\n",n);
	    /*	    bulk_read(arr, (int *global)(A+i), n*sizeof(int)); */
	    MPI_Recv(arr, n, MPI_INT, i, MPI_ANY_TAG,
		     MPI_COMM_WORLD, &stat);
	    for (j=0 ; j<n ; j++)
		fprintf(datafile,"%d ",arr[j]);
	    fprintf(datafile,"\n");
	}

	fclose(datafile);
    }
    else {
      MPI_Send(A, N[MYPROC], MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    free(arr);
}
