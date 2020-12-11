/*                                                                    tab:8
 *
 * select_b.c - Parallel Selection Algorithm
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
 * Filename:            select_b.c
 * History:
 */

#include "select_b.h"

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
    
static int    q[MAXPROCS];
static int    d[MAXPROCS];
static int    src[MAXPROCS];
static int    snk[MAXPROCS];
static int    f_rank_src[MAXPROCS];
static int    l_rank_src[MAXPROCS];
static int    f_rank_snk[MAXPROCS];
static int    l_rank_snk[MAXPROCS];

void all_LoadBalance_b_i(int M, int *A,
			 int *N, int total_n) {
/* Load Balance in place */

    register int
	j, t;
    
    int 
      src_cnt, snk_cnt,
      lastrank,
      keys,
      d_start,
      currank,
      curdis;

#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_LoadBalance_b: %d\n",
	    MYPROC, total_n);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    t  = DIVPROCS(total_n);

    if (t*PROCS < total_n)
	t++;

    for (j=0 ; j<PROCS-1 ; j++) 
	q[j] = t;

    q[PROCS-1] = total_n - (t * (PROCS-1));

    for (j=0 ; j<PROCS ; j++) {
	d[j]   = N[j] - q[j];
	src[j] = (d[j]>0);
	snk[j] = (d[j]<0);
    }

    src_cnt=0;	
    snk_cnt=0;
    for (j=0 ; j<PROCS ; j++) {
	if (src[j]) {
	  f_rank_src[j] = src_cnt + 1;
	  src_cnt += d[j];
	  l_rank_src[j] = src_cnt;
	}
	else {
	  f_rank_src[j] = 0;
	  l_rank_src[j] = 0;
	}

	if (snk[j]) {
	  f_rank_snk[j] = snk_cnt + 1;
	  snk_cnt -= d[j];
	  l_rank_snk[j] = snk_cnt;
	}
	else {
	  f_rank_snk[j] = 0;
	  l_rank_snk[j] = 0;
	}
    }

    keys     = 0;  /*keys assigned to be sent or received, as appropriate */
    
    if (src[MYPROC]) { /* calculate sends */
      curdis   = q[MYPROC];
      currank  = f_rank_src[MYPROC];
      lastrank = l_rank_src[MYPROC];
      d_start  = d[MYPROC];
      for (j=0 ; j<PROCS ; j++) {
	if ((!snk[j]) || (keys == d_start) ||
	    (currank < f_rank_snk[j]) || (currank > l_rank_snk[j])) {
	  sendcounts[j] = 0;
	  sdispls[j]    = 0;
	}
	else {
	  /* We need to send some keys to j */
	  t = min(lastrank, l_rank_snk[j]) - currank + 1;
	  sendcounts[j] = t;
	  sdispls[j]    = curdis;
	  curdis       += t;
	  currank      += t;
	  keys         += t;
	}
      }
    }
    else {
      for (j=0 ; j<PROCS ; j++) {
	sendcounts[j] = 0;
	sdispls[j]    = 0;
      }
    }

    if (snk[MYPROC]) { /* calculate receives */
      curdis   = N[MYPROC];
      currank  = f_rank_snk[MYPROC];
      lastrank = l_rank_snk[MYPROC];
      d_start  = -d[MYPROC];
      for (j=0 ; j<PROCS ; j++) {
	if ((!src[j]) || (keys == d_start) ||
	    (currank < f_rank_src[j]) || (currank > l_rank_src[j])) {
	  recvcounts[j] = 0;
	  rdispls[j]    = 0;
	}
	else {
	  /* We need to recv some keys from j */
	  t = min(lastrank, l_rank_src[j]) - currank + 1;
	  recvcounts[j] = t;
	  rdispls[j]    = curdis;
	  curdis       += t;
	  currank      += t;
	  keys         += t;
	}
      }
    }
    else {
      for (j=0 ; j<PROCS ; j++) {
	recvcounts[j] = 0;
	rdispls[j]    = 0;
      }
    }

#if DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    for (j=0 ; j< PROCS ; j++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==j) {
	int i;
	fprintf(outfile,"PE%3d: q:%6d d:%6d src: %1d snk: %1d\n",
		  MYPROC, q[MYPROC], d[MYPROC], src[MYPROC], snk[MYPROC]);
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

    MPI_Alltoallv(A, sendcounts, sdispls, MPI_INT,
		  A, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);

    for (j=0 ; j<PROCS ; j++)
      N[j] = q[j];

}

void all_LoadBalance_b_d(int M, double *A,
			 int *N, int total_n) {
/* Load Balance in place */

    register int
	j, t;
    
    int 
      src_cnt, snk_cnt,
      lastrank,
      keys,
      d_start,
      currank,
      curdis;

#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_LoadBalance_b: %d\n",
	    MYPROC, total_n);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    t  = DIVPROCS(total_n);

    if (t*PROCS < total_n)
	t++;

    for (j=0 ; j<PROCS-1 ; j++) 
	q[j] = t;

    q[PROCS-1] = total_n - (t * (PROCS-1));

    for (j=0 ; j<PROCS ; j++) {
	d[j]   = N[j] - q[j];
	src[j] = (d[j]>0);
	snk[j] = (d[j]<0);
    }

    src_cnt=0;	
    snk_cnt=0;
    for (j=0 ; j<PROCS ; j++) {
	if (src[j]) {
	  f_rank_src[j] = src_cnt + 1;
	  src_cnt += d[j];
	  l_rank_src[j] = src_cnt;
	}
	else {
	  f_rank_src[j] = 0;
	  l_rank_src[j] = 0;
	}

	if (snk[j]) {
	  f_rank_snk[j] = snk_cnt + 1;
	  snk_cnt -= d[j];
	  l_rank_snk[j] = snk_cnt;
	}
	else {
	  f_rank_snk[j] = 0;
	  l_rank_snk[j] = 0;
	}
    }

    keys     = 0;  /*keys assigned to be sent or received, as appropriate */
    
    if (src[MYPROC]) { /* calculate sends */
      curdis   = q[MYPROC];
      currank  = f_rank_src[MYPROC];
      lastrank = l_rank_src[MYPROC];
      d_start  = d[MYPROC];
      for (j=0 ; j<PROCS ; j++) {
	if ((!snk[j]) || (keys == d_start) ||
	    (currank < f_rank_snk[j]) || (currank > l_rank_snk[j])) {
	  sendcounts[j] = 0;
	  sdispls[j]    = 0;
	}
	else {
	  /* We need to send some keys to j */
	  t = min(lastrank, l_rank_snk[j]) - currank + 1;
	  sendcounts[j] = t;
	  sdispls[j]    = curdis;
	  curdis       += t;
	  currank      += t;
	  keys         += t;
	}
      }
    }
    else {
      for (j=0 ; j<PROCS ; j++) {
	sendcounts[j] = 0;
	sdispls[j]    = 0;
      }
    }

    if (snk[MYPROC]) { /* calculate receives */
      curdis   = N[MYPROC];
      currank  = f_rank_snk[MYPROC];
      lastrank = l_rank_snk[MYPROC];
      d_start  = -d[MYPROC];
      for (j=0 ; j<PROCS ; j++) {
	if ((!src[j]) || (keys == d_start) ||
	    (currank < f_rank_src[j]) || (currank > l_rank_src[j])) {
	  recvcounts[j] = 0;
	  rdispls[j]    = 0;
	}
	else {
	  /* We need to recv some keys from j */
	  t = min(lastrank, l_rank_src[j]) - currank + 1;
	  recvcounts[j] = t;
	  rdispls[j]    = curdis;
	  curdis       += t;
	  currank      += t;
	  keys         += t;
	}
      }
    }
    else {
      for (j=0 ; j<PROCS ; j++) {
	recvcounts[j] = 0;
	rdispls[j]    = 0;
      }
    }

    MPI_Alltoallv(A, sendcounts, sdispls, MPI_DOUBLE,
		  A, recvcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);

    for (j=0 ; j<PROCS ; j++)
      N[j] = q[j];

}


static int meds[MAXPROCS];

/*************************************************************/
int all_select_b(int M, int *A,
		 int *N, int total_n, int i)
/*************************************************************/
{

    int myN,
	mom,
	j, k, l,
	t, v,
	result;


#if PROFILE_DETAILED
    double secs;
#endif

#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_b: total_n: %d  i: %d\n",
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
	    fprintf(outfile,"Seq2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq2:  %9.6f\n",secs);
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
	all_LoadBalance_b_i(M, A, N, total_n);
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"LB2   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB2:   %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    

	myN = N[MYPROC];
	fastsort(A, myN); 

	MPI_Gather(&A[(myN + 1) >> 1], 1, MPI_INT,
		   meds,  1, MPI_INT,
		   0, MPI_COMM_WORLD);
	if (MYPROC==0) fastsort(meds, PROCS);

	mom     = meds[PROCS>>1];
	MPI_Bcast(&mom, 1, MPI_INT, 0, MPI_COMM_WORLD);
	k       = binsearch_gt(A, 0, myN - 1, mom);
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
	    bcopy(A + l, A, myN*sizeof(int));
#else
	    for (j=0 ; j< myN ; j++)
		A[j] = A[j+l];
#endif
	}

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"Sel2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
	MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
	result = all_select_b(M, A, N, total_n, (v ? i : (i-t)));
    }

    return (result);
}

int all_select_median_b(int M, int *A)
{
    int *N;
    int i, t;
    int r;

    N = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(N);

    for (i=0 ; i<PROCS ; i++)
      N[i] = M;

    t = PROCS*M;
    
#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_median_b (%d)\n",MYPROC,t);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif


    r = all_select_b(M, A, N, t, (t+1)>>1);
    free(N);
    return (r);
}

int all_select_median_unbalanced_b(int M, int *A,
				   int *N, int total_n)
{
    return all_select_b(M, A, N, total_n, (total_n+1)>>1);
}

