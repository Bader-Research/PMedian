/*                                                                    tab:8
 *
 * select_r4.c - Parallel Selection Algorithm
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
 * Filename:            select_r4.c
 * History:
 */

#include "select_r4.h"
#include "select_r.h"
#include "select_c.h"
#include "select_b.h"
#include "select_a.h"
#include "select_seq.h"
#include "select_par.h"
#include "select_randsamp.h"
#include <limits.h>

/* DELTA = 2 sqrt(- log((100-z)/100))  where z is % probability */

#define DELTA      1.0        /* 22.12% prob. */

#define DELTA_MAX  100.0

#if 0
#define DELTA      6.0697085  /* 99.99% prob. */
#define DELTA      5.2565218  /* 99.9 % prob. */
#define DELTA      4.2919321  /* 99   % prob. */
#define DELTA      3.0348543  /* 90   % prob. */
#define DELTA      2.5372725  /* 80   % prob. */
#define DELTA      2.1945139  /* 70   % prob. */
#define DELTA      1.9144615  /* 60   % prob. */
#define DELTA      1.5        /* 43.02% prob. */
#define DELTA      1.0        /* 22.12% prob. */
#define DELTA      0.5        /*  6.05% prob. */
#endif

#define SELECT_MIN_MAX_INIT (PROCS<<1)

#define DELTA_MULT 2.25

#define DEBUG   0

#define PROFILE_TIMING   0

#define USE_SORT         0
#define USE_PICKBUND     0

#if PROFILE_TIMING
#include "timing.h"
#endif

static int THRESH_SEL;
/* #define THRESH_SEL_INIT max((PROCS*PROCS),4096) */
#define THRESH_SEL_INIT (PROCS*PROCS)


/* NOTE:
 *  You must call all_select_r4_init() before a call to all_select_r4...()
 */

int all_gather_picks_r4_i(int *A, int A_size, int *B, int B_max) {
  /* gather A_size elements from each processor into array B on
     processor 0, and return B_size = Sum(A_size) on processor 0 */

  int
    i,
    *recv_cnt,
    *recv_off,
    B_size;
  
  recv_cnt = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_cnt);

  recv_off = (int *)malloc(PROCS*sizeof(int));
  assert_malloc(recv_off);


  MPI_Gather(&A_size,  1, MPI_INT,
	     recv_cnt, 1, MPI_INT,
	     0, MPI_COMM_WORLD);

  if (MYPROC==0) {
    B_size     = 0;
    for (i=0 ; i<PROCS ; i++) {
      recv_off[i] = B_size;
      B_size    += recv_cnt[i];
    }
#if 1
    if (B_size > B_max)
      fprintf(stderr,"ERROR: Bsize > Bmax (%12d %12d)\n",
	      B_size, B_max);
#endif
  }

  MPI_Gatherv(A, A_size, MPI_INT,
	      B, recv_cnt, recv_off, MPI_INT,
	      0, MPI_COMM_WORLD);

  free(recv_off);
  free(recv_cnt);

  return(B_size);
}

int all_select_r4_i(int M, int *A, int total_n, int i) {

#define DATA_TYPE     int
#define DATA_TYPE_MPI MPI_INT

  int
    j,
    n, ni,
    c0, c1,
    t[2],
    a[2];

  DATA_TYPE
    s[2],
    result;

  int
    B_size,
    C_size,
    C_max;
    
  DATA_TYPE
    *B, *C,
    *bufIn, *bufOut, *bufTmp;

  double delta, c_sqrt, Smag, k;
  int thresh;
  int done = 0;
  int SELECT_MIN_MAX = SELECT_MIN_MAX_INIT;


#if PROFILE_TIMING
  all_timer_init();
#endif
    
  B = (DATA_TYPE *)malloc(M * sizeof(DATA_TYPE));
  assert_malloc(B);

  C_max = (int)ceil(PICK_MULT * pow((double)total_n, PICK_EPS));
    
  C = (DATA_TYPE *)malloc(C_max * sizeof(DATA_TYPE));
  assert_malloc(C);

  n  = total_n;
  ni = M;

  bufIn  = A;
  bufOut = B;
    
#if PROFILE_TIMING
  all_timer_reset();
#endif

  thresh = THRESH_SEL;
  while ((n > thresh)&&(!done)) {

    if (i >= n-SELECT_MIN_MAX) {
#if DEBUG
      if (MYPROC==0) {
	fprintf(outfile,"r4: select_max n: %12d i: %12d\n",n, i);
      }
#endif
      result = all_select_max_set_i(ni, bufIn, n, i, bufOut, SELECT_MIN_MAX);
      done = 1;
    }
    else if (i <= SELECT_MIN_MAX) {
#if DEBUG
      if (MYPROC==0) {
	fprintf(outfile,"r4: select_min n: %12d i: %12d\n",n, i);
      }
#endif
      result = all_select_min_set_i(ni, bufIn, n, i, bufOut, SELECT_MIN_MAX);
      done = 1;
    }
    else {
      
#if DEBUG
      if (MYPROC==0) {
	fprintf(outfile,"r4: while n: %12d i: %12d\n",n, i);
      }
#endif
#if USE_PICKBUND
      B_size = all_random_pick_exactbund_i(bufIn, bufOut, ni, n);
#else
      B_size = all_random_pick_exact_i(bufIn, bufOut, ni, n);
#endif
#if 0
      if (B_size==0)
	fprintf(stderr,"PE%3d: B_size = 0\n",MYPROC);
#endif

#if PROFILE_TIMING
      all_timer_mark("random pick");
#endif

      C_size = all_gather_picks_r4_i(bufOut, B_size, C, C_max);
#if 1
      if ((MYPROC==0) && (C_size > C_max))
	fprintf(stderr,"ERROR: (PE%3d): C_size > C_max (%12d %12d)\n",
		MYPROC,C_size,C_max);
#endif
	
#if PROFILE_TIMING
      all_timer_mark("gather picks");
#endif
      
      delta = DELTA;
      if (MYPROC==0) {
	Smag   = (double)C_size * (double)i / (double)n;
	c_sqrt = sqrt((double)C_size);
	k      = delta*c_sqrt;
#if DEBUG_K
	fprintf(outfile,"r4 k: %12d  s: %12d\n",(int)k,C_size);
	fflush(outfile);
#endif
	c0   = (int)floor(Smag-k);
	c1   = (int)floor(Smag+k);
	if (c0 <  0)       c0 = 0;
	if (c1 >= C_size)  c1 = C_size-1;
	
#if USE_SORT	  
	fastsort(C, C_size);
	s[0] = C[c0];
	s[1] = C[c1];
#else
	s[0] = select_mom_alloc_i(C, C_size, c0, bufOut);
	s[1] = select_mom_alloc_i(C, C_size, c1, bufOut);
#endif
      }
#if PROFILE_TIMING
      all_timer_mark("one selects s0 s1");
#endif

      MPI_Bcast(s, 2, DATA_TYPE_MPI, 0, MPI_COMM_WORLD);

#if PROFILE_TIMING
      all_timer_mark("bcast s0 s1");
#endif
      partition_with_two_i(bufIn, bufOut, ni, s[0], s[1], &a[0], &a[1]);

#if PROFILE_TIMING
      all_timer_mark("partition with two");
#endif

      MPI_Allreduce(a, t, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#if PROFILE_TIMING
      all_timer_mark("calculate t0 t1");
#endif

#if DEBUG
      if (MYPROC==0) {
	fprintf(outfile,"AAA i: %8d t0: %8d t1: %8d t2: %8d\n",
		i,t[0],t[1],n - t[0] - t[1]);
      }
#endif

      while ((i <= t[0]) || (i > t[0]+t[1])) {
	delta *= DELTA_MULT;
	if (MYPROC==0) {
	  if (delta > DELTA_MAX) {
	    s[0] = INT_MIN;
	    s[1] = INT_MAX;
	  }
	  else {
#if DEBUG
	    fprintf(outfile,"i: %12d  delta: %9.6f\n",i,delta);
	    fflush(outfile);
#endif
	    k  = delta*c_sqrt;
#if DEBUG_K
	    fprintf(outfile,"r4 k: %12d  s: %12d\n",(int)k,C_size);
	    fflush(outfile);
#endif
	    c0 = (int)floor(Smag-k);
	    c1 = (int)floor(Smag+k);
	    if (c0 <  0)       c0 = 0;
	    if (c1 >= C_size)  c1 = C_size-1;

#if USE_SORT
	    s[0] = C[c0];
	    s[1] = C[c1];
#else
	    s[0] = select_mom_alloc_i(C, C_size, c0, bufOut);
	    s[1] = select_mom_alloc_i(C, C_size, c1, bufOut);
#endif
	  }
	}
	MPI_Bcast(s, 2, DATA_TYPE_MPI, 0, MPI_COMM_WORLD);
	partition_with_two_i(bufIn, bufOut, ni, s[0], s[1], &a[0], &a[1]);
	MPI_Allreduce(a, t, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#if DEBUG
	if (MYPROC==0) {
	  fprintf(outfile,"BBB i: %8d t0: %8d t1: %8d t2: %8d\n",
		  i,t[0],t[1],n - t[0] - t[1]);
	}
#endif
      }

#if DEBUG_K
      if (MYPROC==0) {
	fprintf(outfile,"r4 t1: %12d\n",t[1]);
	fflush(outfile);
      }
#endif
#if PROFILE_TIMING
      all_timer_mark("widen search?");
#endif
#if 0
      if ((MYPROC==0) && (delta != DELTA))
	fprintf(outfile,"delta: %9.3f\n",delta);
#endif

      if ((i>t[0])&&(i<=t[0]+t[1])) { /* MID */
	n        = t[1];
	i       -= t[0];
	ni       = a[1];
      }
      
      bufTmp   = bufIn;
      bufIn    = bufOut;
      bufOut   = bufTmp;

#if 1
      if ((t[0]==0)||(n == t[0]+t[1])) {
	result = all_select_gatherv_alloc_i(bufIn, ni, n, i, bufOut);
	done = 1;
      }
#endif
    }
  }
      
  if (!done) {
    result = all_select_gatherv_alloc_i(bufIn, ni, n, i, bufOut);
    
#if PROFILE_TIMING
    all_timer_mark("select_gather");
#endif
  }

  MPI_Bcast(&result, 1, DATA_TYPE_MPI, 0, MPI_COMM_WORLD);

#if PROFILE_TIMING
  all_timer_mark("bcast result");
#endif

#if PROFILE_TIMING
  all_timer_report(outfile,"select_r4");
#endif

  free(C);
  free(B);
  return (result);

#undef DATA_TYPE_MPI
#undef DATA_TYPE
}


double all_select_r4_d(int M, double *A, int total_n, int i) {}



void all_select_r4_init() {
    THRESH_SEL = THRESH_SEL_INIT;
}

int all_select_median_r4_i(int M, int *A) {
    int t;

    all_select_r4_init();
    t = PROCS*M;
    return all_select_r4_i(M, A, t, (t>>1)+1);
}

double all_select_median_r4_d(int M, double *A)
{
    int t;

    all_select_r4_init();
    t = PROCS*M;
    return all_select_r4_d(M, A, t, (t>>1)+1);
}

