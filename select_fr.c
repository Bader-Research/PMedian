/*                                                                    tab:8
 *
 * select_fr.c - Fast Randomized Selection Algorithm
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
 * Filename:            select_fr.c
 * History:
 */

#include "select_fr.h"
#include "select_r.h"
#include "select_c.h"
#include "select_b.h"
#include "select_a.h"
#include "select_seq.h"
#include "select_par.h"
#include "select_randsamp.h"
#include "sort-psrs.h"
#ifndef SUN
#include <strings.h>
#endif
#include <math.h>

#define DEBUG   0

#define NEW2    0

#define PROFILE_TIMING   0

#if PROFILE_TIMING
#include "timing.h"
#endif


/* NOTE:
 *  You must call all_select_fr_init() before a call to all_select_fr...()
 */


int all_sort_i(int *A, int A_size, int *B, int B_max) {
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

  B_size     = 0;
  if (MYPROC==0) {
    for (i=0 ; i<PROCS ; i++) {
      recv_off[i] = B_size;
      B_size    += recv_cnt[i];
#if DEBUG
      fprintf(outfile,"sort_i: cnt[%3d]: %8d off[%3d]: %8d B_size: %8d\n",
	      i, recv_cnt[i], i, recv_off[i], B_size);
#endif
    }

#if 1
    if (B_size > B_max)
      fprintf(stderr,"ERROR: (all_sort_i): B_size > B_max (%12d %12d)\n",
	      B_size,B_max);
#endif
  }
    
  MPI_Gatherv(A, A_size, MPI_INT,
	      B, recv_cnt, recv_off, MPI_INT,
	      0, MPI_COMM_WORLD);

  if (MYPROC==0) 
    fastsort(B, B_size);

  free(recv_off);
  free(recv_cnt);

  return(B_size);
}


void partition_with_two_fr_i(int *A, int *B, int M, int s0, int s1,
			  int *a0, int *a1) {
  int *Aptr, *eptr;
  int *a0ptr, *a1ptr, *a2ptr;
  int *D;

  if (M>0) {
    D = (int *)malloc(M*sizeof(int));
    assert_malloc(D);

    a0ptr = B;
    a1ptr = D;
    a2ptr = B+ M-1;
  
    *a0 = 0;
    Aptr = A;
    eptr = Aptr + M;
    while (Aptr < eptr) {
      if (*Aptr<s0) {
	*a0ptr++ = *Aptr;
      }
      else {
	if (*Aptr<=s1) {
	  *a1ptr++ = *Aptr;
	}
	else {
	  *a2ptr-- = *Aptr;
	}
      }
      Aptr++;
    }

    *a0 = a0ptr-B;
    *a1 = a1ptr-D;

#if 1
    if ((B + *a0 + *a1 - 1) != a2ptr) {
      fprintf(outfile,"a0: %8d a1: %8d a2: %8d sum: %8d  M: %8d\n",
	      *a0, *a1, A+M-a2ptr, *a0+*a1+(A+M-a2ptr), M);
    }
#endif
  
#if 0
    bcopy(D, a0ptr, *a1*sizeof(int));
#else
    bcopy(B,       A            , *a0*sizeof(int));
    bcopy(D,       A + *a0      , *a1*sizeof(int));
    bcopy(a2ptr+1, A + *a0 + *a1, (M - (*a0+*a1))*sizeof(int));
#endif
    

    free(D);
  }
  else {
    *a0 = 0;
    *a1 = 0;
  }
  return;
}


int all_select_fr_i(int M, int *A, int total_n, int i) {

#define DATA_TYPE     int
#define DATA_TYPE_MPI MPI_INT

    int
      l, r, n, ni,
      a[2], t[2],
#if NEW2
      c_tot,
      c_ps_first, c_ps_last,
#endif
      c_less, c_mid;
      

    DATA_TYPE
#if NEW2
      cbuf[2],
#endif
      k[2],
      c0, c1,
      result;

    int
      B_size,
      C_size,
      C_max;
    
    DATA_TYPE
        *B, *C;

    double Smag;
    
    B = (DATA_TYPE *)malloc(M * sizeof(DATA_TYPE));
    assert_malloc(B);

    C_max = (int)ceil(PICK_MULT * pow((double)total_n, PICK_EPS));
    
    C = (DATA_TYPE *)malloc(C_max * sizeof(DATA_TYPE));
    assert_malloc(C);

    n = total_n;
    l = 0;
    r = M-1;
    
    while (n > PROCS*PROCS) {

      
      /* Step 0 */      
      ni = r - l + 1;

#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 0:  M: %8d n: %8d  l: %8d  r: %8d ni: %8d i:%8d\n",
		PROCS,M,n,l,r,ni, i);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

#if DEBUG
#if 1
      MPI_Barrier(MPI_COMM_WORLD);
      {
	int mn, mx, z;
	mx = mn = A[l];
	for (z=1 ; z<ni ; z++) {
	  mn = min(mn, A[l+z]);
	  mx = max(mx, A[l+z]);
	}
	MPI_fprintf(outfile,"Min: %12d  Max: %12d\n",mn,mx);
      }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
      
      /* Step 1: Collect a sample S_i from L_i[l,r] by picking
	 n_i n^\eps / n elements at random on P_i between l and r */

      B_size = all_random_pick_exact_i(A+l, B, ni, n);

#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      /*     if (MYPROC==0) */
      MPI_fprintf(outfile,"(%3d)Step 1:  B_size: %8d \n",
		  PROCS, B_size);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

#if 0
      /* Step 1.5: What if no picks were made? */
      if (B_size == 0) {
	bcopy(A+l, B, ni*sizeof(int));
	B_size = ni;
      }
#endif

#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 1.5:  B_size: %8d \n",
		PROCS,B_size);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      /* Step 2: S = ParallelSort(S_i, p); */

#if NEW2
      C_size = all_sort_psrs_i(B, B_size, C); /* C is distributed */
      /* all_sort_psrs_check_i(C, C_size); */
#else
      C_size = all_sort_i(B, B_size, C, C_max);
#endif
#if 1
      if (C_size > C_max)
	fprintf(stderr,"ERROR: (PE%3d): C_size > C_max (%12d %12d)\n",
		MYPROC,C_size,C_max);
#endif
      
#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 2:  C_size: %8d \n",
		PROCS,C_size);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

#if NEW2
      /* Step 3: Pick k1, k2, from S with ranks... */
      MPI_Allreduce(&C_size, &c_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      Smag = (double)c_tot * (double)i / (double)n;
      c0   = (int)ceil(Smag - sqrt((double)c_tot * log((double)n)));
      c1   = (int)ceil(Smag + sqrt((double)c_tot * log((double)n)));
      if (c0 < 0)       c0 = 0;
      if (c1 >= c_tot)  c1 = c_tot-1;

      MPI_Scan(&C_size, &c_ps_last, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      c_ps_last--;
      c_ps_first = c_ps_last - C_size + 1;
      
      if ((c0 >= c_ps_first) && (c0 <= c_ps_last))
	cbuf[0] = C[c0 - c_ps_first];
      else
	cbuf[0] = 0;
      if ((c1 >= c_ps_first) && (c1 <= c_ps_last))
	cbuf[1] = C[c1 - c_ps_first];
      else
	cbuf[1] = 0;
#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_fprintf(outfile,"(%3d)Step 3: c0: %8d c1:%8d\n",PROCS, c0, c1);
      MPI_fprintf(outfile,"(%3d)Step 3: ps(%8d %8d) cbuf0 cbuf1:(%12d %12d)\n",
		  PROCS,c_ps_first, c_ps_last, cbuf[0],cbuf[1]);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

#else
      /* On P0*/
      if (MYPROC==0) {
	/* Step 3: Pick k1, k2, from S with ranks... */


	Smag = (double)C_size * (double)i / (double)n;
	c0   = (int)ceil(Smag - sqrt((double)C_size * log((double)n)));
	c1   = (int)ceil(Smag + sqrt((double)C_size * log((double)n)));
#if DEBUG_K
	fprintf(outfile,"fr k: %12d  s: %12d\n",
		(int)ceil(sqrt((double)C_size * log((double)n))),C_size);
	fflush(outfile);
#endif
#if DEBUG
	fprintf(outfile,"(%3d)Step 3:  C_size: %8d  (c0 c1): (%8d  %8d)\n",
		PROCS,C_size, c0,c1);
#endif
	if (c0 < 0)       c0 = 0;
	if (c1 >= C_size) c1 = C_size-1;
	k[0] = C[c0];
	k[1] = C[c1];
      }

#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 3:  k[0]: %12d  k[1]: %12d\n",
		PROCS,k[0],k[1]);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

      /* Step 4: Broadcast k1, k2.
	 The rank to be found will be in [k1,k2] w.h.p. */

#if NEW2
      MPI_Allreduce(cbuf, k, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
      MPI_Bcast(k, 2, DATA_TYPE_MPI, 0, MPI_COMM_WORLD);
#endif
	
#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      /* if (MYPROC==0) */
      MPI_fprintf(outfile,"(%3d)Step 4:  k[0]: %12d  k[1]: %12d\n",
		  PROCS,k[0],k[1]);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      /* Step 5: Partition L_i between l and r into < k1, [k1,k2] and >k2
	 to give counts less, middle and high and splitters s1, s2 */
      /* B is just temp space */
      partition_with_two_fr_i(A+l, B, ni, k[0], k[1], &a[0], &a[1]);

#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      /*      if (MYPROC==0) */
      MPI_fprintf(outfile,"(%3d)Step 5:  a[0]: %8d  a[1]: %8d\n",
	      PROCS,a[0],a[1]);
#if 1
      {
	int z;
	for (z=0 ; z<a[0] ; z++)
	  if (A[l+z] >= k[0])
	    fprintf(stderr,"ERROR  A: (PE%3d) A[%8d]: %12d >= k[0]: %12d\n",
		    MYPROC, z, A[l+z], k[0]);
	for (z=a[0] ; z<(a[0]+a[1]) ; z++) {
	  if (A[l+z] < k[0])
	    fprintf(stderr,"ERROR B1: (PE%3d) A[%8d]: %12d <  k[0]: %12d\n",
		    MYPROC, z, A[l+z], k[0]);
	  if (A[l+z] > k[1])
	    fprintf(stderr,"ERROR B2: (PE%3d) A[%8d]: %12d >  k[1]: %12d\n",
		    MYPROC, z, A[l+z], k[1]);
	}
	for (z=(a[0]+a[1]) ; z<ni ; z++)
	  if (A[l+z] <= k[1])
	    fprintf(stderr,"ERROR  C: (PE%3d) A[%8d]: %12d <= k[1]: %12d\n",
		    MYPROC, z, A[l+z], k[1]);
      }
#endif
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      /* Step 6: c_mid  = Combine(middle, add); */
      /* Step 7: c_less = Combine(less, add); */

      MPI_Allreduce(a, t, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      c_mid  = t[1];
      c_less = t[0];

#if DEBUG_K
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0) {
	fprintf(outfile,"fr t0: %12d  t1: %12d  t2: %12d\n",
		c_less, c_mid, n - (c_less+c_mid));
	fflush(outfile);
      }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 6:  c_mid:  %12d\n",
		PROCS,c_mid);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 7:  c_less: %12d\n",
		PROCS,c_less);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      /* Step 8:  If (rank \in (c_less, c_mid]) */
      
      if ((i > c_less) && (i<=c_less+c_mid)) {
	n  = c_mid;
	l += a[0];
	r  = l + a[1] - 1;
	i -= c_less;
#if DEBUG_K
	if (MYPROC==0)
	  fprintf(outfile,"fr \t\t\t MID\n");
#endif
#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYPROC==0)
	  fprintf(outfile,"(%3d)Step 8:   MID\n", PROCS);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
      }
      else {
	if (i <= c_less) {
	  r  = l + a[0] - 1;
	  n  = c_less;
#if DEBUG_K
	  if (MYPROC==0)
	    fprintf(outfile,"fr \t\t\tLESS\n");
#endif
#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYPROC==0)
	  fprintf(outfile,"(%3d)Step 8:  LESS\n", PROCS);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
	else {
	  n -= (c_less+c_mid);
	  l += a[0] + a[1];
	  i -= (c_less+c_mid);
#if DEBUG_K
	  if (MYPROC==0)
	    fprintf(outfile,"fr \t\t\tHIGH\n");
#endif
#if DEBUG
	  MPI_Barrier(MPI_COMM_WORLD);
	  if (MYPROC==0)
	    fprintf(outfile,"(%3d)Step 8:  HIGH\n", PROCS);
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
      }
#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 8\n",PROCS);
      MPI_Barrier(MPI_COMM_WORLD);
#endif

    }

    /* Step 9: L = Gather(L_i[l,r]) */
    ni = r - l + 1;

    result = all_select_gatherv_i(A+l, ni, n, i);

#if DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (MYPROC==0)
	fprintf(outfile,"(%3d)Step 9:  result: %12d\n",
		PROCS,result);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    /* Step 10: result = Broadcast(q) */
    
    MPI_Bcast(&result, 1, DATA_TYPE_MPI, 0, MPI_COMM_WORLD);

#if DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0)
      fprintf(outfile,"(%3d)Step 10: result: %12d\n",
	      PROCS,result);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    free(C);
    free(B);
    return (result);

#undef DATA_TYPE_MPI
#undef DATA_TYPE
}

double all_select_fr_d(int M, double *A, int total_n, int i) {}


void all_select_fr_init() {
}

int all_select_median_fr_i(int M, int *A) {
    int t;

    all_select_fr_init();
    t = PROCS*M;
    return all_select_fr_i(M, A, t, (t>>1)+1);
}

double all_select_median_fr_d(int M, double *A)
{
    int t;

    all_select_fr_init();
    t = PROCS*M;
    return all_select_fr_d(M, A, t, (t>>1)+1);
}



