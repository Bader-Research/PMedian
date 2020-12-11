/*                                                                    tab:8
 *
 * select_b_nas.c - Parallel Selection Algorithm for NAS IS
 *
 * 
 * "Copyright (c) 1995,1996 The Regents of the University of Maryland.
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
 * Filename:            select_b_nas.c
 * History:
 */

#include "select_b_nas.h"
#include "select_a.h"

#define THRESH_SEL  (max((PROCS*PROCS),4096))
#define PROFILE_DETAILED 0

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


static int meds[MAXPROCS];

/*************************************************************/
static int all_select_b_nas(int M, int *A,
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
    

    MPI_Barrier(MPI_COMM_WORLD);
#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_b_nas: total_n: %d  i: %d\n",
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
	all_LoadBalance_b(M, A, N, total_n);
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
	fastsort_19(A, myN); 

	MPI_Gather(&A[(myN + 1) >> 1], 1, MPI_INT,
		   meds,  1, MPI_INT,
		   0, MPI_COMM_WORLD);
	if (MYPROC==0) fastsort_19(meds, PROCS);

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
	    bcopy(A+l , A, myN*sizeof(int));
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
	result = all_select_b_nas(M, A, N, total_n, (v ? i : (i-t)));
    }

    return (result);
}

int all_select_median_b_nas(int M, int *A)
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
	fprintf(outfile,"PE %2d: all_select_median_b_nas (%d)\n",MYPROC,t);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    r = all_select_b_nas(M, A, N, t, (t+1)>>1);
    free(N);
    return (r);
}


