/*                                                                    tab:8
 *
 * select_r.c - Parallel Selection Algorithm
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
 * Filename:            select_r.c
 * History:
 */

#include "select_r.h"
#include "select_c.h"
#include "select_b.h"
#include "select_a.h"
#include "select_seq.h"
#include "select_par.h"

#define OVER          2

#define PROFILE_DETAILED 0
#define PROFILE_GENERAL  0

#define USE_MAX_MIN      1

static int THRESH_SEL;
#define THRESH_SEL_INIT max((PROCS*PROCS),4096)


/* NOTE:
 *  You must call all_select_r_init() before a call to all_select_r...()
 */


#if PROFILE_GENERAL
static double time_select, time_rand, time_trans;
#endif

int all_select_r_i(int M, int *A,
		     int *N, int total_n, int i) {

#define DATA_TYPE int

    int	j, k, 
	t,
        myN;

    DATA_TYPE
	mom,
	result;

    int
      B_size,
      bin_size,
      *bin_counts;
    
    DATA_TYPE
        *B, *C;

#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif


#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_r: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    B_size = OVER*M;
    bin_size = B_size / PROCS;

    B = (DATA_TYPE *)malloc(B_size * sizeof(DATA_TYPE));
    assert_malloc(B);
    C = (DATA_TYPE *)malloc(B_size * sizeof(DATA_TYPE));
    assert_malloc(C);

    bin_counts = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(bin_counts);


    if (total_n < THRESH_SEL) {
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif    

	result = all_select_gather_i(M, A, N, total_n, i);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_select += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
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
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif

	all_randomize_i(A, B, M, bin_size, PROCS, PROCSLOG, bin_counts);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_rand += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"RAN   n: %12d  ",total_n);
	    fprintf(outfile,"Time for RAN:   %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif

	myN = all_transpose_i(B, C, bin_size, bin_counts);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_trans += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"TRANS n: %12d  ",total_n);
	    fprintf(outfile,"Time for TRANS: %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif    

	if (MYPROC==0) mom = select_mom_i(C, myN, myN >> 1);

	MPI_Bcast(&mom, 1, MPI_INT, 0, MPI_COMM_WORLD);
	k       = partition_unk_piv_i(C, myN, mom);
	MPI_Allreduce(&k, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if ((t==total_n)||(t==i)) {
	    result = mom;
#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime() - secs_gen;
	    time_select += secs_gen;
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime() - secs;
	    if (MYPROC==0) {
		fprintf(outfile,"Sel1  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel1:  %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
	}
	else {
	    if (i < t) {
		total_n  = t;
		myN      = k;
		MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
#if USE_MAX_MIN
		result = all_select_max_alloc_i(M, C, N, total_n, i,
						B, t-i);
#else
		result = all_select_c_i(M, C, N, total_n, i);
#endif

	    }
	    else {
		i       -= t;
		total_n -= t;
		myN     -= k;
		MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
#if USE_MAX_MIN
		result = all_select_min_alloc_i(M, C+k, N, total_n, i,
						B, i);
#else
		result = all_select_c_i(M, C+k, N, total_n, i);
#endif
	    }

#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime() - secs_gen;
	    time_select += secs_gen;
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime() - secs;
	    if (MYPROC==0) {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
	}
    }

#if PROFILE_GENERAL
#endif

    free(bin_counts);
    free(C);
    free(B);

    return (result);

#undef DATA_TYPE
}


double all_select_r_d(int M, double *A,
		      int *N, int total_n, int i) {

#define DATA_TYPE double

    int	j, k, 
	t,
        myN;

    DATA_TYPE
	mom,
	result;

    int
      B_size,
      bin_size,
      *bin_counts;
    
    DATA_TYPE
        *B, *C;

#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif


#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_r: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    B_size = OVER*M;
    bin_size = B_size / PROCS;

    B = (DATA_TYPE *)malloc(B_size * sizeof(DATA_TYPE));
    assert_malloc(B);
    C = (DATA_TYPE *)malloc(B_size * sizeof(DATA_TYPE));
    assert_malloc(C);

    bin_counts = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(bin_counts);


    if (total_n < THRESH_SEL) {
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif    

	result = all_select_gather_d(M, A, N, total_n, i);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_select += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"Seq2  n: %12d  ",total_n);
	    fprintf(outfile,"Time for Seq2:  %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
	MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    else {

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif

	all_randomize_d(A, B, M, bin_size, PROCS, PROCSLOG, bin_counts);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_rand += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"RAN   n: %12d  ",total_n);
	    fprintf(outfile,"Time for RAN:   %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif

	myN = all_transpose_d(B, C, bin_size, bin_counts);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_trans += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"TRANS n: %12d  ",total_n);
	    fprintf(outfile,"Time for TRANS: %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif    

	if (MYPROC==0) mom = select_mom_d(C, myN, myN >> 1);

	MPI_Bcast(&mom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	k       = partition_unk_piv_d(C, myN, mom);
	MPI_Allreduce(&k, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if ((t==total_n)||(t==i)) {
	    result = mom;
#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime() - secs_gen;
	    time_select += secs_gen;
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime() - secs;
	    if (MYPROC==0) {
		fprintf(outfile,"Sel1  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel1:  %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
	}
	else {
	    if (i < t) {
		total_n  = t;
		myN      = k;
		MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
#if USE_MAX_MIN
		result = all_select_max_alloc_d(M, C, N, total_n, i,
						B, t-i);
#else
		result = all_select_c_d(M, C, N, total_n, i);
#endif

	    }
	    else {
		i       -= t;
		total_n -= t;
		myN     -= k;
		MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
#if USE_MAX_MIN
		result = all_select_min_alloc_d(M, C+k, N, total_n, i,
						B, i);
#else
		result = all_select_c_d(M, C+k, N, total_n, i);
#endif
	    }

#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime() - secs_gen;
	    time_select += secs_gen;
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime() - secs;
	    if (MYPROC==0) {
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
	}
    }

#if PROFILE_GENERAL
#endif

    free(bin_counts);
    free(C);
    free(B);

    return (result);

#undef DATA_TYPE
}


void all_select_r_init() {
    THRESH_SEL = THRESH_SEL_INIT;
#if USE_MY_ALLTOALLV
    if (!init_Alltoallv_param)
      init_Alltoallv();
#endif
}

int all_select_median_r_i(int M, int *A) {

#define DATA_TYPE int

    int *N;
    DATA_TYPE result;
    int i, t;

    all_select_r_init();

    N = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(N);

    for (i=0 ; i<PROCS ; i++)
      N[i] = M;

    t = PROCS*M;
#if PROFILE_GENERAL
    time_select = 0.0;
    time_rand   = 0.0;
    time_trans  = 0.0;
#endif    
    result = all_select_r_i(M, A, N, t, (t>>1)+1);
#if PROFILE_GENERAL
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) 
	fprintf(outfile,"Select_r i general profile: Select: %9.6f  RAND: %9.6f  TRANS: %9.6f\n",
		time_select, time_rand, time_trans);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    free(N);
    return (result);

#undef DATA_TYPE
}

double all_select_median_r_d(int M, double *A)
{

#define DATA_TYPE double

    int *N;
    DATA_TYPE result;
    int i, t;

    all_select_r_init();

    N = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(N);
    
    for (i=0 ; i<PROCS ; i++)
      N[i] = M;

    t = PROCS*M;
#if PROFILE_GENERAL
    time_select = 0.0;
    time_rand   = 0.0;
    time_trans  = 0.0;
#endif    
    result = all_select_r_d(M, A, N, t, (t>>1)+1);
#if PROFILE_GENERAL
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) 
	fprintf(outfile,"Select_r d general profile: Select: %9.6f  RAND: %9.6f  TRANS: %9.6f\n",
		time_select, time_rand, time_trans);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    free(N);
    return (result);

#undef DATA_TYPE
}

int all_select_median_unbalanced_r_i(int M, int *A,
				       int *N, int total_n) {
    all_select_r_init();
    return all_select_r_i(M, A, N, total_n, (total_n>>1) + 1);

}

double all_select_median_unbalanced_r_d(int M, double *A,
					     int *N, int total_n) {
    all_select_r_init();
    return all_select_r_d(M, A, N, total_n, (total_n>>1) + 1);
}
