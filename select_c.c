/*                                                                    tab:8
 *
 * select_c.c - Parallel Selection Algorithm
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
 * Filename:            select_c.c
 * History:
 */

#include "select_c.h"
#include "select_b.h"
#include "select_a.h"
#include "select_seq.h"

static int THRESH_SEL;
#define THRESH_SEL_INIT max((PROCS*PROCS),4096)
static int SELECT_MIN_MAX;
#define SELECT_MIN_MAX_INIT max(PROCS,PROCS*PROCS)

/* NOTE:
 *  You must call all_select_c_init() before a call to all_select_c...()
 */

#define PROFILE_DETAILED 0
#define PROFILE_GENERAL  0

#define DEBUG            0

#if PROFILE_GENERAL
static double time_select, time_lb;
#endif

static int    meds_i[MAXPROCS];
static double meds_d[MAXPROCS];

int all_select_c_i(int M, int *A,
		     int *N, int total_n, int i) {

#define DATA_TYPE int

    int	k, 
	t,
        myN;

    DATA_TYPE
        g,
	mom,
	result;

#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif


#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
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

	all_LoadBalance_b_i(M, A, N, total_n);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_lb += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
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
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif    

	myN = N[MYPROC];

	g = select_mom_i(A, myN, myN >> 1);

	MPI_Gather(&g,      1, MPI_INT,
		   meds_i,  1, MPI_INT,
		   0, MPI_COMM_WORLD);

	if (MYPROC==0) mom = select_mom_i(meds_i, PROCS, PROCS>>1);

	MPI_Bcast(&mom, 1, MPI_INT, 0, MPI_COMM_WORLD);
	k       = partition_unk_piv_i(A, myN, mom);
	MPI_Allreduce(&k, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if (t==total_n) {
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
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
	}
	else {
	    if (i <= t) {
		total_n  = t;
		myN      = k;
	    }
	    else {
		total_n -= t;
		i       -= t;
		myN     -= k;
		
		bcopy(A+k, A, myN*sizeof(DATA_TYPE));
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
	    MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
	    result = all_select_c_i(M, A, N, total_n, i);
	}
    }

#if PROFILE_GENERAL
#endif

    return (result);

#undef DATA_TYPE
}

double all_select_c_d(int M, double *A,
		      int *N, int total_n, int i) {

#define DATA_TYPE double

    int	k, 
	t,
      myN;

    DATA_TYPE
        g,
	mom,
	result;
    
#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif


#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (total_n <= THRESH_SEL) {
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

	all_LoadBalance_b_d(M, A, N, total_n);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_lb += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
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
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif    

	myN = N[MYPROC];

	g = select_mom_d(A, myN, myN >> 1);

	MPI_Gather(&g,     1, MPI_DOUBLE,
		   meds_d, 1, MPI_DOUBLE,
		   0, MPI_COMM_WORLD);

	if (MYPROC==0) mom = select_mom_d(meds_d, PROCS, PROCS>>1);

	MPI_Bcast(&mom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	k       = partition_unk_piv_d(A, myN, mom);
	MPI_Allreduce(&k, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if (t==total_n) {
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
		fprintf(outfile,"Sel2  n: %12d  ",total_n);
		fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif    
	}
	else {
	    if (i <= t) {
		total_n  = t;
		myN      = k;
	    }
	    else {
		total_n -= t;
		i       -= t;
	        myN     -= k;

		bcopy(A+k, A, myN*sizeof(DATA_TYPE));
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
	    MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
	    result = all_select_c_d(M, A, N, total_n, i);
	}
    }

#if PROFILE_GENERAL
#endif

    return (result);

#undef DATA_TYPE
}

int all_select_min_alloc_i(int M, int *A,
			   int *N, int total_n, int i,
			   int *sublist, int X) {
#define DATA_TYPE int

    int myN,
	j, t;

    DATA_TYPE
	piv,
	result,
        r,
       *seqval;

    myN = N[MYPROC];
    t = myN;

    if (i==1) {
	r = A[0];
	for (j=1 ; j<t ; j++)
	    r = min(r, A[j]);
	MPI_Allreduce(&r, &result, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }
    else{
#if 0
	piv = select_mom_alloc_i(A, t, i, sublist);
	j   = partition_i(A, t, piv);
#else
	/* Need to get the min X elements */
	sort_min_i(A, myN, X); 
#endif

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYPROC==0) {
	  MPI_Status stat;
	  seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	  assert_malloc(seqval);

	  bcopy(A, seqval, X*sizeof(DATA_TYPE));
	  t = X;
	  for (j=1; j<PROCS ; j++) {
	    MPI_Recv(seqval+t, X, MPI_INT, j, MPI_ANY_TAG,
		     MPI_COMM_WORLD, &stat);
	    t += X;
	  }
#if 0
	  result = select_mom_alloc_i(seqval, t, i,sublist);
#else
	  result = select_merge_i(seqval, i, X, PROCS);
#endif
	  free(seqval);
	}
	else {
	  MPI_Send(A, X, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

        MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    return (result);
}

int all_select_max_alloc_i(int M, int *A,
			   int *N, int total_n, int i,
			   int *sublist, int X) {
#define DATA_TYPE int

    int myN,
	j, t,
	offset,
	last;

    DATA_TYPE
	piv,
        r,
	result,
	*seqval;
    

    myN = N[MYPROC];
    t = myN;

    if (i==total_n) {
	r = A[0];
	for (j=1 ; j<t ; j++)
	    r = max(r, A[j]);
	MPI_Allreduce(&r, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    else {
#if 0
	piv = select_mom_alloc_i(A, t, i, sublist);
	j   = partition_i(A, t, piv);
#else
	/* Need to get the max X elements... */
	sort_max_i(A, myN, X);
#endif

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYPROC==0) {
	  MPI_Status stat;
	  seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	  assert_malloc(seqval);

	  bcopy(A + (N[0] - X), seqval, X*sizeof(DATA_TYPE));
	  t = X;
	  for (j=1; j<PROCS ; j++) {
	    MPI_Recv(seqval+t, X, MPI_INT, j, MPI_ANY_TAG,
		     MPI_COMM_WORLD, &stat);
	    t += X;
	  }

	  j = i - (total_n - X*PROCS);
#if 0
	  result = select_mom_alloc_i(seqval, t, j,sublist);
#else
	  result = select_merge_i(seqval, j, X, PROCS);
#endif
	  free(seqval);
	}
	else {
	  MPI_Send(A + (N[MYPROC] - X), X, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

        MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    return (result);

#undef DATA_TYPE
}

double all_select_min_alloc_d(int M, double *A,
			      int *N, int total_n, int i,
			      double *sublist, int X) {
#define DATA_TYPE double

    int myN,
	j, t;

    DATA_TYPE
	piv,
        r,
	result,
	*seqval;
    

    myN = N[MYPROC];
    t = myN;
    
    if (i==1) {
	r = A[0];
	for (j=1 ; j<t ; j++)
	    r = min(r, A[j]);
	MPI_Allreduce(&r, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    else {
#if 0
      piv = select_mom_alloc_d(A, t, i, sublist);
      j   = partition_d(A, t, piv);
#else
      sort_min_d(A, myN, X);
#endif
	MPI_Barrier(MPI_COMM_WORLD);

	if (MYPROC==0) {
	  MPI_Status stat;
	  seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	  assert_malloc(seqval);

	  bcopy(A, seqval, X*sizeof(double));
	  t = X;
	  for (j=1; j<PROCS ; j++) {
	    MPI_Recv(seqval+t, X, MPI_DOUBLE, j, MPI_ANY_TAG,
		     MPI_COMM_WORLD, &stat);
	    t += X;
	  }
#if 0
	  result = select_mom_alloc_d(seqval, t, i,sublist);
#else
	  result = select_merge_d(seqval, i, X, PROCS);
#endif
	  free(seqval);
	}
	else {
	  MPI_Send(A, X, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

        MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    return (result);

#undef DATA_TYPE
}

double all_select_max_alloc_d(int M, double *A,
			   int *N, int total_n, int i,
			   double *sublist, int X) {
#define DATA_TYPE double

    int myN,
	j, t,
	offset,
	last;

    DATA_TYPE
	piv,
        r,
	result,
	*seqval;
    

    myN = N[MYPROC];
    t = myN;

    if (i==total_n) {
	r = A[0];
	for (j=1 ; j<t ; j++)
	    r = max(r, A[j]);
	MPI_Allreduce(&r, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    else {
#if 0
      piv = select_mom_alloc_d(A, t, i, sublist);
      j   = partition_d(A, t, piv);
#else
      sort_max_d(A, myN, X);
#endif
	MPI_Barrier(MPI_COMM_WORLD);

	if (MYPROC==0) {
	  MPI_Status stat;
	  seqval = (DATA_TYPE*)malloc(total_n*sizeof(DATA_TYPE));
	  assert_malloc(seqval);

	  bcopy(A + (N[0] - X), seqval, X*sizeof(double));
	  t = X;
	  for (j=1; j<PROCS ; j++) {
	    MPI_Recv(seqval+t, X, MPI_DOUBLE, j, MPI_ANY_TAG,
		     MPI_COMM_WORLD, &stat);
	    t += X;
	  }

	  j = i - (total_n - X*PROCS);
#if 0
	  result = select_mom_alloc_d(seqval, t, j,sublist);
#else
	  result = select_merge_d(seqval, j, X, PROCS);
#endif
	  free(seqval);
	}
	else {
	  MPI_Send(A + (N[MYPROC] - X), X, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

        MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    return (result);

#undef DATA_TYPE
}

int all_select_c_alloc_i(int M, int *A,
			 int *N, int total_n, int i,
			 int *sublist) {

#define DATA_TYPE int

    int	myN,
        k, 
	t;

    DATA_TYPE
      g,
      mom,
      result;

#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
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
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif

	result = all_select_gather_alloc_i(M, A, N, total_n, i, sublist);

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

	all_LoadBalance_b_i(M, A, N, total_n);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_lb += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"LB2   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB2:   %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

	if (i <= SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime();
#endif    

	    result = all_select_min_alloc_i(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

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
		fprintf(outfile,"MIN   n: %12d  ",total_n);
		fprintf(outfile,"Time for MIN:   %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
	else if (i >= total_n - SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime();
#endif    

	    result = all_select_max_alloc_i(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

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
		fprintf(outfile,"MAX   n: %12d  ",total_n);
		fprintf(outfile,"Time for MAX:   %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif
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

	    myN = N[MYPROC];

	    g = select_mom_alloc_i(A, myN, myN >> 1, sublist);
	    
	    MPI_Gather(&g,     1, MPI_INT,
		       meds_i, 1, MPI_INT,
		       0, MPI_COMM_WORLD);

	    if (MYPROC==0) mom = select_mom_alloc_i(meds_i, PROCS, PROCS>>1,sublist);

	    MPI_Bcast(&mom, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    k       = partition_unk_piv_i(A, myN, mom);
	    MPI_Allreduce(&k, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	    if (t==total_n) {
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
		    fprintf(outfile,"Sel2  n: %12d  ",total_n);
		    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif    
	    }
	    else {
		if (i <= t) {
		    total_n  = t;
		    myN      = k;
		}
		else {
		    total_n -= t;
		    i       -= t;
		    myN     -= k;

		    bcopy(A+k, A, myN*sizeof(DATA_TYPE));
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

		MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
		result = all_select_c_alloc_i(M, A, N, total_n, i,sublist);
	    }
	}
    }

    return (result);

#undef DATA_TYPE
}

double all_select_c_alloc_d(int M, double *A,
			    int *N, int total_n, int i,
			    double* sublist) {

#define DATA_TYPE double

    int	k, 
	t,
        myN;

    DATA_TYPE
        g,
	mom,
	result;
    
#if PROFILE_DETAILED
    double secs;
#endif
#if PROFILE_GENERAL
    double secs_gen;
#endif

#if PROFILE_DETAILED
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) {
	fprintf(outfile,"PE %2d: all_select_c: total_n: %d  i: %d\n",
	    MYPROC, total_n, i);
	fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (total_n <= THRESH_SEL) {
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime();
#endif

	result = all_select_gather_alloc_d(M, A, N, total_n, i, sublist); 

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

	all_LoadBalance_b_d(M, A, N, total_n);

#if PROFILE_GENERAL
	MPI_Barrier(MPI_COMM_WORLD);
	secs_gen = MPI_Wtime() - secs_gen;
	time_lb += secs_gen;
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
#if PROFILE_DETAILED
	MPI_Barrier(MPI_COMM_WORLD);
	secs = MPI_Wtime() - secs;
	if (MYPROC==0) {
	    fprintf(outfile,"LB2   n: %12d  ",total_n);
	    fprintf(outfile,"Time for LB2:   %9.6f\n",secs);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif    

	if (i <= SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime();
#endif    

	    result = all_select_min_alloc_d(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

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
		fprintf(outfile,"MIN   n: %12d  ",total_n);
		fprintf(outfile,"Time for MIN:   %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
	else if (i >= total_n - SELECT_MIN_MAX) {
#if PROFILE_DETAILED
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs = MPI_Wtime();
#endif    
#if PROFILE_GENERAL
	    MPI_Barrier(MPI_COMM_WORLD);
	    secs_gen = MPI_Wtime();
#endif    

	    result = all_select_max_alloc_d(M, A, N, total_n, i,
					    sublist, SELECT_MIN_MAX);

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
		fprintf(outfile,"MAX   n: %12d  ",total_n);
		fprintf(outfile,"Time for MAX:   %9.6f\n",secs);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
#endif
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

	    myN = N[MYPROC];

	    g = select_mom_alloc_d(A, myN, myN >> 1, sublist);

	    MPI_Gather(&g,     1, MPI_DOUBLE,
		       meds_d, 1, MPI_DOUBLE,
		       0, MPI_COMM_WORLD);

	    if (MYPROC==0) mom = select_mom_alloc_d(meds_d, PROCS, PROCS>>1, sublist);
	    
	    MPI_Bcast(&mom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	    k       = partition_unk_piv_d(A, myN, mom);
	    MPI_Allreduce(&k, &t, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	    if (t==total_n) {
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
		    fprintf(outfile,"Sel2  n: %12d  ",total_n);
		    fprintf(outfile,"Time for Sel2:  %9.6f\n",secs);
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif    
	    }
	    else {
		if (i <= t) {
		    total_n  = t;
		    myN      = k;
		}
		else {
		    total_n -= t;
		    i       -= t;
		    myN     -= k;

		    bcopy(A+k, A, myN*sizeof(DATA_TYPE));
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
		MPI_Allgather(&myN, 1, MPI_INT, N, 1, MPI_INT, MPI_COMM_WORLD);
		result = all_select_c_alloc_d(M, A, N, total_n, i,sublist);
	    }
	}
    }

    return (result);

#undef DATA_TYPE
}

void all_select_c_init() {
    THRESH_SEL = THRESH_SEL_INIT;
    SELECT_MIN_MAX = SELECT_MIN_MAX_INIT;
}

int all_select_median_c_i(int M, int *A) {

#define DATA_TYPE int

    int *N;
    DATA_TYPE *sublist;
    DATA_TYPE result;
    int i, t;

    all_select_c_init();

    N = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(N);

    sublist = (DATA_TYPE *)malloc(M*sizeof(DATA_TYPE));
    assert_malloc(sublist);

    for (i=0 ; i<PROCS ; i++)
      N[i] = M;

    t = PROCS*M;
#if PROFILE_GENERAL
    time_select = 0.0;
    time_lb = 0.0;
#endif    
    result = all_select_c_alloc_i(M, A, N, t, (t>>1)+1,sublist);
#if PROFILE_GENERAL
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) 
	fprintf(outfile,"Select_c general profile: Select: %9.6f  LB: %9.6f\n",
		time_select, time_lb);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    free(sublist);
    free(N);
    return (result);

#undef DATA_TYPE
}

double all_select_median_c_d(int M, double *A)
{

#define DATA_TYPE double

    int *N;
    DATA_TYPE *sublist;
    DATA_TYPE result;
    int i, t;

    all_select_c_init();

    N = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(N);
    
    sublist = (DATA_TYPE *)malloc(M*sizeof(DATA_TYPE));
    assert_malloc(sublist);

    for (i=0 ; i<PROCS ; i++)
      N[i] = M;

    t = PROCS*M;
#if PROFILE_GENERAL
    time_select = 0.0;
    time_lb = 0.0;
#endif    
    result = all_select_c_alloc_d(M, A, N, t, (t>>1) + 1,sublist);
#if PROFILE_GENERAL
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYPROC==0) 
	fprintf(outfile,"Select_c general profile: Select: %9.6f  LB: %9.6f\n",
		time_select, time_lb);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    free(sublist);
    free(N);
    return (result);

#undef DATA_TYPE
}

int all_select_median_unbalanced_c_i(int M, int *A,
				       int *N, int total_n) {
    all_select_c_init();
    return all_select_c_i(M, A, N, total_n, (total_n>>1) + 1);

}

double all_select_median_unbalanced_c_d(int M, double *A,
					     int *N, int total_n) {
    all_select_c_init();
    return all_select_c_d(M, A, N, total_n, (total_n>>1) + 1);
}
