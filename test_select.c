/*                                                                    tab:8
 *
 * test_select.c - Selection testing routines
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
 * Creation Date:       May 2, 1995
 * Filename:            test_select.c
 * History:
 */

#include "test_select.h"
#include "select_a.h"
#include "select_b.h"
#include "select_c.h"
#include "select_r.h"
#include "select_r2.h"
#include "select_r3.h"
#include "select_r4.h"
#include "select_r5.h"
#include "select_fr.h"
#include "select_ft.h"
#include "inputs_select.h"
#include "math_func.h"

#define TEST_MACH_LIMIT_LOG 12 /* 20 */ 
#define TEST_MACH_LIMIT (1<<TEST_MACH_LIMIT_LOG)
#define LOG_MINSELECT 12
#define NAS_REPEAT    3
#define WARM_UP_COUNT 4
#define WARM_UP       0

#define DO_NAS        0
#define DO_DOUBLES    0

#define DO_A          0
#define DO_B          0
#define DO_C          0
#define DO_R          0
#define DO_R2         1
#define DO_R3         1
#define DO_R4         0
#define DO_R5         1
#define DO_FR         1
#define DO_FT         1

#define DO_I_I        1
#define DO_I_S        1
#define DO_I_R        1
#define DO_I_N        1
#define DO_I_K        1

#define SEL_A         1
#define SEL_B         2
#define SEL_C         3
#define SEL_R         4
#define SEL_R2        5
#define SEL_R3        6
#define SEL_R4        7
#define SEL_R5        8
#define SEL_FR        9
#define SEL_FT       10

#define INP_I_W       0
#define INP_I_I       1
#define INP_I_S       2
#define INP_I_R       3
#define INP_I_N       4
#define INP_I_K       5

#define INP_D_W       0
#define INP_D_R       1
#define INP_D_G       2
#define INP_D_2G      3
#define INP_D_B       4
#define INP_D_S       5
#define INP_D_Z       6


void test_sel_i(int i, int* A, int alg, int inp) {
  double secs;
  int  result;
  char str_alg[3];
  char str_inp[2];
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  switch(inp) {
  case INP_I_I:  fill_same(i, A);    break;
  case INP_I_S:  fill_linear(i, A);  break;
  case INP_I_W:
  case INP_I_R:  fill_random(i, A);  break;
  case INP_I_N:  all_fill_nas(i, A); break;
  case INP_I_K:  all_fill_skewed(i, A); break;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  secs = MPI_Wtime();
  switch(alg) {
  case SEL_A:  result = all_select_median_a(i, A);  break;
  case SEL_B:  result = all_select_median_b(i, A);  break;
  case SEL_C:  result = all_select_median_c(i, A);  break;
  case SEL_R:  result = all_select_median_r(i, A);  break;
  case SEL_R2: result = all_select_median_r2(i, A); break;
  case SEL_R3: result = all_select_median_r3(i, A); break;
  case SEL_R4: result = all_select_median_r4(i, A); break;
  case SEL_R5: result = all_select_median_r5(i, A); break;
  case SEL_FR: result = all_select_median_fr(i, A); break;
  case SEL_FT: result = all_select_median_ft(i, A); break;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  secs = MPI_Wtime() - secs;

  if (MYPROC==0) {
    switch(alg) {
    case SEL_A:  sprintf(str_alg,"_a"); break;
    case SEL_B:  sprintf(str_alg,"_b"); break;
    case SEL_C:  sprintf(str_alg,"_c"); break;
    case SEL_R:  sprintf(str_alg,"_r"); break;
    case SEL_R2: sprintf(str_alg,"r2"); break;
    case SEL_R3: sprintf(str_alg,"r3"); break;
    case SEL_R4: sprintf(str_alg,"r4"); break;
    case SEL_R5: sprintf(str_alg,"r5"); break;
    case SEL_FR: sprintf(str_alg,"fr"); break;
    case SEL_FT: sprintf(str_alg,"ft"); break;
    }
    switch(inp) {
    case INP_I_I:  sprintf(str_inp,"I");   break;
    case INP_I_S:  sprintf(str_inp,"S");   break;
    case INP_I_R:  sprintf(str_inp,"R");   break;
    case INP_I_W:  sprintf(str_inp,"W");   break;
    case INP_I_N:  sprintf(str_inp,"N");   break;
    case INP_I_K:  sprintf(str_inp,"K");   break;
    }
    if (inp != INP_I_W) { /* WARM UP */
      fprintf(outfile,"Median_%2s %4d %12d [%1s] i %12d \t  %9.6f\n",
	      str_alg,PROCS, i, str_inp, result, secs);
      fflush(outfile);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return;
}

void test_sel_d(int i, double* A, int alg, int inp) {
  double secs;
  double result;
  char str_alg[3];
  char str_inp[4];
  
  switch(inp) {
  case INP_D_W:
  case INP_D_R:  fill_random_d(i, A);     break;
  case INP_D_G:  fill_gaussian_d(i, A);   break;
  case INP_D_2G: fill_g_group_d(i, A, 2); break;
  case INP_D_B:  fill_bucket_d(i, A);     break;
  case INP_D_S:  fill_staggered_d(i, A);  break;
  case INP_D_Z:  fill_zero_d(i, A);       break;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  secs = MPI_Wtime();
  switch(alg) {
  case SEL_C:  result = all_select_median_c_d(i, A);  break;
  case SEL_R:  result = all_select_median_r_d(i, A);  break;
  case SEL_R2: result = all_select_median_r2_d(i, A); break;
  case SEL_R3: result = all_select_median_r3_d(i, A); break;
  case SEL_R4: result = all_select_median_r4_d(i, A); break;
  case SEL_R5: result = all_select_median_r5_d(i, A); break;
  case SEL_FR: result = all_select_median_fr_d(i, A); break;
  case SEL_FT: result = all_select_median_ft_d(i, A); break;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  secs = MPI_Wtime() - secs;

  if (MYPROC==0) {
    switch(alg) {
    case SEL_C:  sprintf(str_alg,"c "); break;
    case SEL_R:  sprintf(str_alg,"r "); break;
    case SEL_R2: sprintf(str_alg,"r2"); break;
    case SEL_R3: sprintf(str_alg,"r3"); break;
    case SEL_R4: sprintf(str_alg,"r4"); break;
    case SEL_R5: sprintf(str_alg,"r5"); break;
    case SEL_FR: sprintf(str_alg,"fr"); break;
    case SEL_FT: sprintf(str_alg,"ft"); break;
    }
    switch(inp) {
    case INP_D_R:  sprintf(str_inp,"R  "); break;
    case INP_D_G:  sprintf(str_inp,"G  "); break;
    case INP_D_2G: sprintf(str_inp,"2-G"); break;
    case INP_D_B:  sprintf(str_inp,"B  "); break;
    case INP_D_S:  sprintf(str_inp,"S  "); break;
    case INP_D_Z:  sprintf(str_inp,"Z  "); break;
    case INP_D_W:  sprintf(str_inp,"W  "); break;
    }
    if (inp != INP_D_W) { /* WARM UP */
      fprintf(outfile,"Median_%2s %4d %12d [%3s] d %15g \t  %9.6f\n",
	      str_alg,PROCS, i, str_inp, result, secs);
      fflush(outfile);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return;
}

void all_test_select_d()
{
    register int
	i, j;

    double result;

    double *A;

    double secs;

    A = (double *)malloc(TEST_MACH_LIMIT*sizeof(double));
    assert_malloc(A);

    RNG_init();
	
#if WARM_UP
    for (i=0 ; i<WARM_UP_COUNT ; i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      j = TEST_MACH_LIMIT;
      fill_random_d(j, A);
      MPI_Barrier(MPI_COMM_WORLD);
      secs = MPI_Wtime();
#if DO_C
      result = all_select_median_c_d(j, A); 
#endif
      MPI_Barrier(MPI_COMM_WORLD);
      secs = MPI_Wtime() - secs;
      if (MYPROC==0) {
	fprintf(outfile,"dis%d: %9.6f %g\n",i,secs,result);
	fflush(outfile);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

#if 1
    for (i = (1<<LOG_MINSELECT) ; i <= TEST_MACH_LIMIT ; i = (i<<1) ) {

	MPI_Barrier(MPI_COMM_WORLD);

#if DO_C
/* Test median random double */
	test_sel_d(i, A, SEL_C, INP_D_R);
#endif
	
#if DO_R
/* Test median random double */
	test_sel_d(i, A, SEL_R, INP_D_R);
#endif

#if DO_R2
/* Test median random double */
	test_sel_d(i, A, SEL_R2, INP_D_R);
#endif

#if DO_R3
/* Test median random double */
	test_sel_d(i, A, SEL_R3, INP_D_R);
#endif

#if DO_R4
/* Test median random double */
	test_sel_d(i, A, SEL_R4, INP_D_R);
#endif

#if DO_R5
/* Test median random double */
	test_sel_d(i, A, SEL_R5, INP_D_R);
#endif

#if DO_FR
/* Test median random double */
	test_sel_d(i, A, SEL_FR, INP_D_R);
#endif

#if DO_FT
/* Test median random double */
	test_sel_d(i, A, SEL_FT, INP_D_R);
#endif

#if DO_C
/* Test median gaussian double */
	test_sel_d(i, A, SEL_C, INP_D_G);
#endif
	
#if DO_R
/* Test median gaussian double */
	test_sel_d(i, A, SEL_R, INP_D_G);
#endif

#if DO_R2
/* Test median gaussian double */
	test_sel_d(i, A, SEL_R2, INP_D_G);
#endif

#if DO_R3
/* Test median gaussian double */
	test_sel_d(i, A, SEL_R3, INP_D_G);
#endif

#if DO_R4
/* Test median gaussian double */
	test_sel_d(i, A, SEL_R4, INP_D_G);
#endif

#if DO_R5
/* Test median gaussian double */
	test_sel_d(i, A, SEL_R5, INP_D_G);
#endif

#if DO_FR
/* Test median gaussian double */
	test_sel_d(i, A, SEL_FR, INP_D_G);
#endif

#if DO_FT
/* Test median gaussian double */
	test_sel_d(i, A, SEL_FT, INP_D_G);
#endif

#if 0
#if DO_C
/* Test median 2-g double */
	test_sel_d(i, A, SEL_C, INP_D_2G);
#endif
#if DO_R
/* Test median 2-g double */
	test_sel_d(i, A, SEL_R, INP_D_2G);
#endif
#if DO_R2
/* Test median 2-g double */
	test_sel_d(i, A, SEL_R2, INP_D_2G);
#endif
#if DO_R3
/* Test median 2-g double */
	test_sel_d(i, A, SEL_R3, INP_D_2G);
#endif
#if DO_R4
/* Test median 2-g double */
	test_sel_d(i, A, SEL_R4, INP_D_2G);
#endif
#if DO_R5
/* Test median 2-g double */
	test_sel_d(i, A, SEL_R5, INP_D_2G);
#endif
#if DO_FR
/* Test median 2-g double */
	test_sel_d(i, A, SEL_FR, INP_D_2G);
#endif
#if DO_FT
/* Test median 2-g double */
	test_sel_d(i, A, SEL_FT, INP_D_2G);
#endif
#endif
	
#if 0
#if DO_C
/* Test median bucket double */
	test_sel_d(i, A, SEL_C, INP_D_B);
#endif
#if DO_R
/* Test median bucket double */
	test_sel_d(i, A, SEL_R, INP_D_B);
#endif
#if DO_R2
/* Test median bucket double */
	test_sel_d(i, A, SEL_R2, INP_D_B);
#endif
#if DO_R3
/* Test median bucket double */
	test_sel_d(i, A, SEL_R3, INP_D_B);
#endif
#if DO_R4
/* Test median bucket double */
	test_sel_d(i, A, SEL_R4, INP_D_B);
#endif
#if DO_R5
/* Test median bucket double */
	test_sel_d(i, A, SEL_R5, INP_D_B);
#endif
#if DO_FR
/* Test median bucket double */
	test_sel_d(i, A, SEL_FR, INP_D_B);
#endif
#if DO_FT
/* Test median bucket double */
	test_sel_d(i, A, SEL_FT, INP_D_B);
#endif
#endif

#if DO_C
/* Test median staggered double */
	test_sel_d(i, A, SEL_C, INP_D_S);
#endif
#if DO_R
/* Test median staggered double */
	test_sel_d(i, A, SEL_R, INP_D_S);
#endif
#if DO_R2
/* Test median staggered double */
	test_sel_d(i, A, SEL_R2, INP_D_S);
#endif
#if DO_R3
/* Test median staggered double */
	test_sel_d(i, A, SEL_R3, INP_D_S);
#endif
#if DO_R4
/* Test median staggered double */
	test_sel_d(i, A, SEL_R4, INP_D_S);
#endif
#if DO_R5
/* Test median staggered double */
	test_sel_d(i, A, SEL_R5, INP_D_S);
#endif
#if DO_FR
/* Test median staggered double */
	test_sel_d(i, A, SEL_FR, INP_D_S);
#endif
#if DO_FT
/* Test median staggered double */
	test_sel_d(i, A, SEL_FT, INP_D_S);
#endif

#if 0
#if DO_C
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_C, INP_D_Z);
#endif	
#if DO_R
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_R, INP_D_Z);
#endif	
#if DO_R2
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_R2, INP_D_Z);
#endif	
#if DO_R3
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_R3, INP_D_Z);
#endif	
#if DO_R4
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_R4, INP_D_Z);
#endif	
#if DO_R5
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_R5, INP_D_Z);
#endif	
#if DO_FR
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_FR, INP_D_Z);
#endif	
#if DO_FT
/* Test median zero entropy double */
	test_sel_d(i, A, SEL_FT, INP_D_Z);
#endif	
#endif	

    }
#endif
    
    free(A);
}


void all_test_select() 
{
    register int
	i, j;

    int n;
    int *N;
    int *A;

    double secs;

#if DO_DOUBLES
    all_test_select_d();
#endif

    N = (int *)malloc(PROCS*sizeof(int));
    assert_malloc(N);

    A = (int *)malloc(TEST_MACH_LIMIT*sizeof(int));
    assert_malloc(A);
    
    RNG_init();

#if WARM_UP
    for (j=0 ; j<WARM_UP_COUNT ; j++) {
      MPI_Barrier(MPI_COMM_WORLD);

      i = TEST_MACH_LIMIT;
      MPI_Barrier(MPI_COMM_WORLD);
#if DO_A
      test_sel_i(i, A, SEL_A, INP_I_W);
#endif
#if DO_B
      test_sel_i(i, A, SEL_B, INP_I_W);
#endif
#if DO_C
      test_sel_i(i, A, SEL_C, INP_I_W);
#endif
#if DO_R
      test_sel_i(i, A, SEL_R, INP_I_W);
#endif
#if DO_R2
      test_sel_i(i, A, SEL_R2, INP_I_W);
#endif
#if DO_R3
      test_sel_i(i, A, SEL_R3, INP_I_W);
#endif
#if DO_R4
      test_sel_i(i, A, SEL_R4, INP_I_W);
#endif
#if DO_R5
      test_sel_i(i, A, SEL_R5, INP_I_W);
#endif
    }




#endif

#if DO_NAS
    MPI_Barrier(MPI_COMM_WORLD);

    if (TEST_MACH_LIMIT_LOG + PROCSLOG >= 23) {
	int *B;

	i = DIVPROCS(1<<23);

	B = (int *)malloc(i*sizeof(int));
	assert_malloc(B);

	all_fill_nas(i, B);

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));

	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_B
	  j = all_select_median_b(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;

	  if (MYPROC==0) {
	    fprintf(outfile,"Median_b  %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));

	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_B
	  j = all_select_median_b_nas(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;

	  if (MYPROC==0) {
	    fprintf(outfile,"Median_b_nas  %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_C
	  j = all_select_median_c(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_c  %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_R
	  j = all_select_median_r(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_r  %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_R2
	  j = all_select_median_r2(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_r2 %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_R3
	  j = all_select_median_r3(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_r3 %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_R4
	  j = all_select_median_r4(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_r4 %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_R5
	  j = all_select_median_r5(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_r5 %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_FR
	  j = all_select_median_FR(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_fr %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}
	
	for (n=0 ; n<NAS_REPEAT ; n++) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  bcopy(B, A, i*sizeof(int));
	
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime();
#if DO_FT
	  j = all_select_median_FT(i, A);
#endif
	  MPI_Barrier(MPI_COMM_WORLD);
	  secs = MPI_Wtime() - secs;
	  
	  if (MYPROC==0) {
	    fprintf(outfile,"Median_ft %4d %12d NAS  %9d  %9.6f\n",
		    PROCS, i, j,secs);
	    fflush(outfile);
	  }
	}

	free(B);
      }
    MPI_Barrier(MPI_COMM_WORLD);
#endif    

#if 1
    for (i = (1<<LOG_MINSELECT) ; i <= TEST_MACH_LIMIT ; i = (i<<1) ) {

	MPI_Barrier(MPI_COMM_WORLD);

#if DO_I_I
#if DO_A
/* Test median same */
	test_sel_i(i, A, SEL_A, INP_I_I);
#endif

#if DO_B
/* Test median same */
	test_sel_i(i, A, SEL_B, INP_I_I);
#endif

#if DO_C
/* Test median same */
	test_sel_i(i, A, SEL_C, INP_I_I);
#endif
	
#if DO_R
/* Test median same */
	test_sel_i(i, A, SEL_R, INP_I_I);
#endif

#if DO_R2
/* Test median same */
	test_sel_i(i, A, SEL_R2, INP_I_I);
#endif

#if DO_R3
/* Test median same */
	test_sel_i(i, A, SEL_R3, INP_I_I);
#endif

#if DO_R4
/* Test median same */
	test_sel_i(i, A, SEL_R4, INP_I_I);
#endif

#if DO_R5
/* Test median same */
	test_sel_i(i, A, SEL_R5, INP_I_I);
#endif

#if DO_FR
/* Test median same */
	test_sel_i(i, A, SEL_FR, INP_I_I);
#endif

#if DO_FT
/* Test median same */
	test_sel_i(i, A, SEL_FT, INP_I_I);
#endif
#endif

#if DO_I_S	
#if DO_A
/* Test median linear */
	test_sel_i(i, A, SEL_A, INP_I_S);
#endif

#if DO_B
/* Test median linear */
	test_sel_i(i, A, SEL_B, INP_I_S);
#endif

#if DO_C
/* Test median linear */
	test_sel_i(i, A, SEL_C, INP_I_S);
#endif
	
#if DO_R
/* Test median linear */
	test_sel_i(i, A, SEL_R, INP_I_S);
#endif

#if DO_R2
/* Test median linear */
	test_sel_i(i, A, SEL_R2, INP_I_S);
#endif

#if DO_R3
/* Test median linear */
	test_sel_i(i, A, SEL_R3, INP_I_S);
#endif

#if DO_R4
/* Test median linear */
	test_sel_i(i, A, SEL_R4, INP_I_S);
#endif

#if DO_R5
/* Test median linear */
	test_sel_i(i, A, SEL_R5, INP_I_S);
#endif

#if DO_FR
/* Test median linear */
	test_sel_i(i, A, SEL_FR, INP_I_S);
#endif

#if DO_FT
/* Test median linear */
	test_sel_i(i, A, SEL_FT, INP_I_S);
#endif
#endif

#if DO_I_R
#if DO_A
/* Test median random */
	test_sel_i(i, A, SEL_A, INP_I_R);
#endif

#if DO_B
/* Test median random */
	test_sel_i(i, A, SEL_B, INP_I_R);
#endif

#if DO_C
/* Test median random */
	test_sel_i(i, A, SEL_C, INP_I_R);
#endif
	
#if DO_R
/* Test median random */
	test_sel_i(i, A, SEL_R, INP_I_R);
#endif
	
#if DO_R2
/* Test median random */
	test_sel_i(i, A, SEL_R2, INP_I_R);
#endif

#if DO_R3
/* Test median random */
	test_sel_i(i, A, SEL_R3, INP_I_R);
#endif

#if DO_R4
/* Test median random */
	test_sel_i(i, A, SEL_R4, INP_I_R);
#endif

#if DO_R5
/* Test median random */
	test_sel_i(i, A, SEL_R5, INP_I_R);
#endif

#if DO_FR
/* Test median random */
	test_sel_i(i, A, SEL_FR, INP_I_R);
#endif

#if DO_FT
/* Test median random */
	test_sel_i(i, A, SEL_FT, INP_I_R);
#endif
#endif

#if DO_I_N
#if DO_A
/* Test median NAS */
	test_sel_i(i, A, SEL_A, INP_I_N);
#endif

#if DO_B
/* Test median NAS */
	test_sel_i(i, A, SEL_B, INP_I_N);
#endif

#if DO_C
/* Test median NAS */
	test_sel_i(i, A, SEL_C, INP_I_N);
#endif
	
#if DO_R
/* Test median NAS */
	test_sel_i(i, A, SEL_R, INP_I_N);
#endif
	
#if DO_R2
/* Test median NAS */
	test_sel_i(i, A, SEL_R2, INP_I_N);
#endif

#if DO_R3
/* Test median NAS */
	test_sel_i(i, A, SEL_R3, INP_I_N);
#endif

#if DO_R4
/* Test median NAS */
	test_sel_i(i, A, SEL_R4, INP_I_N);
#endif

#if DO_R5
/* Test median NAS */
	test_sel_i(i, A, SEL_R5, INP_I_N);
#endif

#if DO_FR
/* Test median NAS */
	test_sel_i(i, A, SEL_FR, INP_I_N);
#endif

#if DO_FT
/* Test median NAS */
	test_sel_i(i, A, SEL_FT, INP_I_N);
#endif
#endif

#if DO_I_K
#if DO_A
/* Test median skewed */
	test_sel_i(i, A, SEL_A, INP_I_K);
#endif

#if DO_B
/* Test median skewed */
	test_sel_i(i, A, SEL_B, INP_I_K);
#endif

#if DO_C
/* Test median skewed */
	test_sel_i(i, A, SEL_C, INP_I_K);
#endif
	
#if DO_R
/* Test median skewed */
	test_sel_i(i, A, SEL_R, INP_I_K);
#endif
	
#if DO_R2
/* Test median skewed */
	test_sel_i(i, A, SEL_R2, INP_I_K);
#endif

#if DO_R3
/* Test median skewed */
	test_sel_i(i, A, SEL_R3, INP_I_K);
#endif

#if DO_R4
/* Test median skewed */
	test_sel_i(i, A, SEL_R4, INP_I_K);
#endif

#if DO_R5
/* Test median skewed */
	test_sel_i(i, A, SEL_R5, INP_I_K);
#endif

#if DO_FR
/* Test median skewed */
	test_sel_i(i, A, SEL_FR, INP_I_K);
#endif

#if DO_FT
/* Test median skewed */
	test_sel_i(i, A, SEL_FT, INP_I_K);
#endif
#endif
	
    }
#endif

    free(N);
    free(A);
}


