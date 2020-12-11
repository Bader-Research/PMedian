/*                                                                    tab:8
 *
 * select_randsamp.c - Parallel Selection Algorithm
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
 * Filename:            select_randsamp.c
 * History:
 */

#include "select_randsamp.h"

#define DEBUG         0

#define BUNDLE        4


int all_random_pick_exact_i(int *A, int *B, int M, int n) {

  double mscale;
  int *B_ptr, *E_ptr;
  int i, j, k;

  if (M>0) {
    k = (int)ceil((double)M * pow((double)n, PICK_EPS-1.0));
    B_ptr = B;
    E_ptr = B_ptr + k;
    mscale = (double)M/2147483648.0;
    
    /*     for (i=0 ; i<k ; i++) */
    while (B_ptr < E_ptr) 
    {
      j = (int)floor((double)rrandom()*mscale);
#if 0
      if ((j >= M) || (j < 0))
	fprintf(stderr,"ERROR: all_random_pick_exact_i  j: %d M: %d\n",j,M);
#endif
      *B_ptr++ = A[j];
    }
  }
  else
    k = 0;

  return (k);
}

int all_random_pick_exact_d(double *A, double *B, int M, int n) {

  double mscale;
  double *B_ptr, *E_ptr;
  int i, j, k;

  if (M>0) {
    B_ptr = B;
    k = (int)ceil((double)M * pow((double)n, PICK_EPS-1.0));
    mscale = (double)M/2147483648.0;
    
    /* for (i=0 ; i<k ; i++) */
    E_ptr = B_ptr + k;
    while (B_ptr < E_ptr)
    {
      j = (int)floor((double)rrandom()*mscale);
#if 0
      if ((j >= M) || (j < 0))
	fprintf(stderr,"ERROR: all_random_pick_exact_d  j: %d M: %d\n",j,M);
#endif
      *B_ptr++ = A[j];
    }
  }
  else
    k = 0;

  return (k);
}


int all_random_pick_exactbund_i(int *A, int *B, int M, int n) {

  double mscale;
  int *B_ptr;
  int i, j, k;
  int times, rem, ebs;

  if (M>0) {
    k = (int)ceil((double)M * pow((double)n, PICK_EPS-1.0));

    times = k / BUNDLE;
    rem   = k - (times*BUNDLE);
    ebs   = BUNDLE*sizeof(int);

    B_ptr = B;
    mscale = (double)M/2147483648.0;

    for (i=0 ; i<times ; i++) {
      j = (int)floor((double)rrandom()*mscale);
      if (j+BUNDLE-1 >= M) {
	j = M - BUNDLE;
#if 0
	fprintf(stderr,"ERROR: all_random_pick_exact_i  j: %d M: %d\n",j,M);
#endif
      }
      memcpy(B_ptr,  A+j, ebs);
      B_ptr += BUNDLE;
    }

    if (rem>0) {
      j = (int)floor((double)rrandom()*mscale);
      if (j+BUNDLE-1 >= M) {
	j = M - BUNDLE;
#if 0
	fprintf(stderr,"ERROR: all_random_pick_exact_i  j: %d M: %d\n",j,M);
#endif
      }
      memcpy(B_ptr,  A+j, rem*sizeof(int));
    }

#if 0
    if (k != B_ptr+rem - B)
      fprintf(stderr,"ERROR: k: %d  B_ptr-rem-B: %d\n",k, B_ptr+rem-B);
#endif

  }
  else
    k = 0;

  
  return (k);
}

int all_random_pick_exactbund_d(double *A, double *B, int M, int n) {

  double mscale;
  double *B_ptr;
  int i, j, k;
  int times, rem, ebs;

  if (M>0) {
    k = (int)ceil((double)M * pow((double)n, PICK_EPS-1.0));

    times = k / BUNDLE;
    rem   = k - (times*BUNDLE);
    ebs   = BUNDLE*sizeof(double);

    B_ptr = B;
    mscale = (double)M/2147483648.0;

    for (i=0 ; i<times ; i++) {
      j = (int)floor((double)rrandom()*mscale);
      if (j+BUNDLE-1 >= M) {
	j = M - BUNDLE;
#if 0
	fprintf(stderr,"ERROR: all_random_pick_exact_d  j: %d M: %d\n",j,M);
#endif
      }
      memcpy(B_ptr,  A+j, ebs);
      B_ptr += BUNDLE;
    }

    if (rem>0) {
      j = (int)floor((double)rrandom()*mscale);
      if (j+BUNDLE-1 >= M) {
	j = M - BUNDLE;
#if 0
	fprintf(stderr,"ERROR: all_random_pick_exact_d  j: %d M: %d\n",j,M);
#endif
      }
      memcpy(B_ptr,  A+j, rem*sizeof(double));
    }

#if 0
    if (k != B_ptr+rem - B)
      fprintf(stderr,"ERROR: k: %d  B_ptr-rem-B: %d\n",k, B_ptr+rem-B);
#endif

  }
  else
    k = 0;

  
  return (k);
}
