/*                                                                    tab:8
 *
 * select_c.h - Parallel Selection Algorithm
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
 * Filename:            select_c.h
 * History:
 */

#ifndef _SELECT_C_H
#define _SELECT_C_H

#include "main.h"
#include "sorting.h"

#define all_select_median_c(a,b) all_select_median_c_i(a,b)
int all_select_median_c_i(int M, int *A);
double all_select_median_c_d(int M, double *A);

#define all_select_median_unbalanced_c(a,b,c,d) \
        all_select_median_unbalanced_c_i(a,b,c,d) 
int all_select_median_unbalanced_c_i(int M, int *A,
				     int *N, int total_n);
double all_select_median_unbalanced_c_d(int M, double *A,
					int *N, int total_n);

#define all_select_c(a,b,c,d,e) all_select_c_i(a,b,c,d,e)
int all_select_c_i(int M, int *A,
		   int *N, int total_n, int i);
double all_select_c_d(int M, double *A,
		      int *N, int total_n, int i);

#define all_select_c_alloc(a,b,c,d,e,f) all_select_c_alloc_i(a,b,c,d,e,f)
int all_select_c_alloc_i(int M, int *A,
			 int *N, int total_n, int i,
			 int* sublist);
double all_select_c_alloc_d(int M, double *A,
			    int *N, int total_n, int i,
			    double* sublist);

int all_select_max_alloc_i(int M, int *A,
			   int *N, int total_n, int i,
			   int *sublist, int X);
double all_select_min_alloc_d(int M, double *A,
			      int *N, int total_n, int i,
			      double *sublist, int X);
double all_select_max_alloc_d(int M, double *A,
			   int *N, int total_n, int i,
			   double *sublist, int X);
int all_select_min_alloc_i(int M, int *A,
			   int *N, int total_n, int i,
			   int *sublist, int X);

#endif
