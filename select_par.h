/*                                                                    tab:8
 *
 * select_par.h - Parallel Selection Algorithm
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
 * Filename:            select_par.h
 * History:
 */

#ifndef _SELECT_PAR_H
#define _SELECT_PAR_H

#include "main.h"
#include "sorting.h"

#define DEBUG_K 0

#define USE_MY_ALLTOALLV 0

void all_randomize_i(int *A, int *B, int M, int bin_size, int bins,
		     int logbins, int *bin_counts);
void all_randomize_d(double *A, double *B, int M, int bin_size, int bins,
		     int logbins, int *bin_counts);

int all_transpose_i(int *A, int *B, int bin_size, int *bin_counts);
int all_transpose_d(double *A, double *B, int bin_size, int *bin_counts);

void all_gather_picks_i(int    *A, int A_size, int    **B, int *B_size);
void all_gather_picks_d(double *A, int A_size, double **B, int *B_size);

int    all_select_gatherv_i(int    *A, int ni, int total_n, int i);
double all_select_gatherv_d(double *A, int ni, int total_n, int i);

int    all_select_gatherv_alloc_i(int    *A, int ni, int total_n, int i,
				  int    *B);
double all_select_gatherv_alloc_d(double *A, int ni, int total_n, int i,
				  double *B);

void partition_with_two_i(int *A, int *B, int M, int s0, int s1,
			  int *a0, int *a1);
void partition_with_two_d(double *A, double *B, int M, double s0, double s1,
			  int *a0, int *a1);

int all_select_min_set_i(int ni, int *A, int total_n, int i,
			   int *sublist, int X);
int all_select_max_set_i(int ni, int *A, int total_n, int i,
			 int *sublist, int X);
double all_select_min_set_d(int ni, double *A, int total_n, int i,
			      double *sublist, int X);
double all_select_max_set_d(int ni, double *A, int total_n, int i,
			   double *sublist, int X);

#endif
