/*                                                                    tab:8
 *
 * select_b.h - Parallel Selection Algorithm
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
 * Filename:            select_b.h
 * History:
 */

#ifndef _SELECT_B_H
#define _SELECT_B_H

#include "main.h"
#include "sorting.h"

int all_select_median_b(int M, int *A);

int all_select_median_unbalanced_b(int M, int *A,
				   int *N, int total_n);

int all_select_b(int M, int *A,
		 int *N, int total_n, int i);


#define all_LoadBalance_b(a,b,c,d) all_LoadBalance_b_i(a,b,c,d)

void all_LoadBalance_b_i(int M, int *A,
			 int *N, int total_n);

void all_LoadBalance_b_d(int M, double *A,
			 int *N, int total_n);

#endif
