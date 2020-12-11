/*                                                                    tab:8
 *
 * select_seq.h - Sequential Selection Algorithm
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
 * Filename:            select_seq.h
 * History:
 */

#ifndef _SELECT_SEQ_H
#define _SELECT_SEQ_H

#include "main.h"

#define partition(a,b,c) partition_i(a,b,c)

#define partition_unk_piv(a,b,c) partition_unk_piv_i(a,b,c)
#if 0
int partition_unk_piv_i(int *, int, int);
int partition_unk_piv_d(double *, int, double);
#endif

#define select_mom(a,b,c) select_mom_i(a,b,c)
int    select_mom_i(int *, int, int);
double select_mom_d(double *, int, int);

#define select_mom_alloc(a,b,c,d) select_mom_alloc_i(a,b,c,d)
int    select_mom_alloc_i(int *, int, int, int *);
double select_mom_alloc_d(double *, int, int, double *);

int select_merge_i(int *A, int i, int bin_size, int bins);
double select_merge_d(double *A, int i, int bin_size, int bins);

void sort_min_i(int *A, int n, int X);
void sort_max_i(int *A, int n, int X);
void sort_min_d(double *A, int n, int X);
void sort_max_d(double *A, int n, int X);

_INLINE int partition_i(int *A, int n, int piv) {

#define DATA_TYPE int

    register int done = 0;
    
    register DATA_TYPE
	*d1, *d2,
	item;
    
    d1 = A;
    d2 = A + n-1;
    do {
	while (*d1 <= piv) d1++;
	while (*d2 >  piv) d2--;
	if (d1 < d2) {
	    item = *d1;
	    *d1 = *d2;
	    *d2 = item;
	}
	else done = 1;
    }
    while (!done);
    
    return((d2-A) + 1);

#undef DATA_TYPE
}

_INLINE int partition_d(double *A, int n, double piv) {

#define DATA_TYPE double

    register int done = 0;
    
    register DATA_TYPE
	*d1, *d2,
	item;
    
    d1 = A;
    d2 = A + n-1;
    do {
	while (*d1 <= piv) d1++;
	while (*d2 >  piv) d2--;
	if (d1 < d2) {
	    item = *d1;
	    *d1 = *d2;
	    *d2 = item;
	}
	else done = 1;
    }
    while (!done);
    
    return((d2-A) + 1);

#undef DATA_TYPE
}

_INLINE int partition_unk_piv_i(int *A, int n, int piv) {

#define DATA_TYPE int

    register int done = 0;
    
    register DATA_TYPE
	*d1, *d2,
	item;
    
    d1 = A;
    d2 = A + n-1;
    do {
	while ((d1<=d2)&&(*d1 <= piv)) d1++;
	while ((d2>=d1)&&(*d2 >  piv)) d2--;
	if (d1 < d2) {
	    item = *d1;
	    *d1 = *d2;
	    *d2 = item;
	}
	else done = 1;
    }
    while (!done);

    return ((d2-A) + 1);

#undef DATA_TYPE
}

_INLINE int partition_unk_piv_d(double *A, int n, double piv) {

#define DATA_TYPE double

    register int done = 0;
    
    register DATA_TYPE
	*d1, *d2,
	item;
    
    d1 = A;
    d2 = A + n-1;
    do {
	while ((d1<=d2)&&(*d1 <= piv)) d1++;
	while ((d2>=d1)&&(*d2 >  piv)) d2--;
	if (d1 < d2) {
	    item = *d1;
	    *d1 = *d2;
	    *d2 = item;
	}
	else done = 1;
    }
    while (!done);
    
    return((d2-A) + 1);

#undef DATA_TYPE
}

#endif
