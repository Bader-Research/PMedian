/*									tab:8
 *
 * sorting.h - fast sequential sorting routines
 *
 * 
 * "Copyright (c) 1994,1995 The Regents of the University of Maryland.
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
 * Authors: 		David A. Bader   <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:		1.0
 * Creation Date:	October 20, 1994
 * Filename:		sorting.h
 * History:
 */


#ifndef _SORTING_H
#define _SORTING_H

#include "main.h"

#define RADIXSORT_INT_BREAKPT 100

int intcompare(int *, int *);
void seq_radixsort(int *, int);
void seq_radixsort_19(int *, int);

_INLINE unsigned bits(unsigned x, int k, int j) {
/* Return the j bits which appear k bits from the right in x */
    return (x>>k) & ~(~0<<j);
}

#define insertsort(a,b) insertsort_i(a,b)

_INLINE void insertsort_i(int *A, int n) {

#define DATA_TYPE int
    
    register DATA_TYPE item;
    register int i,j;
    
    for (i=1 ; i<n ; i++) {
	item = A[i];
	j = i-1;
	while ((j>=0)&&(item < A[j])) {
	    A[j+1] = A[j];
	    j--;
	}
	A[j+1] = item;
    }

#undef DATA_TYPE    
}

_INLINE void insertsort_d(double *A, int n) {

#define DATA_TYPE double
    
    register DATA_TYPE item;
    register int i,j;
    
    for (i=1 ; i<n ; i++) {
	item = A[i];
	j = i-1;
	while ((j>=0)&&(item < A[j])) {
	    A[j+1] = A[j];
	    j--;
	}
	A[j+1] = item;
    }

#undef DATA_TYPE    
}


_INLINE void fastsort(int* arr, int nel) {

    if (nel>=RADIXSORT_INT_BREAKPT)
	seq_radixsort(arr,nel); 
    else
	insertsort(arr,nel);
}

_INLINE void fastsort_19(int* arr, int nel) {

    if (nel>=RADIXSORT_INT_BREAKPT)
	seq_radixsort_19(arr,nel); 
    else
	insertsort_i(arr,nel);
}



#endif


