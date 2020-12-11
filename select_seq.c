/*                                                                    tab:8
 *
 * select_seq.c - Sequential Selection Algorithm
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
 * Filename:            select_seq.c
 * History:
 */

#include "select_seq.h"
#include "sorting.h"
#include <math.h>

#define DEFAULT_R       7
#define SELECT_CUTOFF  22

#define DEBUG 1

int select_mom_i(int* list, int n, int k) {

#define DATA_TYPE int
    
    register DATA_TYPE
	*s1, *s2;

    register
	blocks,
	left;

    DATA_TYPE
	result,
	*sublist,
	*listptr,
	*listptr_med;

    if (k==0) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = min(result,*list);
	list++;
      }
      return(result);
    }

    if (k==n) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = max(result,*list);
	list++;
      }
      return(result);
    }
    
    if (n < SELECT_CUTOFF) {
	insertsort_i(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;
    
    sublist = (DATA_TYPE *)malloc(blocks*sizeof(DATA_TYPE));
    assert_malloc(sublist);
    
    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_i(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_i(sublist,blocks,blocks/2);
    
    free(sublist);
    
    left = partition_i(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_i(list,left,k);
	else
	    if (k > left) 
		result = select_mom_i(list+left,n-left,k-left);
    }
    else {
	insertsort_i(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}

double select_mom_d(double* list, int n, int k) {

#define DATA_TYPE double
    
    register DATA_TYPE
	*s1, *s2;

    register int
	blocks,
	left;

    DATA_TYPE
	result,
	*sublist,
	*listptr,
	*listptr_med;
    
    if (k==0) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = min(result,*list);
	list++;
      }
      return(result);
    }

    if (k==n) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = max(result,*list);
	list++;
      }
      return(result);
    }
    
    if (n < SELECT_CUTOFF) {
	insertsort_d(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;

    sublist = (DATA_TYPE *)malloc(blocks*sizeof(DATA_TYPE));
    assert_malloc(sublist);

    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_d(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_d(sublist,blocks,blocks/2);
    
    free(sublist);
    
    left = partition_d(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_d(list,left,k);
	else
	    if (k > left) 
		result = select_mom_d(list+left,n-left,k-left);
    }
    else {
	insertsort_d(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}


int select_mom_alloc_i(int* list, int n, int k, int* sublist) {

#define DATA_TYPE int
    
    register DATA_TYPE
	*s1, *s2;

    register int
	blocks,
	left;

    DATA_TYPE
	result,
	*listptr,
	*listptr_med;
    
    if (k==0) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = min(result,*list);
	list++;
      }
      return(result);
    }

    if (k==n) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = max(result,*list);
	list++;
      }
      return(result);
    }
    
    if (n < SELECT_CUTOFF) {
	insertsort_i(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;
    
    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_i(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_alloc_i(sublist,blocks,blocks/2,s2);
    
    left = partition_i(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_alloc_i(list,left,k,sublist);
	else
	    if (k > left) 
		result = select_mom_alloc_i(list+left,n-left,k-left,sublist);
    }
    else {
	insertsort_i(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}

double select_mom_alloc_d(double* list, int n, int k, double* sublist) {

#define DATA_TYPE double
    
    register DATA_TYPE
	*s1, *s2;

    register int
	blocks,
	left;

    DATA_TYPE
	result,
	*listptr,
	*listptr_med;
    
    if (k==0) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = min(result,*list);
	list++;
      }
      return(result);
    }

    if (k==n) {
      s1 = list+n;
      result = *list;
      list++;
      while (list < s1) {
	result = max(result,*list);
	list++;
      }
      return(result);
    }
    
    if (n < SELECT_CUTOFF) {
	insertsort_d(list,n);
	return(list[k-1]);
    }
    
    blocks = n/DEFAULT_R;
    
    listptr = list;
    listptr_med = list + DEFAULT_R/2;
    s1 = sublist;
    s2 = sublist + blocks;
    while (s1 < s2) {
	insertsort_d(listptr, DEFAULT_R);
	listptr += DEFAULT_R;
	*s1++ = *listptr_med;
	listptr_med += DEFAULT_R;
    }

    result = select_mom_alloc_d(sublist,blocks,blocks/2,s2);
    
    left = partition_d(list,n,result);
    
    if (left<n) {
	if (k < left) 
	    result = select_mom_alloc_d(list,left,k,sublist);
	else
	    if (k > left) 
		result = select_mom_alloc_d(list+left,n-left,k-left,sublist);
    }
    else {
	insertsort_d(list,n);
	result = list[k-1];
    }
    
    return (result);

#undef DATA_TYPE
}

void sort_min_i(int *A, int n, int X) {

  register int melem, elem;
  register int i,j;

#if DEBUG
  if (X > n) {
    fprintf(stderr,"PE%3d: ERROR: in sort_min_i() [X>n] X: %d n: %d\n",
	    MYPROC, X, n);
  }
#endif
  
  for (i=1 ; i<X ; i++) {
    j = i-1;
    elem = A[i];
    while ((j>=0)&&(elem < A[j])) {
      A[j+1] = A[j];
      j--;
    }
    A[j+1] = elem;
  }

  melem = A[X-1];
  
  for (i=X ; i<n ; i++) {
    elem = A[i];
    if (elem < melem) {
      j = X-1;
      while ((j>=0)&&(elem < A[j])) {
	A[j+1] = A[j];
	j--;
      }
      A[j+1] = elem;
      melem  = A[X-1];
    }
  }
}

void sort_max_i(int *A, int n, int X) {

  register int melem, elem;
  register int i,j;
   
#if DEBUG
  if (X > n) {
    fprintf(stderr,"PE%3d: ERROR: in sort_max_i() [X>n] X: %d n: %d\n",
	    MYPROC, X, n);
  }
#endif
  
  for (i=n-2 ; i>n-X-1 ; i--) {
    j = i+1;
    elem = A[i];
    while ((j<n)&&(elem > A[j])) {
      A[j-1] = A[j];
      j++;
    }
    A[j-1] = elem;
  }

  melem = A[n-X];
  
  for (i=n-X-1 ; i>=0 ; i--) {
    elem = A[i];
    if (elem > melem) {
      j = n-X;
      while ((j<n)&&(elem > A[j])) {
	A[j-1] = A[j];
	j++;
      }
      A[j-1] = elem;
      melem  = A[n-X];
    }
  }
}


void sort_min_d(double *A, int n, int X) {

  register double melem, elem;
  register int i,j;
   
#if DEBUG
  if (X > n) {
    fprintf(stderr,"ERROR: in sort_min_d()\n");
  }
#endif
  
  for (i=1 ; i<X ; i++) {
    j = i-1;
    elem = A[i];
    while ((j>=0)&&(elem < A[j])) {
      A[j+1] = A[j];
      j--;
    }
    A[j+1] = elem;
  }

  melem = A[X-1];
  
  for (i=X ; i<n; i++) {
    elem = A[i];
    if (elem < melem) {
      j = X-1;
      while ((j>=0)&&(elem < A[j])) {
	A[j+1] = A[j];
	j--;
      }
      A[j+1] = elem;
      melem  = A[X-1];
    }
  }
}

void sort_max_d(double *A, int n, int X) {

  register double melem, elem;
  register int i,j;
   
#if DEBUG
  if (X > n) {
    fprintf(stderr,"ERROR: in sort_max_d()\n");
  }
#endif
  
  for (i=n-2 ; i>n-X-1 ; i--) {
    j = i+1;
    elem = A[i];
    while ((j<n)&&(elem > A[j])) {
      A[j-1] = A[j];
      j++;
    }
    A[j-1] = elem;
  }

  melem = A[n-X];
  
  for (i=n-X-1 ; i>=0 ; i--) {
    elem = A[i];
    if (elem > melem) {
      j = n-X;
      while ((j<n)&&(elem > A[j])) {
	A[j-1] = A[j];
	j++;
      }
      A[j-1] = elem;
      melem  = A[n-X];
    }
  }
}



/************************* p-way merge ***********************/

typedef struct i_record {
  int val;
  int start;
} i_record_t;

int select_merge_pway_i(int *A, int idx, int bin_size, int bins) {
  /* Assume bins is a power of two */

  int result;
  register int v0, v1;
  int c_val, t_val;
  int **Source, **source_ptr;
  
  int i,j,k,
    c_start,t_start;
  i_record_t *root_ptr,*tree_ptr, *Tree, *TempT;

  Tree = (i_record_t *)malloc(bins*sizeof(i_record_t));
  assert_malloc(Tree); 

  TempT = (i_record_t *)malloc(bins*sizeof(i_record_t));
  assert_malloc(TempT); 

  Source = (int **)malloc(bins*sizeof(int *));
  assert_malloc(Source); 
   
   
  Source[0] = A;
  for (i=1 ; i<bins ; i++) 
    Source[i] = Source[i-1] + bin_size;
   
  source_ptr = Source;
  j          = 0;
  for (i=(bins>>1) ; i<bins ; i++) {
    v0 = **source_ptr;
    v1 = **(source_ptr + 1);
    if (v1 > v0) {
      Tree[i].val    = v1;
      Tree[i].start  = j + 1;
      TempT[i].val   = v0;
      TempT[i].start = j;
    }
    else {
      Tree[i].val    = v0;
      Tree[i].start  = j;
      TempT[i].val   = v1;
      TempT[i].start = j + 1;
    }
    source_ptr    += 2;
    j             += 2;
  }
   
  for (i=(bins/2 - 1) ; i>0 ; i--) {
    j = (i<<1);
    k = j + 1;
    if (TempT[k].val > TempT[j].val) {
      Tree[i].val    = TempT[k].val;
      Tree[i].start  = TempT[k].start;
      TempT[i].val   = TempT[j].val;
      TempT[i].start = TempT[j].start;
    }
    else {
      Tree[i].val    = TempT[j].val;
      Tree[i].start  = TempT[j].start;
      TempT[i].val   = TempT[k].val;
      TempT[i].start = TempT[k].start;
    }
  }
   
  result          = TempT[1].val; 
  c_start         = TempT[1].start;
  root_ptr = Tree + 1;

  for (j=1 ; j<idx ; j++) {
#if DEBUG
    if ((Source[c_start] - (A + (c_start*bin_size))) > bin_size) {
      fprintf(stderr,"ERROR: Out of range\n");
    }
#endif
    c_val = *(++Source[c_start]); 
    i = c_start + bins;
    while (i > 3) {
      if (c_val > Tree[i>>=1].val) {
	t_val           = c_val;
	t_start         = c_start;
	tree_ptr        = Tree+i;
	c_val           = tree_ptr->val;
	c_start         = tree_ptr->start;
	tree_ptr->val   = t_val;
	tree_ptr->start = t_start;
      }
    }
    if (c_val > root_ptr->val) {
      t_start         = c_start;
      c_start         = root_ptr->start;
      root_ptr->start = t_start;
      result          = root_ptr->val;
      root_ptr->val   = c_val;
    }
    else 
      result          = c_val;
  }

  free(Tree);
  free(TempT);
  free(Source);

  return (result);
}


int select_merge_linear_i(int *A, int i, int bin_size, int bins) {
  register int j, k, b;
  int **list_ptr;
  int m;

  list_ptr = (int **)malloc(bins*sizeof(int *));
  assert_malloc(list_ptr);

  list_ptr[0] = A;
  for (j=1 ; j<bins ; j++)
    list_ptr[j] = list_ptr[j-1] + bin_size;

  for (j=0 ; j<i ; j++) {
    b = 0;
    m = *list_ptr[0];
    for (k=1 ; k<bins ; k++) 
      if (*list_ptr[k] < m) {
	b = k;
	m = *list_ptr[k];
      }
    list_ptr[b]++;
  }
  free(list_ptr);
  return (m);
}

int select_merge_i(int *A, int i, int bin_size, int bins) {
  if (i <= (bins<<3))
    return select_merge_linear_i(A, i, bin_size, bins);
  else 
    return select_merge_pway_i(A, i, bin_size, bins);
}


/************************* p-way merge ***********************/

typedef struct d_record {
  double val;
  int start;
} d_record_t;

double select_merge_pway_d(double *A, int idx, int bin_size, int bins) {
  /* Assume bins is a power of two */

  double result;
  register double v0, v1;
  double c_val,t_val;
  double **Source, **source_ptr;
  
  int i,j,k,
    c_start,t_start;
  d_record_t *root_ptr,*tree_ptr, *Tree, *TempT;

  Tree = (d_record_t *)malloc(bins*sizeof(d_record_t));
  assert_malloc(Tree); 

  TempT = (d_record_t *)malloc(bins*sizeof(d_record_t));
  assert_malloc(TempT); 

  Source = (double **)malloc(bins*sizeof(double *));
  assert_malloc(Source); 
   
   
  Source[0] = A;
  for (i=1 ; i<bins ; i++) 
    Source[i] = Source[i-1] + bin_size;
   
  source_ptr = Source;
  j          = 0;
  for (i=(bins>>1) ; i<bins ; i++) {
    v0 = **source_ptr;
    v1 = **(source_ptr + 1);
    if (v1 > v0) {
      Tree[i].val    = v1;
      Tree[i].start  = j + 1;
      TempT[i].val   = v0;
      TempT[i].start = j;
    }
    else {
      Tree[i].val    = v0;
      Tree[i].start  = j;
      TempT[i].val   = v1;
      TempT[i].start = j + 1;
    }
    source_ptr    += 2;
    j             += 2;
  }
   
  for (i=(bins/2 - 1) ; i>0 ; i--) {
    j = (i<<1);
    k = j + 1;
    if (TempT[k].val > TempT[j].val) {
      Tree[i].val    = TempT[k].val;
      Tree[i].start  = TempT[k].start;
      TempT[i].val   = TempT[j].val;
      TempT[i].start = TempT[j].start;
    }
    else {
      Tree[i].val    = TempT[j].val;
      Tree[i].start  = TempT[j].start;
      TempT[i].val   = TempT[k].val;
      TempT[i].start = TempT[k].start;
    }
  }
   
  result          = TempT[1].val; 
  c_start         = TempT[1].start;
  root_ptr = Tree + 1;

  for (j=1 ; j<idx ; j++) {
#if DEBUG
    if ((Source[c_start] - (A + (c_start*bin_size))) > bin_size) {
      fprintf(stderr,"ERROR: Out of range\n");
    }
#endif
    c_val = *(++Source[c_start]); 
    i = c_start + bins;
    while (i > 3) {
      if (c_val > Tree[i>>=1].val) {
	t_val           = c_val;
	t_start         = c_start;
	tree_ptr        = Tree+i;
	c_val           = tree_ptr->val;
	c_start         = tree_ptr->start;
	tree_ptr->val   = t_val;
	tree_ptr->start = t_start;
      }
    }
    if (c_val > root_ptr->val) {
      t_start         = c_start;
      c_start         = root_ptr->start;
      root_ptr->start = t_start;
      result          = root_ptr->val;
      root_ptr->val   = c_val;
    }
    else 
      result          = c_val;
  }

  free(Tree);
  free(TempT);
  free(Source);

  return (result);
}


double select_merge_linear_d(double *A, int i, int bin_size, int bins) {
  register int j, k, b;
  double **list_ptr;
  double m;

  list_ptr = (double **)malloc(bins*sizeof(double *));
  assert_malloc(list_ptr);

  list_ptr[0] = A;
  for (j=1 ; j<bins ; j++)
    list_ptr[j] = list_ptr[j-1] + bin_size;

  for (j=0 ; j<i ; j++) {
    b = 0;
    m = *list_ptr[0];
    for (k=1 ; k<bins ; k++) 
      if (*list_ptr[k] < m) {
	b = k;
	m = *list_ptr[k];
      }
    list_ptr[b]++;
  }
  free(list_ptr);
  return (m);
}

double select_merge_d(double *A, int i, int bin_size, int bins) {
  if (i <= (bins<<3))
    return select_merge_linear_d(A, i, bin_size, bins);
  else
    return select_merge_pway_d(A, i, bin_size, bins);
}

