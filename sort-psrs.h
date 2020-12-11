#ifndef _SORT_PSRS_H
#define _SORT_PSRS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int  all_sort_psrs_i(int *A, int A_size, int *B);
void all_sort_psrs_check_i(int *A, int A_size);

#endif
