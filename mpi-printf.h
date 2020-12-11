#ifndef _MPI_PRINTF_H
#define _MPI_PRINTF_H

#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

void MPI_printf(char *fmt, ...);
void MPI_fprintf(FILE *fp, char *fmt, ...);

#endif
