#include "nas_rand.h"

#define _nas_rand_A 1220703125.00

/******************************************************************************
C
      FUNCTION RANDLC (X, A)
C
C   This routine returns a uniform pseudorandom double precision number in the
C   range (0, 1) by using the linear congruential generator
C
C   x_{k+1} = a x_k  (mod 2^46)
C
C   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
C   before repeating.  The argument A is the same as 'a' in the above formula,
C   and X is the same as x_0.  A and X must be odd double precision integers
C   in the range (1, 2^46).  The returned value RANDLC is normalized to be
C   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
C   the new seed x_1, so that subsequent calls to RANDLC using the same
C   arguments will generate a continuous sequence.
C
C   This routine should produce the same results on any computer with at least
C   48 mantissa bits in double precision floating point data.  On Cray systems,
C   double precision should be disabled.
C
C   David H. Bailey     October 26, 1990
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      SAVE KS, R23, R46, T23, T46
      DATA KS/0/
C
C   If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
C   T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
C   by merely using the ** operator, in order to insure that the results are
C   exact on all systems.  This code assumes that 0.5D0 is represented exactly.
C
******************************************************************************/
double randlc(double *X)
{
      static int	KS;
      static double	R23, R46, T23, T46;
      double		T1, T2, T3, T4;
      double		A1;
      double		A2;
      double		X1;
      double		X2;
      double		Z;
      int		i, j;

      if (KS == 0) 
      {
        R23 = 1.0;
        R46 = 1.0;
        T23 = 1.0;
        T46 = 1.0;
    
        for (i=1; i<=23; i++)
        {
          R23 = 0.50 * R23;
          T23 = 2.0 * T23;
        }
        for (i=1; i<=46; i++)
        {
          R46 = 0.50 * R46;
          T46 = 2.0 * T46;
        }
        KS = 1;
      }
/*
|   Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.
*/
      T1 = R23 * _nas_rand_A;
      j  = T1;
      A1 = j;
      A2 = _nas_rand_A - T23 * A1;
/*
|   Break X into two parts such that X = 2^23 * X1 + X2, compute
|   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
|   X = 2^23 * Z + A2 * X2  (mod 2^46).
*/
      T1 = R23 * *X;
      j  = T1;
      X1 = j;
      X2 = *X - T23 * X1;
      T1 = A1 * X2 + A2 * X1;
      
      j  = R23 * T1;
      T2 = j;
      Z = T1 - T23 * T2;
      T3 = T23 * Z + A2 * X2;
      j  = R46 * T3;
      T4 = j;
      *X = T3 - T46 * T4;
      return(R46 * *X);
} 

static double	_nas_rand_seed;
#define _NAS_BITS 19

void nas_srand() {
    _nas_rand_seed = 314159265.00;
}

int nas_rand() {
    
    register int    k;
    double  x;

    k    = (1<<_NAS_BITS) / 4;

    x  = randlc(&_nas_rand_seed);
    x += randlc(&_nas_rand_seed);
    x += randlc(&_nas_rand_seed);
    x += randlc(&_nas_rand_seed);

    return (k*x);
}

