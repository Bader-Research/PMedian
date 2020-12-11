/*									tab:8
 *
 * math_func.h - math functions
 * 
 * "Copyright (c) 1994,1995,1996 The Regents of the University of Maryland.
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
 * Creation Date:	February 15, 1995
 * Filename:		math_func.h
 * History:
 */

#ifndef _MATHFUNC_H
#define _MATHFUNC_H

#include "math.h"

#define RAND_NORM_COUNT 12

_INLINE void RNG_init()
{
    srandom(MYPROC*(int)time(0) + 17*MYPROC*(int)floor(MPI_Wtime()*317.0));
}

_INLINE double rand_uniform()
/* random() is an integer in [0, 2^31) */
/* To get uniform in [0,1), divide random() by 2^31 */
{
  /*    return (double)random() / (double)(1<<31) ; */
    return (double)random() / (double)pow(2.0,31.0);
}

_INLINE double rand_normal(register double mean, register double sigma)
/* add together > 10 uniform random numbers (-1,+1),
   return mean + sigma*sum*sqrt(3.0/count).  */

{
    register int i;
    register double sum = 0;

    for (i = 0 ; i < RAND_NORM_COUNT ; i ++) 
	sum += (double)random();

    sum /= (double)(1<<30);                /* range is [0,2*count) */
    sum -= (double)RAND_NORM_COUNT;        /* range is [-count,+count) */
    
    sum *= sigma*sqrt(3.0/(double)RAND_NORM_COUNT); /* change variance */
    sum += mean;                           /* move mean */

    return(sum);
}

/***********************************************************************/
/*  Random bits                                                        */
/***********************************************************************/

#define random_initialize() srandom (1), random_count = random ()
#define random_bitstring_size 31

static long random_bitstring;
static int  random_count;

_INLINE int random_bit (void)
{
    long r;
    
    if ( random_count >= random_bitstring_size ) {
        random_count = 0;
        random_bitstring = random ();
    }
    r = (random_bitstring & 0x1);
    random_bitstring >> 1;
    random_count ++;

    return (int)r;
}

void srrandom(unsigned int);
long rrandom();

#endif



