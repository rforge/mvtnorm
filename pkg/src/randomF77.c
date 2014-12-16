/* $Id$
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include "mvtnorm.h"

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }
double F77_SUB(sqrtqchisqint)(int *n, double *p) { return(sqrt(qchisq(p[0], (double) n[0], 0, 0))); }

double F77_SUB(phid)(double *x){ return pnorm(*x, 0.0, 1.0, 1, 0); }
double F77_SUB(studnt)(int *nu, double *x){ return pt(*x, *nu, 1, 0); }
    