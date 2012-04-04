#include <R.h>
#include <R_ext/Random.h>
#include <Rmath.h>
  
void F77_CALL(gammarand)(int *n, double * param, double *x)
{
    int i;
    double a=param[0], b=param[1] ;
    GetRNGstate();
    for (i = 0; i < *n; i++) x[i] = rgamma(a,b);
    PutRNGstate();
}
