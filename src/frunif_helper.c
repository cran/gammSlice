#include <R.h>
#include <R_ext/Random.h>
  
void F77_CALL(urand)(int *n, double *x)
{
    int i;
    GetRNGstate();
    for (i = 0; i < *n; i++) x[i] = unif_rand();
    PutRNGstate();
}
