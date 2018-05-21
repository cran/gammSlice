/*
 *  Part of R package gammSlice
 *  Copyright (C) 2018  T.H. Pham and M.P. Wand
 *
 *  Unlimited use and distribution.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void F77_CALL(gammarand)(int *n, double * param, double *x);

void F77_SUB(gslcmc)(double *y, double *Cspl, double *CsplTy, double *Ahyp, 
                    int *ibkstt, int *ibkend, int *nMCMC, int *numObs,
                    int *ncCspl, int *idPres, int *idnum, int *numGps,
                    int * idStt, int *idEnd, int *nVec, double *xnu,
                    double *sigsq, double *xnucur, int *lennu, int *lenbki, 
                    int *lenssq,int *lnumn1, int *ifam, double *gparm, 
                    double *xncnoj, double *buSpl, double *betauO, int *idxSpl, 
                    double *uSbj, double *uSbjO, int *idxSbj, double *aux, 
                    double *ssqnu, int *msgCod);

void F77_SUB(lgunds)(int *j, double *xnucrj, double *xncnoj, int *lennu, 
                     int *lnumn1, double *Cspl, double *y, double *CsplTy, 
                     int *ncCspl, int *numObs, double *ssqnu, int *idPres, 
                     int *nVec, int *idnum, int *numGps, int *idStt, int *idEnd, 
                     int *ifam, double *buSpl, double *betauO, int *idxSpl, 
                     double *uSbj, double *uSbjO, int *idxSbj, double *ans);

void F77_CALL(urand)(int *n, double *x);

static const R_FortranMethodDef FortEntries[] = {
    {"gammarand", (DL_FUNC) &F77_SUB(gammarand), 3},
    {"gslcmc",    (DL_FUNC) &F77_SUB(gslcmc), 34},
    {"lgunds",    (DL_FUNC) &F77_SUB(lgunds), 25},
    {"urand",     (DL_FUNC) &F77_SUB(urand), 2},
    {NULL, NULL, 0}
};

void R_init_gammSlice(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
