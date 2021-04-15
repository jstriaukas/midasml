#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


extern void F77_NAME(panelsglfitf)(int *nf, int *T, double *gamma, int *ngroups, int *gindex, int *nobs, int *nvars, double *x, double *y, double *pf, int *dfmax, int *pmax, int *nlam, double *flmin, double *ulam, double *eps, double *peps, int *isd, int *intr, int *maxit, int *nalam, double *b0, double *a0, double *beta, int *ibeta, int *nbeta, double *alam, int *npass, int *jerr);
extern void F77_NAME(sglfitf)(double *gamma, int *ngroups, int *gindex, int *nobs, int *nvars, double *x, double *y, double *pf, int *dfmax, int *pmax, int *nlam, double *flmin, double *ulam, double *eps, double *peps, int *isd, int *intr, int *maxit, int *nalam, double *b0, double *beta, int *ibeta, int *nbeta, double *alam, int *npass, int *jerr);

static const R_FortranMethodDef FortranEntries[] = {
  {"panelsglfitf", (DL_FUNC) &F77_NAME(panelsglfitf), 29},
  {"sglfitf",      (DL_FUNC) &F77_NAME(sglfitf),      26},
  {NULL, NULL, 0}
};

void R_init_midasml(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}