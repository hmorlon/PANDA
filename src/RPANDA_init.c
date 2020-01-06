#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP relToAbs(SEXP lambda, SEXP parents, SEXP length);
extern SEXP relToAbsSum(SEXP lambda, SEXP parents, SEXP length);
extern SEXP absToRel(SEXP lambda, SEXP parents, SEXP length);
extern SEXP absToRelSum(SEXP lambda, SEXP parents, SEXP length);
extern SEXP loglik(SEXP lambda, SEXP lambda2, SEXP sigma, SEXP alpha, SEXP times, SEXP internalAndRoots, SEXP nNodes, SEXP root_depth);

static const R_CallMethodDef CallEntries[] = {
    {"relToAbs",                    (DL_FUNC) &relToAbs,                     3},
    {"relToAbsSum",                 (DL_FUNC) &relToAbsSum,                  3},
    {"absToRel",                    (DL_FUNC) &absToRel,                     3},
    {"absToRelSum",                 (DL_FUNC) &absToRelSum,                  3},
    {"loglik",                      (DL_FUNC) &loglik,                       8},
    {NULL, NULL, 0}
};

void R_init_RPANDA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
