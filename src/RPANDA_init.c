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
extern SEXP fitnessFunction(SEXP X, SEXP x, SEXP r, SEXP alpha, SEXP Ncol, SEXP D, SEXP dSpace, SEXP Xs, SEXP Ys, SEXP I);
extern SEXP C_panda_covar_ou_fixed(SEXP A, SEXP alpha, SEXP sigma);
extern SEXP C_panda_covar_ou_random(SEXP A, SEXP alpha, SEXP sigma);
extern SEXP C_panda_weights (SEXP nterm, SEXP epochs, SEXP lambda, SEXP S, SEXP S1, SEXP beta, SEXP root);

static const R_CallMethodDef CallEntries[] = {
    {"relToAbs",                    (DL_FUNC) &relToAbs,                     3},
    {"relToAbsSum",                 (DL_FUNC) &relToAbsSum,                  3},
    {"absToRel",                    (DL_FUNC) &absToRel,                     3},
    {"absToRelSum",                 (DL_FUNC) &absToRelSum,                  3},
    {"loglik",                      (DL_FUNC) &loglik,                       8},
    {"fitnessFunction",             (DL_FUNC) &fitnessFunction,             10},
    {"C_panda_covar_ou_fixed",      (DL_FUNC) &C_panda_covar_ou_fixed,       3},
    {"C_panda_covar_ou_random",     (DL_FUNC) &C_panda_covar_ou_random,      3},
    {"C_panda_weights",             (DL_FUNC) &C_panda_weights,              7},
    {NULL, NULL, 0}
};

void R_init_RPANDA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
