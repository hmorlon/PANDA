#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP absToRel(SEXP, SEXP, SEXP);
extern SEXP loglik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP relToAbs(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"absToRel", (DL_FUNC) &absToRel, 3},
  {"loglik",   (DL_FUNC) &loglik,   8},
  {"relToAbs", (DL_FUNC) &relToAbs, 3},
  {NULL, NULL, 0}
};

void R_init_RPANDA(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
