/* mvmorph.h 2015-01-01 */
/* Julien Clavel        */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#include <complex.h>

static SEXP makearray (int rank, int *dim) {
  int nprotect = 0;
  int *dimp, k;
  SEXP dimx, x;
  PROTECT(dimx = NEW_INTEGER(rank)); nprotect++;
  dimp = INTEGER(dimx); 
  for (k = 0; k < rank; k++) dimp[k] = dim[k];
  PROTECT(x = allocArray(REALSXP,dimx)); nprotect++;
  UNPROTECT(nprotect);
  return x;
}


