#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void permute(void *, void *, void *, void *, void *, void *, void *, void *);
extern void permuteKendall(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
	{"permute",				(DL_FUNC) &permute,     		8},
	{"permuteKendall",    	(DL_FUNC) &permuteKendall,     	8},
    {NULL, NULL, 0}
};

void R_init_RPANDA(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
