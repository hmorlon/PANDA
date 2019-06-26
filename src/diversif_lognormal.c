// C functions Odile
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


SEXP relToAbs(SEXP lambda, SEXP parents, SEXP length){
    int len, i;
    
    
    len=INTEGER(length)[0];
    PROTECT(lambda = coerceVector(lambda,REALSXP));
    PROTECT(parents = coerceVector(parents,INTSXP));
    SEXP lambda2 = PROTECT(allocVector(REALSXP,len));
    
    // pointers pour les objets SEXP
    double *lamb2 = REAL(lambda2), *lamb = REAL(lambda);
    int *paren = INTEGER(parents);
    
    lamb2[0]=lamb[0];
    for(i = 1; i < len; i ++){
        lamb2[i]=lamb2[paren[i-1]]*lamb[i];
    }
    
    UNPROTECT(3);
    return lambda2;
    
}


SEXP absToRel(SEXP lambda, SEXP parents, SEXP length){
    int len, i;
    
    
    len=INTEGER(length)[0];
    PROTECT(lambda = coerceVector(lambda,REALSXP));
    PROTECT(parents = coerceVector(parents,INTSXP));
    SEXP lambda2 = PROTECT(allocVector(REALSXP,len));
    
    // pointers pour les objets SEXP
    double *lamb2 = REAL(lambda2), *lamb = REAL(lambda);
    int *paren = INTEGER(parents);
    
    lamb2[0]=lamb[0];
    for(i = 1; i < len; i ++){
        lamb2[i]=lamb[i]/lamb[paren[i-1]];
    }
    
    UNPROTECT(3);
    return lambda2;
     
}

SEXP loglik(SEXP lambda, SEXP lambda2, SEXP sigma, SEXP alpha, SEXP times, SEXP internalAndRoots, SEXP nNodes, SEXP root_depth){
    int n, i;
    double pi=3.141593;
    
    n=INTEGER(nNodes)[0];
    PROTECT(lambda = coerceVector(lambda,REALSXP));
    PROTECT(lambda2 = coerceVector(lambda2,REALSXP));
    PROTECT(sigma = coerceVector(sigma,REALSXP));
    PROTECT(alpha = coerceVector(alpha,REALSXP));
    PROTECT(times = coerceVector(times,REALSXP));
    PROTECT(root_depth = coerceVector(root_depth,REALSXP));
    PROTECT(internalAndRoots = coerceVector(internalAndRoots,INTSXP));
    SEXP tot = PROTECT(allocVector(REALSXP,1));
    
    // pointers to SEXP extractors
    double *lamb2 = REAL(lambda2), *lamb = REAL(lambda), *sig = REAL(sigma), *al = REAL(alpha), *t = REAL(times), *total = REAL(tot), *rd = REAL(root_depth);
    int *intAndRoot = INTEGER(internalAndRoots);
    
    total[0]=0;
    
    // contribution de la speciation
    for(i = 0; i < (n-1); i ++){
        total[0]=total[0]+log(lamb2[intAndRoot[i]]);
    }
    
    // contribution de la longueur des branches
    for(i = 0; i < (2*n); i ++){
        total[0]=total[0]-lamb2[i+1]*t[i];
    }
    
    // contribution des changements de taux (ici loi lognormale * alpha)
    total[0]=total[0]-lamb2[0]*rd[0];
    for(i = 0; i < (2*n); i ++){
        total[0]=total[0]-pow((log(lamb[i+1])-log(al[0]))/sig[0],2)/2;
    }
    total[0]=total[0]-n*log(2*pi)-2*n*log(sig[0]);
    
    UNPROTECT(8);
    
    return tot;
    
}
