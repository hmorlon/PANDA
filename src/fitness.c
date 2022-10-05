#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP fitnessFunction(SEXP X, SEXP x, SEXP r, SEXP alpha, SEXP Ncol, SEXP D, SEXP dSpace, SEXP Xs, SEXP Ys, SEXP I){
    
    int n, i, j, d, pos;
    
    n=INTEGER(Ncol)[0];
    d=INTEGER(D)[0];
    pos=INTEGER(I)[0];
    PROTECT( X=coerceVector(X,REALSXP));
    PROTECT(x = coerceVector(x,REALSXP));
    PROTECT(r=coerceVector(r,REALSXP));
    PROTECT( alpha=coerceVector(alpha,REALSXP));
    PROTECT( Xs=coerceVector(Xs,REALSXP));
	PROTECT( Ys=coerceVector(Ys,REALSXP));
	PROTECT( dSpace=coerceVector(dSpace,REALSXP));
	//PROTECT( fitness=coerceVector(fitness,REALSXP));
    SEXP fitness = PROTECT(allocVector(REALSXP,n));
    SEXP sumTrait = PROTECT(allocVector(REALSXP,1));
    SEXP Dist1 = PROTECT(allocVector(REALSXP,1));
    SEXP Dist2 = PROTECT(allocVector(REALSXP,1));
    
    // pointers to SEXP extractors
    double *Xt1 = REAL(X), *Xt2 = REAL(x), *R = REAL(r), *al = REAL(alpha), *fit = REAL(fitness), *st = REAL(sumTrait);
    double *xs = REAL(Xs), *ys = REAL(Ys), *dS = REAL(dSpace), *dist1=REAL(Dist1), *dist2=REAL(Dist2);

    
    // calcul des fitness
    if(*al<0){
    	for(i = 0; i < (n); i ++){
    		st[0]=0;
    		for(j = 0; j < (d); j ++){
    			st[0]=st[0]-(Xt1[j*n+i]-Xt2[j])*(Xt1[j*n+i]-Xt2[j]);
    		}
    		fit[i]= 1. /( *R - 1)+1-exp( st[0] * ( *al * *al / 2 ));
        }
    }
    else{
       	for(i = 0; i < (n); i ++){
    		st[0]=0;
    		for(j = 0; j < (d); j ++){
    			st[0]=st[0]-(Xt1[j*n+i]-Xt2[j])*(Xt1[j*n+i]-Xt2[j]);
    		}
        	fit[i]=1./( *R - 1.) + exp( st[0] * ( *al * *al /2.));
        }
    }
    
    //printf("%lf",exp( -1. /(2. * *dS * *dS)));
    for(i = 0; i < (n); i ++){
        	if(fit[i]<(1. /( *R - 1))){
        	fit[i]= 1. /( *R - 1);    
        }
        *dist1 = xs[pos-1] - xs[i];
        *dist2 = ys[pos-1] - ys[i];
        if( *dist1 <0 ){
        	*dist1 = -1. * *dist1;
        }
        if( *dist1 >(xs[n-1] / 2.)){
        	*dist1= *dist1 - xs[n-1];
        }
        if( *dist2 <0 ){
        	*dist2 = -1. * *dist2;
        }
        if( *dist2 >(ys[n-1] / 2.)){
        	*dist2= *dist2 - ys[n-1];
        }
        fit[i]=fit[i] * exp( -1. * ( *dist1 * *dist1 + *dist2 * *dist2)/(2. * *dS * *dS));
        //printf("%lf",exp( -1. * ( *dist1 * *dist1 + *dist2 * *dist2)/(2. * *dS * *dS)));
        }
    
    
    UNPROTECT(11);
    
    return fitness;
    
}

