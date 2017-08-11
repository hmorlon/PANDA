################################################################################
##                                                                            ##
##  RPANDA : Generalized Information Criterion (GIC) for high-dimensional     ##
##           Penalized Likelihood Comparative Methods                         ##
##                                                                            ##
##   Julien Clavel - 01-08-2017                                               ##
##   require: mvMORPH, glassoFast                                             ##
##                                                                            ##
################################################################################



gic_criterion <- function(Y, tree, model="BM", method=c("RidgeAlt","RidgeArch","LASSO","ML"), targM=c("null","Variance","unitVariance"), param=NULL, tuning=0, REML=TRUE, ...){
  
  # ellipsis for additional arguments
  par <- list(...)
  if(is.null(par[["fit"]])){ fit <- FALSE}else{ fit <- par$fit}
  
  # Select the method
  method <- match.arg(method)[1]
  
  # Select the method
  targM <- match.arg(targM)[1]
  
  # check for parameters
  if(is.null(param) & model!="BM") stop("please provide a parameter value for the evolutionary model!!")
  if(method=="ML" & p>=n) warning("The covariance matrix is singular, the log-likelihood (and the GIC) is unreliable!!")
  
  # Choose the model
  switch(model,
         "lambda"={
           phyTrans <- function(phy, lambda)
           {
             rootOrig <- max(branching.times(phy))
             tips <- match(c(1:Ntip(phy)), phy$edge[,2])
             phy$edge.length <- phy$edge.length * lambda
             phy$edge.length[tips] <- phy$edge.length[tips] + (rootOrig * (1-lambda))
             return(phy)
           }
           
           mod.par=1
         },
         "OU"={
           phyTrans<-function(phy,alpha)
           {
             if(alpha<=1e-8) return(phy)
             times <- branching.times(phy)
             names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             Tmax<-times[1]
             phy2<-phy
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[which(names(times) == phy$edge[i, 1])]
               t1 <- max(times) - age
               t2 <- t1+bl
               phy2$edge.length[i] <- (1/(2*alpha))*exp(-2*alpha * (Tmax-t2)) * (1 - exp(-2 * alpha * t2)) -
                 (1/(2*alpha))*exp(-2*alpha * (Tmax-t1)) * (1 - exp(-2 * alpha * t1))
             }
             phy <- phy2
             return(phy)
           }
           
           mod.par=1
         },
         "EB"={
           phyTrans <- function (phy, beta)
           {
             if(abs(beta)<=.Machine$double.eps) return(phy)
             
             times <- branching.times(phy);
             
             names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age = times[which(names(times) == phy$edge[i, 1])]
               t1 = max(times) - age
               t2 = t1+bl
               phy$edge.length[i] = (exp(beta*t2)-exp(beta*t1))/(beta)
             }
             return(phy)
           }
           
           mod.par=1
           
         },
         "BM"={
           phyTrans <- function(phy, beta) return(phy)
           mod.par = 0
         })
  n <- nrow(Y)
  p <- ncol(Y)
  nC <- n
  if(REML==TRUE) n <- n-1
  
  # construct the covariance matrix
  if(fit==FALSE){
    tr <- phyTrans(tree, param)
  }else{
    tr <- tree
  }
  
  D <- pruning(tr)$sqrtMat
  Yi <- D%*%Y
  X <- D%*%matrix(1,nrow(Y))
  beta <- pseudoinverse(X)%*%Yi
  Yk <- Yi-X%*%beta
  S <- crossprod(Yk)/n
  
  # Determinant for the phylogenetic tree
  var_pic <- pruning(tr)
  var_root <- var_pic$varRoot
  var_contr <- var_pic$varNode
  
  if(REML==TRUE){
    Ccov <- sum(log(var_contr))
  }else{
    Ccov <- sum(log(c(var_root,var_contr)))
  }
  
  # Switch between targets
  if(method=="RidgeAlt"){
    switch(targM,
           "null"={Target <- matrix(0,p,p)},
           "Variance"={Target <- diag(1/diag(S))},
           "unitVariance"={Target <- diag(1/mean(diag(S)),p)})
  }else if(method=="RidgeArch"){
    switch(targM,
           "null"={Target <- matrix(0,p,p)},
           "Variance"={Target <- diag(diag(S))},
           "unitVariance"={Target <- diag(mean(diag(S)),p)})
  }
  
  # Construct the penalty term
  switch(method,
         "RidgeAlt"={
           P <- .makePenalty(S,tuning,Target,targM)
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%diag(1/d)%*%t(V)
         },
         "RidgeArch"={
           P <- (1-tuning)*S + tuning*Target
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%diag(1/d)%*%t(V)
         },
         "LASSO"={
           LASSO <- glassoFast(S,tuning)
           P <- LASSO$w
           Pi <- LASSO$wi
         },
         "ML"={
           P <- S
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%diag(1/d)%*%t(V)
         }
  )
  
  # GIC score
  
  if(method=="RidgeArch"){
    # 1) Hessian matrix (for Ridge)
    H <- (1/(0.5*(kronecker(d,d))))
    
    # 2) First derivative of the functional
    T1 <- sapply(1:nC, function(i){
      Sk <- Yk[i,]%*%t(Yk[i,]) ;
      VSV <- .vec(t(V)%*%(0.5*(P - (1-tuning)*Sk - tuning*Target))%*%(V))
      VSV2 <- .vec(t(V)%*%(0.5*(P - (1-tuning)*Sk - tuning*Target))%*%(V)) # if we use the patterned matrix ONE
      sum(VSV * (H*VSV2))
    })
    
    df = sum(T1)/n
    sigma_df <- df
    
  }else if(method=="LASSO" | method=="ML"){
    # LASSO or ML
    Tf2 <- function(S, P) {
      I <- ifelse(P==0,0,1) ;
      t(.vec(S*I))%*%.vec(P%*%(S*I)%*%P)
    }
    
    sigma_df <- (1/(2*n))*sum(sapply(1:nC, function(i){
      Sk <- Yk[i,]%*%t(Yk[i,]) ;
      Tf2(Sk, Pi)})) - (1/2)*Tf2(S,Pi)
    
  }else if(method=="RidgeAlt"){
    # Alternative Ridge
    H <- (1/(0.5*(kronecker(d,d)+tuning)))
    
    # 2) First derivative of the functional
    T1 <- sapply(1:nC, function(i){
      Sk <- Yk[i,]%*%t(Yk[i,]) ;
      VSV <- .vec(t(V)%*%(0.5*(P - (Sk - tuning*Target) - tuning*Pi))%*%(V))
      VSV2 <- .vec(t(V)%*%(0.5*(P - Sk))%*%(V))
      sum(VSV * (H*VSV2))
    })
    
    df = sum(T1)/n
    sigma_df <- df
  }
  
  # Number of parameters for the root state:
  # The Information matrix from the Hessian
  H2 <- pseudoinverse(Pi)
  gradient <- Pi%*%t(Yk)
  J2 <- tcrossprod(gradient)/n
  beta_df <- .tr(H2%*%J2)
  
  # LogLikelihood (minus)
  DP <- as.numeric(determinant(P)$modulus)
  llik <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*DP + n*.tr(S%*%Pi))
  GIC <- 2*llik + 2*(sigma_df+beta_df+mod.par)
  
  # return the results
  results <- list(LogLikelihood=-llik, GIC=GIC, p=p, n=n, bias=sigma_df+beta_df+mod.par)
  class(results) <- "gic.rpanda"
  return(results)
}

# Extractor for fit_pl.rpanda 'class'
GIC <- function(object){
  if(!inherits(object,"fit_pl.rpanda")) stop("only works with \"fit_pl.rpanda\" class objects")
  gic_criterion(Y=object$Y, tree=object$scaled_tree, model=object$model, method=object$method, targM=object$targM, param=object$model.par, tuning=object$gamma, REML=object$REML, fit=TRUE)
}
