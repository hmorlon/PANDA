################################################################################
##                                                                            ##
##  RPANDA : Leave-One-Out Cross-Validation for High-dimentional              ##
##           Penalized Likelihood Comparative Methods                         ##
##                                                                            ##
##   Julien Clavel - 01-08-2017                                               ##
##   require: mvMORPH, glassoFast                                             ##
##                                                                            ##
################################################################################


# Y        = traits matrix (columns: variables, rows: species)
# tree     = phylogenetic tree (an object of class 'phylo')
# model    = either "BM", "OU", "EB", or "lambda"; the model of traits evolution
# method   = "RidgeArch": Archetype (linear) Ridge penalty, "RidgeAlt": Quadratic Ridge penalty, "LASSO": Least Absolute Selection and Shrinkage Operator. "RidgeAltapprox" and "LASSOapprox" are fast approximations of the LOOCV for the Ridge quadratic and LASSO penalties
# targM    = "null", "Variance" for a diagonal unequal variance target, "unitVariance" for an equal diagonal target. Only works with "RidgeArch","RidgeAlt", and "RidgeAltapprox" methods.
# REML     = TRUE (default) or FALSE. The likelihood method used to estimate the parameters. REML must be preferred with small sample size in order to obtain unbiased estimates.
# up       = upper bound for model parameter search
# low      = lower bound for model parameter search
# tol      = lower bound tolerance for the regularizer (tuning parameter)
# starting = starting values for the parameter search. Must be a vector likes: c(model, regularization)


require(mvMORPH)    # >= 1.0.9
require(glassoFast) # https://github.com/JClavel/glassoFast

loocvPhylo <- function(Y, tree, model=c("BM","OU","EB","lambda"), method=c("RidgeAlt","RidgeArch","RidgeAltapprox","LASSO","LASSOapprox"), targM=c("null","Variance","unitVariance"), REML=TRUE, up=NULL, low=NULL, tol=1e-6, starting=NULL){
  
  # Checks
  if(missing(tree)) stop("Please provide a phylogenetic tree of class \"phylo\" ")
  
  # Select the model
  model <- match.arg(model)[1]
  
  # Select the method
  method <- match.arg(method)[1]
  
  # Select the method
  targM <- match.arg(targM)[1]
  
  # Bounds for models
  if(is.null(up)){
    switch(model,
           "EB"={up <- 0},
           "OU"={up <- 10},
           "lambda"={up <- 1.1})
  }
  
  if(is.null(low)){
    switch(model,
           "EB"={low <- -10},
           "OU"={low <- 0},
           "lambda"={low <- 1e-5})
  }
  
  # Parameters
  n <- nO <- nrow(Y)
  nC <- n-1
  if(REML==TRUE) n <- n-1
  
  p <- ncol(Y)
  Yest <- apply(Y,2,function(i) pic(i,tree))
  # Empirical covariance (contrasts) for starting the algorithm
  S <- t(Yest)%*%Yest / n
  
  # Identity matrix
  I <- diag(p)
  
  # Default penalty is Ridge "null"
  target <- matrix(0,p,p)
  
  # vague starting values for the regularization
  if(method=="RidgeArch"){
    tuning = abs(1 - (1/log(p)))
  }else{
    tuning = log(n)+log(p)
  }
  
  ## ---- Transform phylo
  switch(model,
         "lambda"={
           phyTrans <- function(phy, lambda) {
             if(lambda==1) return(phy)
             rootOrig <- max(branching.times(phy))
             tips <- match(c(1:Ntip(phy)), phy$edge[,2])
             phy$edge.length <- phy$edge.length * lambda
             phy$edge.length[tips] <- phy$edge.length[tips] + (rootOrig * (1-lambda))
             return(phy)
           }
           if(is.null(starting)){
             start <- c(0.5, tuning)
             }else{
             start <- starting
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,log(1e6))
             lowerBound <- c(low,log(tol))
           }else{
             upperBound <- c(up,1)
             lowerBound <- c(low,tol)
           }
           
           idx1 <- 1
           idx2 <- 2
         },
         "OU"={
           phyTrans<-function(phy,alpha){
             if(alpha<=.Machine$double.eps) return(phy) # reduce to BM
             
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
           if(is.null(starting)){
             start <- c(log(log(2)/max(branching.times(tree))*5), tuning)
           }else{
             start <- starting
           }
           
           transform <- function(x) exp(x)
           
           if(method!="RidgeArch"){
             upperBound <- c(log(up),log(1e6))
             lowerBound <- c(log(low),log(tol))
           }else{
             upperBound <- c(log(up),1)
             lowerBound <- c(log(low),tol)
           }
           
           idx1 <- 1
           idx2 <- 2
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
           if(is.null(starting)){
             start <- c(-log(2)/(max(branching.times(tree))*5), tuning)
           }else{
             start <- starting
           }
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,log(1e6))
             lowerBound <- c(low,log(tol))
           }else{
             upperBound <- c(up,1)
             lowerBound <- c(low,tol)
           }
           
           idx1 <- 1
           idx2 <- 2
         },
         "BM"={
           phyTrans <- function(phy, beta) return(phy)
           transform <- function(x) x
           
           if(is.null(starting)){
             start <- tuning
           }else{
             start <- starting
           }
           
           if(method!="RidgeArch"){
             upperBound <- log(1e6)
             lowerBound <- log(tol)
           }else{
             upperBound <- 1
             lowerBound <- tol
           }
           
           idx1 <- idx2 <- 1
         })
  ## ------ Leave-One-Out Cross-validation / Penalty
  
  switch(method,
         "LASSOapprox"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             #  Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             # LASSO penalty
             LASSO <- glassoFast(Sk, alpha, maxIt=5000)
             G <- LASSO$w
             Gi <- LASSO$wi
             
             # LOO cross-validated log-likelihood
             LOObias <- function(Semp, Iridge, Sridge, Y, lambda){
               # Indices matrice
               Ind <- ifelse(Iridge==0,0,1);
               
               Tk <- sapply(1:nC, function(i){
                 Sk <- Y[i,]%*%t(Y[i,]) ;
                 A <-.vec((Sridge - Sk)*Ind)
                 BC <- Iridge%*%((Semp - Sk)*Ind)%*%Iridge
                 sum(A*BC)
               })
               
               bias <- (1/(2*n*(n-1))) * sum(Tk)
               
               return(bias)
             }
             
             klbias <- LOObias(Sk, Gi, G, Yk, alpha)
             
             ll <- -(1/n)*(-0.5 * (n*p*log(2*pi) + p*Ccov + n*determinant(G)$modulus + n*sum(diag(Gi%*%Sk)))) + klbias
             if (!is.finite(ll)) return(1e6)
             return(ll)
           }
         },
         "LASSO"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           require(glassoFast)
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             # log-lik
             llik <- sapply(1:nC, function(x){
               Sk <- crossprod(Yk[-x,])/(n-1)
               LASSO <- glassoFast(Sk, alpha, maxIt=500)
               G <- LASSO$w
               Gi <- LASSO$wi
               
               Swk <- Yk[x,]%*%t(Yk[x,])
               rk <- sum(diag(Swk%*%Gi))
               determinant(G)$modulus + rk
             })
             
             # det of the phylo matrix
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*mean(llik))
             if (!is.finite(ll)) return(1e6)
             return(ll)
           }
         },
         "RidgeArch"={
           # transform for the tuning parameter
           transAlpha <- function(x) (x)
           
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)
             Yk <-  apply(Y,2,function(i) pic(i,tr))
             Sk <- crossprod(Yk)/n
             if(any(!is.finite(Sk))) return(1e6)
             
             # Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             if(targM=="Variance"){
               target <- diag(diag(Sk))
             }else if(targM=="unitVariance"){
               target <- I*mean(diag(Sk))
             }
             
             # Hoffbeck & Landgrebe 1996 parameterization
             beta <- (1 - alpha)/(n - 1)
             G <- n*beta*Sk + alpha * target
             Gi <- solve(G)
             
             # log-lik
             llik <- sapply(1:nC, function(x){
               # log-lik form of Hoffbeck & Landgrebe 1996
               rk <- t(Yk[x,])%*%Gi%*%Yk[x,]
               log(1 - beta*rk) + (rk/(1 - beta*rk))
             })
             
             ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*determinant(G)$modulus + n*mean(llik))
             if (!is.finite(ll)) return(1e6)
             return(ll)
           }
           
         },
         "RidgeAltapprox"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             # Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             # Choose an updated target matrix if it's not the usual Ridge
             if(targM=="Variance"){
               target <- diag(1/diag(Sk))
             }else if(targM=="unitVariance"){
               target <- I*(1/mean(diag(Sk)))
             }
             
             # Ridge penalty
             G <- .makePenalty(Sk,alpha,target,targM)
             if (any(!is.finite(G))) return(1e6)
             eig <- eigen(G)
             V <- eig$vectors
             d <- eig$values
             Gi <- V%*%diag(1/d)%*%t(V)
             
             # LOO cross-validated log-likelihood
             LOObias <- function(Semp, Iridge, Sridge, Y, lambda){
               # need to be computed only once
               H <- (1/(kronecker(d,d)+lambda))
               
               Tk <- sapply(1:nC, function(i){
                 Sk <- Y[i,]%*%t(Y[i,]) ;
                 VSV <- .vec(t(V)%*%((Sridge - (Sk - lambda*target) - lambda*Iridge))%*%(V))
                 sum(VSV * H*.vec(t(V)%*%(-I%*%(Sk - Semp))%*%(V)))
               })
               
               bias <- (1/(2*n*(n-1))) * sum(Tk)
               
               return(bias)
             }
             
             klbias <- LOObias(Sk, Gi, G, Yk, alpha)
             
             ll <- -(1/n)*(-0.5 * (n*p*log(2*pi) + p*Ccov + n*determinant(G)$modulus + n*sum(diag(Gi%*%Sk)))) + klbias
             if(!is.finite(ll)) return(1e6)
             return(ll)
           }
         },
         "RidgeAlt"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             # Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr)) 
             }else{
               Ccov <- sum(log(c(var_root,var_contr))) 
             }
             
             # Choose an updated target matrix if it's not the usual Ridge
             if(targM=="Variance"){
               target <- diag(1/diag(Sk))
             }else if(targM=="unitVariance"){
               target <- I*(1/mean(diag(Sk)))
             }
             
             # log-lik
             llik <- sapply(1:nC, function(x){
               #Sk = (n/(n-1)) * Sk - (1/(n-1))*(Y[x,]%*%t(Y[x,]))
               Sk <- crossprod(Yk[-x,])/(n-1)
               G <- .makePenalty(Sk,alpha,target,targM)
               Gi <- solve(G)
               
               Swk <- Yk[x,]%*%t(Yk[x,])
               rk <- sum(diag(Swk%*%Gi))
               determinant(G)$modulus + rk
             })
             
             # det of the phylo matrix
             ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*mean(llik))
             if(!is.finite(ll)) return(1e6)
             return(ll)
           }
           
         }     
  )
  
  # Optimization of the cross-validated likelihood
  estimModel <- optim(start, fn = loocv, method="L-BFGS-B", upper=upperBound, lower=lowerBound)
  
  # Compute the scaled tree
  phy_estim <- phyTrans(tree, transform(estimModel$par[idx1]))
  
  # Estimated value for the model parameter
  if(model=="BM"){
      model.par <- 0
  }else{
      model.par <- transform(estimModel$par[idx1])
  }
  
  # Estimated value for the regularization parameter
  gamma <- transAlpha(estimModel$par[idx2])
  # Compute R
  Ytransform <-  apply(Y,2,function(i) pic(i,phy_estim))
  Snew <- crossprod(Ytransform)/n
  matMeth <- method
  if(method == "LASSOapprox"){
      matMeth <- "LASSO"
  }else if(method == "RidgeAltapprox"){
      matMeth <- "RidgeAlt"
  }
  regularizedEstimates <- .covPenalized(S=Snew, method=matMeth, targM=targM, tuning=gamma)
  
  # return the results
  results <- list(loocv=estimModel$value, model.par=model.par, gamma=gamma, scaled_tree=phy_estim, model=model, method=method, p=p, n=nO, targM=targM, R=regularizedEstimates, REML=REML, Y=Y)
  class(results) <- "fit_pl.rpanda"
  return(results)
}

## -------------- Miscellaneous functions

# Alternative penalty of van Wieringen & Peeters 2016 - Computational Statistics and Data Analysis
# see also Witten & Tibshirani 2009
.makePenalty <- function(S,lambda,target,targM){
  
  switch(targM,
        "Variance"={
          D <- (S - lambda * target)
          Alt <- D/2 + .sqM((D %*% D)/4 + lambda * diag(nrow(S)))
         },
        "unitVariance"={
          eig  <- eigen(S, symmetric = TRUE)
          Q <- eig$vectors
          d <- eig$values - lambda*target[1]
          D <- diag(sqrt(lambda + d^2/4) + d/2)
          Alt <- Q %*% D %*% t(Q)
        },
        "null"={
          eig  <- eigen(S, symmetric = TRUE)
          Q <- eig$vectors
          d <- eig$values
          D <- diag(sqrt(lambda + d^2/4) + d/2)
          Alt <- Q %*% D %*% t(Q)
        }
  )
  
  return(Alt)
}

# Matrix square root
.sqM <- function(x){
  if(!all(is.finite(x))) return(Inf)
  eig <- eigen(x)
  sqrtM <- eig$vectors %*% diag(sqrt(eig$values)) %*% solve(eig$vectors)
  return(sqrtM)
}

# Compute the Regularized covariance and it's inverse
.covPenalized <- function(S, method, targM="null", tuning=0){
  
  # dim of S
  p = ncol(S)
  
  # init the target matrix
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
  estimate <- list(R=P, Rinv=Pi)
  return(estimate)
}


# Vec operator
.vec <- function(x) as.numeric(x)

# Trace operator
.tr <- function(x) sum(diag(x))

