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

loocvPhylo <- function(Y, tree, model=c("BM","OU","EB","lambda"), method=c("RidgeAlt","RidgeArch","RidgeAltapprox","LASSO","LASSOapprox"), targM=c("null","Variance","unitVariance"), REML=TRUE, up=NULL, low=NULL, tol=NULL, starting=NULL, SE=NULL, scale.height=TRUE){
    
    # Preliminary checks
    if(missing(tree)) stop("Please provide a phylogenetic tree of class \"phylo\" ")
    if(!is.ultrametric(tree)) stop("The method is not currently working with non-ultrametric trees")
    if(nrow(Y) != Ntip(tree)) stop("Length of phenotypic and phylogenetic data do not match")
    if (all(rownames(Y) %in% tree$tip.label)){
        Y <- Y[tree$tip.label,]
    }else{
        warning("rownames in Y are missing. It is assumed that they are in the same order as in the tree.")
    }
    
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
    if(ncol(Y)==1) stop("Only works with multivariate datasets")
    n <- nO <- nrow(Y)
    nC <- n-1
    if(REML==TRUE) n <- n-1
    
    # Scale the tree to unit length
    if(scale.height==TRUE) tree$edge.length <- tree$edge.length/max(branching.times(tree))
    
    p <- ncol(Y)
    Yest <- apply(Y,2,function(i) pic(i,tree))
    # Empirical covariance (contrasts) for starting the algorithm
    S <- t(Yest)%*%Yest / n
    
    # Identity matrix
    I <- diag(p)
    
    # Default penalty is Ridge "null"
    target <- matrix(0,p,p)
    
    # Identifying tips values
    tipsIndices <- which(  tree$edge[, 2] <= Ntip(tree))
    if(SE==FALSE) SE <- NULL
    
    
    # Default tolerance for the parameter search
    if(is.null(tol)){
        if(method=="RidgeArch"){
            tol = 1e-8
        }else{
            tol = 0
        }
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
        idx3 <- 3
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
        idx3 <- 3
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
        idx3 <- 3
    },
    "BM"={
        phyTrans <- function(phy, beta) return(phy)
        transform <- function(x) x
        
        
        if(method!="RidgeArch"){
            upperBound <- log(1e6)
            lowerBound <- log(tol)
        }else{
            upperBound <- 1
            lowerBound <- tol
        }
        
        idx1 <- idx2 <- 1
        idx3 <- 2
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
            
            # estimate noise
            if(!is.null(SE)){
                error = par[idx3]*par[idx3]
                tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
            }
            
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
                    Sk <- tcrossprod(Y[i,]); #Y[i,]%*%t(Y[i,]) ;
                    A <-.vec((Sridge - Sk)*Ind);
                    BC <- Iridge%*%((Semp - Sk)*Ind)%*%Iridge;
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
            
            # estimate noise
            if(!is.null(SE)){
                error = par[idx3]*par[idx3]
                tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
            }
            
            Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
            Sk <- crossprod(Yk)/n                     # Compute it in C
            if(any(!is.finite(Sk))) return(1e6)
            
            var_pic <- pruning(tr)
            var_root <- var_pic$varRoot
            var_contr <- var_pic$varNode
            
            # log-lik
            llik <- sapply(1:nC, function(x){
                Sk <- crossprod(Yk[-x,])/(n-1);
                LASSO <- glassoFast(Sk, alpha, maxIt=500);
                G <- LASSO$w;
                Gi <- LASSO$wi;
                
                Swk <- tcrossprod(Yk[x,]); #Yk[x,]%*%t(Yk[x,])
                rk <- sum(diag(Swk%*%Gi));
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
            
            # estimate noise
            if(!is.null(SE)){
                error = par[idx3]*par[idx3]
                tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
            }
            
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
            Gi <- try(chol(G), silent = TRUE)
            if(inherits(Gi, 'try-error')) return(1e6)
            
            # log-lik
            llik <- sapply(1:nC, function(x){
                # log-lik form of Hoffbeck & Landgrebe 1996
                rk <- sum(backsolve(Gi, Yk[x,], transpose = TRUE)^2)
                log(1 - beta*rk) + (rk/(1 - beta*rk))
            })
            
            ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*sum(2*log(diag(Gi))) + n*mean(llik))
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
            
            # estimate noise
            if(!is.null(SE)){
                error = par[idx3]*par[idx3]
                tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
            }
            
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
            G <- .makePenalty(Sk,alpha,target,targM)$S
            if (any(!is.finite(G))) return(1e6)
            eig <- eigen(G, symmetric=TRUE) # we can return the vectors and inverse from .makePenalty directly
            V <- eig$vectors
            d <- eig$values
            Gi <- V%*%diag(1/d)%*%t(V)
            H <- (1/(kronecker(d,d)+alpha))
            
            
            # LOO cross-validated log-likelihood
            LOObias <- function(Semp, Iridge, Sridge, Y, lambda){
                Tk <- sapply(1:nC, function(i){
                    Sk <- tcrossprod(Y[i,]); # Y[i,]%*%t(Y[i,]) ;
                    VSV <- .vec(t(V)%*%((Sridge - (Sk - lambda*target) - lambda*Iridge))%*%(V));
                    sum(VSV * H*.vec(t(V)%*%(-I%*%(Sk - Semp))%*%(V)))
                })
                
                bias <- (1/(2*n*(n-1))) * sum(Tk)
                
                return(bias)
            }
            
            klbias <- LOObias(Sk, Gi, G, Yk, alpha)
            
            ll <- -(1/n)*(-0.5 * (n*p*log(2*pi) + p*Ccov + n*sum(log(d)) + n*sum(diag(Gi%*%Sk)))) + klbias
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
            
            if(!is.null(SE)){
                error = par[idx3]*par[idx3]
                tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
            }
            
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
            llik <- try(sapply(1:nC, function(x){
                Sk <- crossprod(Yk[-x,])/(n-1)
                pen <- .makePenalty(Sk,alpha,target,targM)
                Gi <- pen$P
                detG <- sum(log(pen$ev))
                Swk <- tcrossprod(Yk[x,]) # Yk[x,]%*%t(Yk[x,])
                rk <- sum(diag(Swk%*%Gi))
                detG + rk
            }), silent = TRUE)
            
            if(inherits(llik, 'try-error')) return(1e6)
            # det of the phylo matrix
            ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*mean(llik))
            if(!is.finite(ll)) return(1e6)
            return(ll)
        }
        
    }
    )
    
    # Starting values over range of parameters for quadratic ridge and LASSO (because the tuning value is between 0->Inf)
    # computationally intensive but maybe better to ensure good starting values
    
    if(is.null(starting)){
        
        # Here we can use the forking to spread the calculus over the grid on several cores
        # we can use a randomized search to speed up the computations of very complex models
        
        if(!is.null(SE)){
            # various errors
            guess <- c(0.001,0.01,0.1,1,10)
            error_guess = sqrt(guess)
            lowerBound = c(lowerBound,0)
            upperBound = c(upperBound,Inf)
        }
        
        message("Initialization via grid search. Please wait...")
        if(method=="RidgeArch"){
            range_val <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9)
        }else{
            range_val <- log(c(1e-12, 1e-9, 1e-6, 0.01, 0.1, 1, 10, 100, 1000, 10000))
        }
        switch(model,
        "lambda"={
            mod_val <- c(0.2,0.5,0.8)
            if(!is.null(SE)){
                brute_force <- expand.grid(mod_val,range_val,error_guess)
            }else{
                brute_force <- expand.grid(mod_val,range_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[2]
        },
        "OU"={
            mod_val <- log(log(2)/(max(branching.times(tree))/c(0.1,0.5,1.5,3,8)))
            if(!is.null(SE)){
                brute_force <- expand.grid(mod_val,range_val,error_guess)
            }else{
                brute_force <- expand.grid(mod_val,range_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[2]
        },
        "EB"={
            mod_val <- -log(2)/(max(branching.times(tree))/c(0.1,0.5,1.5,3,8))
            if(!is.null(SE)){
                brute_force <- expand.grid(mod_val,range_val,error_guess)
            }else{
                brute_force <- expand.grid(mod_val,range_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[2]
        },
        "BM"={
            mod_val = NULL
            if(!is.null(SE)){
                brute_force <- expand.grid(range_val,error_guess)
            }else{
                brute_force <- expand.grid(range_val)
            }
            
            start <- brute_force[which.min(apply(brute_force,1,loocv)),]
            tuning <- start[1]
        })
        if(method=="RidgeArch"){
            cat("Best starting for the tuning: ",as.numeric(tuning))
        }else{
            cat("Best starting for the tuning: ",as.numeric(exp(tuning)))
        }
    }else{
        start <- .starting_val(starting, model, method, SE)
    }
    # Initial guesses found we start the optimization
    message("Start optimization. Please wait...")
    # Optimization of the cross-validated likelihood
    estimModel <- optim(start, fn = loocv, method="L-BFGS-B", upper=upperBound, lower=lowerBound)
    
    # Compute the scaled tree
    phy_estim <- phyTrans(tree, transform(estimModel$par[idx1]))
    if(!is.null(SE)){
        SE = (estimModel$par[idx3])*(estimModel$par[idx3])
        phy_estim$edge.length[tipsIndices] <- phy_estim$edge.length[tipsIndices] + SE
    }
    
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
    
    # End
    message("Done in ", estimModel$count[1]," iterations.")
    
    # return the results
    results <- list(loocv=estimModel$value, model.par=model.par, gamma=gamma, scaled_tree=phy_estim, model=model, method=method, p=p, n=nO, targM=targM, R=regularizedEstimates, REML=REML, Y=Y, SE=SE)
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
        D2 <- D %*% D
        sqrtM <- .sqM(D2/4 + lambda * diag(nrow(S)))
        Alt <- D/2 + sqrtM
        AltInv <- (1/lambda)*(Alt - D)
        evalues <- eigen(Alt, symmetric=TRUE, only.values = TRUE)$values
    },
    "unitVariance"={
        eig  <- eigen(S, symmetric = TRUE)
        Q <- eig$vectors
        d <- eig$values - lambda*target[1]
        evalues <- sqrt(lambda + d^2/4) + d/2
        D1 <- diag(evalues)
        D2 <- diag(1/evalues) # Inverse
        Alt <- Q %*% D1 %*% t(Q)
        AltInv <- Q %*% D2 %*% t(Q)
    },
    "null"={
        eig  <- eigen(S, symmetric = TRUE)
        Q <- eig$vectors
        d <- eig$values
        evalues <- sqrt(lambda + d^2/4) + d/2
        D1 <- diag(evalues)
        D2 <- diag(1/evalues)
        Alt <- Q %*% D1 %*% t(Q)
        AltInv <- Q %*% D2 %*% t(Q)
    }
    )
    pen <- list(S=Alt, P=AltInv, ev=evalues)
    return(pen)
}

# Matrix square root
.sqM <- function(x){
    if(!all(is.finite(x))) return(Inf)
    eig <- eigen(x, symmetric = TRUE)
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
        pen <- .makePenalty(S,tuning,Target,targM)
        P <- pen$S
        Pi <- pen$P
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

# Transform starting values
.starting_val <- function(starting, model, method, SE=NULL){
    switch(model,
    "OU"={
        if(method=="RidgeArch"){
            start <- c(log(starting[1]), starting[2])
        }else{
            start <- c(log(starting[1]), log(starting[2]))
        }
    },
    "EB"={
        if(method=="RidgeArch"){
            start <- c(starting[1], starting[2])
        }else{
            start <- c(starting[1], log(starting[2]))
        }
    },
    "BM"={
        if(method=="RidgeArch"){
            start <- starting[1]
        }else{
            start <- log(starting[1])
        }
    },
    "lambda"={
        if(method=="RidgeArch"){
            start <- c(starting[1], starting[2])
        }else{
            start <- c(starting[1], log(starting[2]))
        }
    })
    
    if(!is.null(SE) & model!="BM") start <- c(start, sqrt(starting[3]))
    if(!is.null(SE) & model=="BM") start <- c(start, sqrt(starting[2]))
    return(start)
}
