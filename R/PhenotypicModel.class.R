setClass(
    Class = "PhenotypicModel",
    representation = representation(
        name= "character",
        period = "numeric",
        aAGamma = "function",
        numbersCopy = "numeric",
        numbersPaste = "numeric",
        initialCondition = "function",
        paramsNames = "character",
        constraints = "function",
        params0 = "numeric",
        comment = "character"
    ),
    prototype=prototype(
        name = "BMtest",
        period = c(0,1,2,3,4,5,6),
        aAGamma = function(i, params){
            functiona <- function(t){
                return(rep(0,i+1))
            }
            matrixA <- diag(0, i+1)
            functionGamma <- function(t){
                return(diag(params[1], i+1))
            }
            return(list(a=functiona, A=matrixA, Gamma=functionGamma))
        },
        numbersCopy = c(1, 1, 2, 1, 2, 5),
        numbersPaste = c(2, 3, 4, 5, 6, 7),
        initialCondition = function(params){
            return(list(mean=c(0,0), var=c(0,0)))
        },
        paramsNames = c("sigma"),
        constraints = function(params){
            return(params[1] > 0)
        },
        params0 = c(1),
        comment = "Toy model defined by defaut"
    ),
    validity=function(object){
        if( length(object@numbersCopy) != length(object@numbersPaste) ){
            stop("[PhenotypicModel : validation] The sequence of positions of branching lineages and the sequence of new positions for the traits in the newly born lineages should have the same length.")
        }
        if( length(object@numbersCopy) != length(object@period) ){
            stop("[PhenotypicModel : validation] The sequence of positions of branching lineages and the sequence of time periods should have the same length.")
        }
        if( length(object@params0) != length(object@paramsNames) ){
            stop("[PhenotypicModel : validation] There should be the same number of defaut parameters and parameter names.")
        }
        return(TRUE)
    }
)


###################################
#    Getters and setters
###################################

setMethod(
    f="[",
    signature="PhenotypicModel",
    definition=function(x,i,j,drop){
        switch( EXPR=i,
                "name"={return(x@name)},
                "period"={return(x@period)},
                "aAGamma"={return(x@aAGamma)},
                "numbersCopy"={return(x@numbersCopy)},
                "numbersPaste"={return(x@numberPaste)},
                "initialCondition"={return(x@initialCondition)},
                "paramsNames"={return(x@paramsNames)},
                "constraints"={return(x@constraints)},
                "params0"={return(x@params0)},
                "comment"={return(x@comment)},
                stop("This variable name does not exist !")
        )
    }
)

setReplaceMethod(
    f="[",
    signature="PhenotypicModel",
    definition=function(x,i,j,value){
        switch( EXPR=i,
                "name"={x@name <- value},
                "period"={x@period <- value},
                "aAGamma"={x@aAGamma <- value},
                "numbersCopy"={x@numbersCopy <- value},
                "numbersPaste"={x@numberPaste <- value},
                "initialCondition"={x@initialCondition <- value},
                "paramsNames"={x@paramsNames <- value},
                "constraints"={x@constraints <- value},
                "params0"={x@params0 <- value},
                "comment"={x@comment <- value},
                stop("This variable name does not exist !")
        )
        validObject(x)
        return(x)
    }
)

###################################
#    Affichage
###################################

setMethod(
    f="print",
    signature="PhenotypicModel",
    definition=function(x, ...){
        cat("****************************************************************\n")
        cat("*** Object of Class PhenotypicModel *** \n")
        cat("*** Name of the model : ")
        print(x@name)
        cat("*** Parameters of the model : ")
        print(x@paramsNames)
        cat("*** Description : ")
        cat(x@comment)
        cat(paste("\n*** Periods : the model is cut into ", length(x@period), " parts. \n"))
        print(x@period)
        cat("*** Lineages branching (to be copied at the end of the corresponding period) :\n")
        print(x@numbersCopy)
        cat("*** Positions of the new trait at the end of each period :\n")
        print(x@numbersPaste)
        cat("*** Initial condition :\n")
        print(x@initialCondition)
        cat("*** Vectors a_i, A_i, Gamma_i on each period i : \n")
        print(x@aAGamma)
        cat("*** Constraints on the parameters : \n")
        print(x@constraints)
        cat("*** Defaut parameter values : ")
        print(x@params0)
        cat("****************************************************************\n")
    }
)

setMethod(
    f="show",
    signature="PhenotypicModel",
    definition=function(object){
        cat("****************************************************************\n")
        cat("*** Object of Class PhenotypicModel *** \n")
        cat("*** Name of the model : ")
        print(object@name)
        cat("*** Parameters of the model : ")
        print(object@paramsNames)
        cat("*** Description : ")
        cat(object@comment)
        cat(paste("\n*** Periods : the model is cut into ", length(object@period), " parts. \n"))
        cat("For more details on the model, call : print(PhenotypicModel)\n")
        cat("****************************************************************\n")
    }
)

###################################
#    Distribution
###################################

getPositionVectSigma <- function(n,i,j){
    # finds the position of Cov(X^i,X^j) in the one-D vector of covariances
    return((n+1-i/2)*(i-1)+(j-i+1))
}

getOldPosition <- function(i, copy, paste){
    # We are in position (i,j) of the new matrix Sigma : which value of the old Sigma, with indices (k,l), do we put there ?
    if( i < paste ){
        k <- i
    }
    else if( i == paste){
        k <- copy
    }
    else if( i > paste ){
        k <- i-1
    }
    return(k)
}

updateBranchingVectSigma <- function(sigma, copy, paste){
    # copy of a branching lineage in the one-D vector of covariances 'sigma'
    n = (1/2)*(-1+sqrt(1+8*length(sigma)))
    L = (n+1)*(n+2)/2
    newsigma <- rep(0, times=L)
    for(i in 1:(n+1)){
        for(j in i:(n+1)){
            k <- getOldPosition(i, copy, paste)
            l <- getOldPosition(j, copy, paste)
            newsigma[getPositionVectSigma(n+1,i,j)] <- sigma[getPositionVectSigma(n,k,l)]
        }
    }
    return(newsigma)
}

updateBranchingMatrixSigma <- function(Sigma, copy, paste){
    # copy of a branching lineage in the matrix of covariances 'Sigma'
    n = length(Sigma[1,])
    newSigma <- diag(0, n+1)
    for(i in 1:(n+1)){
        for(j in 1:(n+1)){
            k <- getOldPosition(i, copy, paste)
            l <- getOldPosition(j, copy, paste)
            newSigma[i,j] <- Sigma[k,l]
        }
    }
    return(newSigma)
}

setGeneric(
    name="getTipDistribution2",
    def=function(object="PhenotypicModel", params="numeric", v="boolean"){standardGeneric("getTipDistribution2")}
)

setMethod(
    f="getTipDistribution2",
    signature="PhenotypicModel",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through ODE resolution (Riccati equation: dX = AX +XA +B) ***\n(Method working for any model)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        # Sur chaque periode [t_i, t_i+1[ :
        for(i in 1:(length(object@period)-1)){

            # If there is a branching event at the beginning of the period, we update the mean and covariances
            if(object@numbersPaste[i] != 0){
                # update of the vector of means
                if( object@numbersPaste[i] <= length(mean) ){
                    mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)] )
                }else{
                    mean <- c( mean, mean[object@numbersCopy[i]] )
                }
                # update of the matrix of covariances
                Sigma <- updateBranchingMatrixSigma(Sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params)
            ai <- aAGammai$a
            Ai <- aAGammai$A
            Gammai <- aAGammai$Gamma
            n = length(mean)
            # We now need to build the ODE system such that dSigma/dt = -A Sigma - Sigma A + Gamma
            derivativeSigma <- function(t,y,params){
                X = matrix(y,nrow=n)
                dX <- -Ai %*% X - X %*% Ai + Gammai(t)
                return(list(dX))
            }

            # And we build a second ODE system such that dm/dt = -Ai m + ai
            derivativemean <- function(t,y,params){
                return(list(-Ai %*% y + ai(t)))
            }

            # We update the vectors of means and covariances through their ODE system resolution
            times <- c(object@period[i], object@period[i+1])
            if((object@period[i+1]-object@period[i])> 1e-15 ){
                mean  <- ode(mean, times, derivativemean)[2, 2:(n+1)]
                sigma <- ode(as.vector(Sigma), times, derivativeSigma)[2, 2:(n*n+1)]
            }
            Sigma = matrix(sigma,nrow=n)
        }
	mean <- matrix(data=mean, ncol=1)

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)

setGeneric(
    name="getTipDistribution",
    def=function(object="PhenotypicModel", params="numeric", v="boolean"){standardGeneric("getTipDistribution")}
)

setMethod(
    f="getTipDistribution",
    signature="PhenotypicModel",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through ODE resolution ***\n(Method working for any model)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        # Sigma is a matrix, we need to transform it into a one-D vector of covariances
        n <- length(Sigma[,1])
        sigma <- c()
        for(i in 1:n){
            for(j in i:n){
                sigma <- c(sigma, Sigma[i,j])
            }
        }

        # Sur chaque periode [t_i, t_i+1[ :
        for(i in 1:(length(object@period)-1)){

            # If there is a branching event at the beginning of the period, we update the mean and covariances
            if(object@numbersPaste[i] != 0){
                # update of the vector of means
                if( object@numbersPaste[i] <= length(mean) ){
                    mean <- c( mean[1:(object@numbersPaste[i]-1)], mean[object@numbersCopy[i]], mean[object@numbersPaste[i]:length(mean)] )
                }else{
                    mean <- c( mean, mean[object@numbersCopy[i]] )
                }
                # update of the vector of covariances
                sigma <- updateBranchingVectSigma(sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params)
            ai <- aAGammai$a
            Ai <- aAGammai$A
            Gammai <- aAGammai$Gamma
            n = length(Ai[,1])
            L =  n*(n+1)/2

            # We now need to build the ODE system such that dsigma/dt = Msigma + d
            #print(paste("la longeur de sigma est de ", length(sigma), " tandis que la longueur determinee a partir de a est ", L))
            #print(paste("La longeur de mean est de ", length(mean), "tandis que la longueur de A est", n))
            #print(paste("la taille de Gamma est ", length(Gammai(0)[,1]), "par ", length(Gammai(0)[1,]) ))
            # M does not depend on time t because Ai is a constant matrix on the period
            M = diag(0, L)
            for(k in 1:n){
                for(l in k:n){
                    # On the line concerning Cov(k,l) :
                    p1 <- getPositionVectSigma(n,k,l)
                    for(m in 1:n){                        
                        p2 <- getPositionVectSigma(n, min(m,l), max(m,l))
                        M[p1, p2] <- M[p1, p2] - Ai[k,m]
                        
                        p3 <- getPositionVectSigma(n, min(m,k), max(m,k))
                        M[p1, p3] <- M[p1, p3] - Ai[l,m]
                    }
                }
            }
            # d is a function of time, because Gammai is a function of time
            d <- function(t){
                d = rep(0, times=L)
                for(k in 1:n){
                    for(l in k:n){
                        p1 <- getPositionVectSigma(n,k,l)
                        for(m in 1:n){
                            d[p1] <- d[p1] + Gammai(t)[k,m]*Gammai(t)[l,m]
                        }
                    }
                }
                return(d)
            }
            derivativesigma <- function(t,y,params){
                return(list(M %*% y + d(t)))
            }

            # And we build a second ODE system such that dm/dt = -Ai m + ai
            derivativemean <- function(t,y,params){
                return(list(-Ai %*% y + ai(t)))
            }

            # We update the vectors of means and covariances through their ODE system resolution
            times <- c(object@period[i], object@period[i+1])
            if((object@period[i+1]-object@period[i])> 1e-15 ){
                sigma <- ode(sigma, times, derivativesigma)[2, 2:(L+1)]
                mean  <- ode(mean, times, derivativemean)[2, 2:(n+1)]
            } 
        }

        # At the end, we get m and the list of covariances sigma. We thus reconstruct a variance matrix before returning it (Sigma is a matrix, sigma is a vector)
        n = length(mean)
        Sigma <- diag(0, n)
        for(k in 1:n){
            for(l in 1:n){
                p <- getPositionVectSigma(n, min(k,l), max(k,l))
                Sigma[k,l] <- sigma[p]
            }
        }
	mean <- matrix(data=mean, ncol=1)
        #names(mean) <- c()

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)

setGeneric(
    name="getDataLikelihood",
    def=function(object="PhenotypicModel", data="numeric", params="numeric", v="boolean"){standardGeneric("getDataLikelihood")}
)

setMethod(
    f="getDataLikelihood",
    signature="PhenotypicModel",
    definition=function(object, data, params, v=FALSE){
        if(v){
            cat("*** Computing -log( likelihood ) of tip trait data under a given set of parameters ***\n")
        }

        if(object@constraints(params)){
            n <- length(data)
            tipdistribution <- getTipDistribution(object, params)

            dataminusXT <- matrix(data - tipdistribution$mean, nrow=1)
            dataminusX <- matrix(data - tipdistribution$mean, ncol=1)

            ProdVectoriel = dataminusXT %*% solve( tipdistribution$Sigma ) %*% dataminusX
            deter = det( tipdistribution$Sigma )

            calcul <-  (ProdVectoriel + log( deter ) + n*log(2*pi)) /2
            #if(deter <= 0 || ProdVectoriel < 0){
            #    calcul <- -Inf
            #}else{
            #    calcul <-  (ProdVectoriel + log( deter ) + n*log(2*pi)) /2
            #}
        }else{
            calcul <- -Inf
        }
        return(as.numeric(calcul))
    }
)


###################################
#    Parameter inferences
###################################

setGeneric(
    name="fitTipData",
    def=function(object="PhenotypicModel", data="numeric", params0="numeric"){standardGeneric("fitTipData")}
)

setMethod(
    f="fitTipData",
    signature="PhenotypicModel",
    definition=function(object, data, params0=NULL){
        cat("*** Fit of tip trait data ***\n")
        cat("Finding the maximum likelihood estimator of the parameters, before returning the likelihood and the inferred parameters...\n")
        beginning <- Sys.time()

        # If params0 is not given, we use the 'params0' value contained in the model
        if(is.null(params0)){
            params0 <- object@params0
        }

        n <- length(data)
        
        # computing the mean vector and variance matrix for the model, returns -log(likelihood) (a real number)
        toBeOptimized <- function(params){
            return(getDataLikelihood(object, data, params))
        }

        # looking for the argmin of -log(likelihood) (i.e. argmax of likelihood)
        optimisation <- optim(params0, toBeOptimized)
        inferredParams <- optimisation$par
        names(inferredParams) <- object@paramsNames

        end <- Sys.time()
        cat("Computation time :", format(end-beginning), "\n")

        return(list(value = optimisation$value, inferredParams = inferredParams))
    }
)

setGeneric(
    name="modelSelection",
    def=function(object="PhenotypicModel", data="numeric"){standardGeneric("modelSelection")}
)

setMethod(
    f="modelSelection",
    signature="PhenotypicModel",
    definition=function(object, data){
        cat("*** Model selection with tip trait data ***\n")
        cat("For each model in \"object\", fits the model and returns its AIC value in a recap table...\n")

        aic <- c()
        names <- c()
        for(model in object){
            fit <- fitTipData(model, data)
            aic <- c(aic, 2*length(model@params0)+2*fit$value )
            names <- c(names, model@name)
        }
        names(aic) <- names

        return(sort(aic))
    }
)

###################################
#    Simulation
###################################

setGeneric(
    name="simulateTipData",
    def=function(object="PhenotypicModel", params="numeric", method="numeric"){standardGeneric("simulateTipData")}
)

setMethod(
    f="simulateTipData",
    signature="PhenotypicModel",
    definition=function(object, params, method=3){
        cat("*** Simulation of tip trait values ***\n")
        if( method == 1 ){
            cat("Computing first the tip distribution, and returning a simulated dataset drawn in this distribution...\n")
            tipdistribution <- getTipDistribution(object, params)
            X <- rmvnorm(1, tipdistribution$mean, tipdistribution$Sigma)

        }else if( method == 2 ){        
            cat("Simulating step-by-step the whole trajectory of a realization of the model and plotting the whole trajectory, before returning the tip data...\n")

        }else{
            cat("Simulating step-by-step the whole trajectory of a realization of the model, before returning only the tip data...\n")

            initialCondition <- object@initialCondition(params)
            X <- rnorm(length(initialCondition$mean), initialCondition$mean, initialCondition$var)
            dt <- 0.001
            sqrtdt <- sqrt(dt)
        
            for(i in 1:(length(object@period)-1)){
                # If there is a branching event at the beginning of the period
                if(object@numbersPaste[i] != 0){
                    if( object@numbersCopy[i] < length(X) ){
                        X <- c( X[1:object@numbersCopy[i]], X[object@numbersCopy[i]], X[(object@numbersCopy[i]+1):length(X)] )
                    }else{
                        X <- c( X, X[object@numbersCopy[i]] )
                    }

                }
                # The time period is sliced
                time <- seq( from = object@period[i], to = object@period[i+1], by = dt )
                for(t in time){
                    aAGammai <- object@aAGamma(i, params)  
                    X <- X + (aAGammai$a(t) - aAGammai$A %*% X)*dt + sqrtdt* aAGammai$Gamma(t) %*% rnorm(length(X), 0, 1)
                }
            }
        }
	X <- matrix(data=X, ncol=1)
        return(X)
    }
)
