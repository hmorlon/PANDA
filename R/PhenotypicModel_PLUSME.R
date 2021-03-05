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
        tipLabels = "character",
        tipLabelsSimu = "character",
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
        tipLabels = c("A", "B", "C", "D", "E", "F", "G"),
        tipLabelsSimu = c("A", "B", "C", "D", "E", "F", "G"),
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
                "tipLabels"={return(x@tipLabels)},
                "tipLabelsSimu"={return(x@tipLabelsSimu)},
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
                "tipLabels"={x@tipLabels <- value},
                "tipLabelsSimu"={x@tipLabelsSimu <- value},
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
        cat("*** Tip labels : \n")
        print(x@tipLabels)
        cat("*** Tip labels for simulations : \n")
        print(x@tipLabelsSimu)
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

updateBranchingMatrixSigma <- function(Sigma, copy, paste){
    # copy of a branching lineage in the matrix of covariances 'Sigma'
    n = length(Sigma[1,])
    newSigma <- diag(0, n+1)
    
    newSigma[1:(paste-1),1:(paste-1)] <- Sigma[1:(paste-1),1:(paste-1)]
    newSigma[paste,1:(paste-1)] <- Sigma[copy,1:(paste-1)]
    newSigma[1:(paste-1),paste] <- Sigma[1:(paste-1),copy]
    newSigma[paste,paste] <- Sigma[copy,copy]

    if(paste < n+1){
        newSigma[(paste+1):(n+1),1:(paste-1)] <- Sigma[paste:n,1:(paste-1)]
        newSigma[(paste+1):(n+1),paste] <- Sigma[paste:n,copy]
        newSigma[1:(paste-1),(paste+1):(n+1)] <- Sigma[1:(paste-1),paste:n]
        newSigma[paste,(paste+1):(n+1)] <- Sigma[copy,paste:n]
        newSigma[(paste+1):(n+1),(paste+1):(n+1)] <- Sigma[paste:n, paste:n]
    }

    return(newSigma)
}

setGeneric(
    name="getTipDistribution",
    def=function(object="PhenotypicModel", params="numeric", v="boolean"){standardGeneric("getTipDistribution")}
)

setMethod(
    f="getTipDistribution",
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
                Sigma = matrix(y,nrow=n)
                dSigma <- -Ai %*% Sigma - t(Sigma) %*% t(Ai) + Gammai(t) %*% t(Gammai(t))
                return(list(dSigma))
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
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = mean, Sigma = Sigma))
    }
)


setGeneric(
    name="getDataLikelihood",
    def=function(object="PhenotypicModel", data="numeric", error="numeric",params="numeric", v="boolean"){standardGeneric("getDataLikelihood")}
)

setMethod(
    f="getDataLikelihood",
    signature="PhenotypicModel",
    definition=function(object, data, error, params, v=FALSE){
        if(v){
            cat("*** Computing -log( likelihood ) of tip trait data under a given set of parameters ***\n")
        }

        if(object@constraints(params)){
            n <- length(data)
            tipdistribution <- getTipDistribution(object, params)
            error<-error[rownames(tipdistribution$Sigma)]
            V<- tipdistribution$Sigma + diag(error^2) + diag(rep(exp(params[length(params)]),n))
			data<-data[rownames(V)]
	
  			op <- getOption("show.error.messages")
  			options(show.error.messages=FALSE)
			IV=try(solve(V))
  			options(show.error.messages=op)
  			if(class(IV)=="try-error"){
    			IV=pseudoinverse(V)
  				if(max(IV)==0){return(Inf)}
  			}

            dataminusXT <- matrix(data - tipdistribution$mean, nrow=1)
            dataminusX <- matrix(data - tipdistribution$mean, ncol=1)

            ProdVectoriel = dataminusXT %*% IV %*% dataminusX

            calcul <-  (ProdVectoriel + determinant(V)$modulus + n*log(2*pi)) /2
			if(is.na(calcul) | is.infinite(calcul)){calcul=-1000000}
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
    def=function(object="PhenotypicModel", data="numeric", error="numeric", params0="numeric", GLSstyle="logical"){standardGeneric("fitTipData")}
)

setMethod(
    f="fitTipData",
    signature="PhenotypicModel",
    definition=function(object, data, error, params0=NULL, GLSstyle=TRUE){
        cat("*** Fit of tip trait data ***\n")
        cat("Finding the maximum likelihood estimator of the parameters, before returning the likelihood and the inferred parameters...\n")
        beginning <- Sys.time()

        n <- length(data)

        # If params0 is not given, we use the 'params0' value contained in the model
        if(is.null(params0)){
            params0 <- object@params0
        }
        # In "GLS-style" mode, there is an analytical expression for the first parameter, namely m0
        if(GLSstyle){
            params0 <- params0[2:length(params0)]
        }

        # computing the mean vector and variance matrix for the model, returns -log(likelihood) (a real number)
        toBeOptimized <- function(params){

            if(GLSstyle){paramsPrVerif <- c(0, params)}else{paramsPrVerif <- params}
            if(object@constraints(paramsPrVerif)){

                if(GLSstyle){
                    tipdistribution <- getTipDistribution(object, c(0,params))
            		error<-error[rownames(tipdistribution$Sigma)]
            		V<- tipdistribution$Sigma + diag(error^2) + diag(rep(exp(params[length(params)]),n))
		            data<-data[rownames(V)]
		  			op <- getOption("show.error.messages")
		  			options(show.error.messages=FALSE)
					IV=try(solve(V))
		  			options(show.error.messages=op)
		  			if(class(IV)=="try-error"){
		    			IV=pseudoinverse(V) 
		  				if(max(IV)==0){return(Inf)}
		  			}
		
                    I<-matrix(rep(1,n))

                    m0 <-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%as.matrix(data)[,1]

                    dataminusXT <- matrix(data - rep(m0, times=n), nrow=1)
                    dataminusX <- matrix(data - rep(m0, times=n), ncol=1)

                    ProdVectoriel = dataminusXT %*% IV %*% dataminusX

                    calcul <-  (ProdVectoriel + determinant(V)$modulus+ n*log(2*pi)) /2
                    params <- c(m0, params)

                }else{
                    calcul <- getDataLikelihood(object, data, error=NULL, params)
                    print(calcul)
                }

            }else{
                calcul <- -Inf
            }

            return(calcul)
        }

        # looking for the argmin of -log(likelihood) (i.e. argmax of likelihood)
        optimisation <- optim(params0, toBeOptimized,control=list(maxit=2000))
        inferredParams <- optimisation$par
        # In GLS-style, we got all parameters except the first one, 'm0' that we compute through a last call to getTipDistribution
        if(GLSstyle){
            tipdistribution <- getTipDistribution(object, c(0,inferredParams))
            error<-error[rownames(tipdistribution$Sigma)]
            V<- tipdistribution$Sigma + diag(error^2) + diag(rep(exp(inferredParams[length(inferredParams)]),n))
		  	op <- getOption("show.error.messages")
		  	options(show.error.messages=FALSE)
			IV=try(solve(V))
		  	options(show.error.messages=op)
		  	if(class(IV)=="try-error"){
		    	IV=pseudoinverse(V) 
		  		if(max(IV)==0){return(Inf)}
		  	}
		    data<-data[rownames(V)]
            I<-matrix(rep(1,n))

            m0 <-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%as.matrix(data)[,1]
            
            inferredParams <- c(m0, inferredParams)
        }
        names(inferredParams) <- object@paramsNames

        end <- Sys.time()
        cat("Computation time :", format(end-beginning), "\n")

        return(list(value = optimisation$value, inferredParams = inferredParams,convergence=optimisation$convergence))
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

