#The superclass
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
      stop("[PhenotypicModel : validation] \n The sequence of positions of branching lineages \n and the sequence of new positions for the traits in the newly born lineages\n should have the same length.")
    }
    if( length(object@numbersCopy) != length(object@period) ){
      stop("[PhenotypicModel : validation] \n The sequence of positions of branching lineages \n and the sequence of time periods \n should have the same length.")
    }
    if( length(object@params0) != length(object@paramsNames) ){
      stop("[PhenotypicModel : validation] \n There should be the same number of defaut parameters \n and parameter names.")
    }
    return(TRUE)
  }
)


###################################
#    Subclasses
###################################
#a new subclass has been defined for each class of models for which the "getTipDistribution" function had been optimized.

setClass(
  Class = "PhenotypicACDC",
  representation = representation(
    matrixCoalescenceTimes="matrix"
  ),
  contains="PhenotypicModel"
)

setClass(
  Class = "PhenotypicADiag",
  representation = representation(),
  contains="PhenotypicModel"
)

setClass(
  Class = "PhenotypicBM",
  representation = representation(
    matrixCoalescenceTimes="matrix"
  ),
  contains="PhenotypicModel"
)

setClass(
  Class = "PhenotypicDD",
  representation = representation(
    matrixCoalescenceJ="matrix",
    nLivingLineages="numeric"
  ),
  contains="PhenotypicModel"
)

setClass(
  Class = "PhenotypicGMM",
  representation = representation(
    n1="numeric",
    n2="numeric"
  ),
  contains="PhenotypicModel"
)

setClass(
  Class = "PhenotypicOU",
  representation = representation(
    matrixCoalescenceTimes="matrix"
  ),
  contains="PhenotypicModel"
)

setClass(
  Class = "PhenotypicPM",
  representation = representation(),
  contains="PhenotypicModel"
)


###################################
#    Basic methods
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
            "numbersPaste"={return(x@numbersPaste)},
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
            "numbersPaste"={x@numbersPaste <- value},
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
    cat(paste("\n*** Epochs : the model is cut into ", length(x@period), " parts. \n"))
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



setGeneric(
    name="fitTipData",
    def=function(object="PhenotypicModel", data="numeric", params0="numeric", GLSstyle="logical", v="logical"){standardGeneric("fitTipData")}
)

setMethod(
    f="fitTipData",
    signature="PhenotypicModel",
    definition=function(object, data, params0=NULL, GLSstyle=FALSE, v=FALSE){
        if(v){
            cat("*** Fit of tip trait data ***\n")
            cat("Finds the maximum likelihood estimators of the parameters, \nreturns the likelihood and the inferred parameters.\n")
            cat("**WARNING** : This function uses the standard R optimizer \"optim\".\nIt may not always converge well.\nPlease double check the convergence by trying\ndistinct parameter sets for the initialisation.\n")
            beginning <- Sys.time()
        }

        n <- length(data)
        if(!is.null(rownames(data))){
            data <- data[object@tipLabels,]
        }

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
                  
		            V<-tipdistribution$Sigma
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
                    calcul <- getDataLikelihood(object, data, params)
                }

            }else{
                calcul <- -Inf
            }

            return(calcul)
        }

        # looking for the argmin of -log(likelihood) (i.e. argmax of likelihood)
        optimisation <- optim(params0, toBeOptimized)
        inferredParams <- optimisation$par
        # In GLS-style, we got all parameters except the first one, 'm0' that we compute through a last call to getTipDistribution
        if(GLSstyle){
            tipdistribution <- getTipDistribution(object, c(0,inferredParams))

		    V<-tipdistribution$Sigma
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

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(value = optimisation$value, inferredParams = inferredParams, convergence = optimisation$convergence))
    }
)
