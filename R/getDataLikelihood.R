setGeneric(
    name="getDataLikelihood",
    def=function(object="PhenotypicModel", data="numeric", error="numeric", params="numeric", v="logical"){standardGeneric("getDataLikelihood")}
)

setMethod(
    f="getDataLikelihood",
    signature="PhenotypicModel",
    definition=function(object, data, error=NULL, params, v=FALSE){
        if(v){
            message("*** Computing -log( likelihood ) of tip trait data under a given set of parameters ***\n")
            beginning <- Sys.time()
        }

        if(!is.null(rownames(data))){
            data <- data[object@tipLabels,]
        }

        if(object@constraints(params)){
            n <- length(data)
            tipdistribution <- getTipDistribution(object, params)
            V<-tipdistribution$Sigma
            if(!is.null(error)){
            	error<-error[rownames(tipdistribution$Sigma)]
            	V<- V + diag(error^2) + diag(rep(exp(params[length(params)]),n))
			}
            
			data<-data[rownames(V)]
	
  			op <- getOption("show.error.messages")
            on.exit( options(show.error.messages=op) )
  			options(show.error.messages=FALSE)
			IV=try(solve(V))
  
  			if(inherits(IV,"try-error")){
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

        if(v){
            end <- Sys.time()
            message("Computation time :", format(end-beginning), "\n")
        }
        return(as.numeric(calcul))
    }
)
