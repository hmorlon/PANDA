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

