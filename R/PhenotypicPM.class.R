setClass(
    Class = "PhenotypicPM",
    representation = representation(),
    contains="PhenotypicModel"
)

setMethod(
    f="getTipDistribution",
    signature="PhenotypicPM",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Analytical computation of tip traits distribution ***\n(Method working for models with a constant, A = aI + bU, and Gamma constant)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        # Works for models with or without the "OU" part
        if( object@aAGamma(1, params)$OU == TRUE ){
            theta <-params[3]
            psi   <-params[4]
            S     <-params[5]
            sigma <-params[6]
        }else{
            theta <-0
            psi   <-0
            S     <-params[3]
            sigma <-params[4]
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
                # update of the matrix of covariances
                Sigma <- updateBranchingMatrixSigma(Sigma, object@numbersCopy[i], object@numbersPaste[i])
            }

            # On the considered period, the model is determined by
            aAGammai <- object@aAGamma(i, params) 
            vectorU <- aAGammai$u
            delta_t <- object@period[i+1]-object@period[i]

            n <- sum(vectorU)
            d <- exp(-(psi+S)*delta_t*vectorU)      # d is a vector
            SigmaU <- as.vector(Sigma%*%vectorU)    # SigmaU is a vector. Note that SigmaU = USigma because A is symmetric.
            f <- (exp(-psi*delta_t) - exp(-(psi+S)*delta_t))/n      # f is a scalar
            Sigma <- outer(d,d)*Sigma + outer(f*d*SigmaU, vectorU) + outer(f*vectorU, SigmaU*d) + outer(f^2 * c(vectorU%*%SigmaU) * vectorU,vectorU) 
            mean <- f*( vectorU%*%mean )*vectorU + d*mean

            # The exponential or the Taylor expansion depending on the parameter values
            if( sigma != 0 ){
                if (abs((psi+S)*delta_t) > 1e-3){
                    integral1 <- (1-exp(-2*(psi+S)*delta_t))/(2*(psi+S))
                }else{
                    integral1 <- delta_t - (psi+S)*delta_t^2
                } 
                if (abs(psi*delta_t) > 1e-3) {
                    integral2 <- (1-exp(-2*psi*delta_t))/(2*psi)
                }else{
                    integral2 <-  delta_t - dpsi*delta_t^2
                }
                Sigma <- Sigma + diag(sigma^2 * integral1 *vectorU) + outer(sigma^2 * (integral2 - integral1)/n * vectorU, vectorU)
            }

            # No problem for the integral associated with a_i and the evolution of the mean
            mean <- mean + theta * (1-exp(-psi*delta_t)) * vectorU
        }

        mean <- matrix(mean,ncol=1)
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

