setClass(
    Class = "PhenotypicGMM",
    representation = representation(
        n1="numeric",
        n2="numeric"
    ),
    contains="PhenotypicModel"
)

setMethod(
    f="getTipDistribution",
    signature="PhenotypicGMM",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Analytical computation of tip traits distribution ***\n(Method working for the GMM model only)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        d1 <-params[3]
        d2 <-params[4]
        S <-params[5]
        sigma <-params[6]

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
            delta_t <- object@period[i]-object@period[i+1]
            n1 <- object@n1[i]
            n2 <- object@n2[i]

            exp_SDelta <- exp(S*delta_t)
            alpha <- (exp_SDelta - 1)^2 / 2
            beta <- ( 1 - exp_SDelta^2 ) / 2
            top_part <- matrix( c(rep(alpha/n1, times=n1^2), rep(beta/n2, times=n1*n2)), nrow=n1 )
            bottom_part <- matrix( c(rep(beta/n1, times=n1*n2), rep(alpha/n2, times=n2*n2)), nrow=n2 )
            exp_DeltaA <- diag(exp_SDelta, n1+n2) + rbind(top_part, bottom_part)

            mean <- exp_DeltaA %*% mean + c( rep(-S*delta_t*(d1+d2)/2 + (1- exp_SDelta^2)*(d1-d2)/4, times=n1), rep(-S*delta_t*(d1+d2)/2 + (1- exp_SDelta^2)*(d2-d1)/4,times=n2) )

            calcul_jaune <- sigma^2*(1/n1 + 1/n2)*((1 - exp_SDelta^4)/(16*S) - (1-exp_SDelta^2)/(4*S) - delta_t/4)
            calcul_orange <- sigma^2*(1/n1 + 1/n2)*((exp_SDelta^4 - 1)/(16*S) - delta_t/4)
            top_part <- matrix( c(rep(calcul_jaune, times=n1^2), rep(calcul_orange, times=n1*n2)), nrow=n1 )
            bottom_part <- matrix( c(rep(calcul_orange, times=n1*n2), rep(calcul_jaune, times=n2*n2)), nrow=n2 )
            integrale_AtA <- diag(sigma^2*(1-exp_SDelta^2)/(2*S), n1+n2) + rbind(top_part, bottom_part)

            Sigma <- exp_DeltaA %*% Sigma %*% t(exp_DeltaA) + integrale_AtA
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

