setClass(
    Class = "PhenotypicADiag",
    representation = representation(),
    contains="PhenotypicModel"
)

setMethod(
    f="getTipDistribution",
    signature="PhenotypicADiag",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through integrated formula ***\n(Method working for models with a constant, A diagonalizable, and Gamma constant)\n")
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
            ai <- aAGammai$a(0)
            Ai <- aAGammai$A
            Gammai <- aAGammai$Gamma(0)
            n = length(Ai[,1])

            # The mean and covariances evolve through the period with the following formula
            e <- eigen(Ai, symmetric=TRUE)
            Q <- e$vectors
            Qt <- t(Q)

            delta_t <- object@period[i+1]-object@period[i]
            R1 <- diag( exp( -delta_t * e$values ) )
            r2 <- rep(0, n)
            r3 <- rep(0, n)
            for(k in 1:n){ 
                if(abs(delta_t*e$values[k]) > 1e-3 ){
                    r2[k] <- (1 - exp( -delta_t*e$values[k] ))/e$values[k]
                    r3[k] <- (1 - exp( -delta_t*2*e$values[k] ))/(2*e$values[k])
                }else{ # DL en 0 ordre 3 
                    r2[k] <- delta_t + 1/2*e$values[k]*delta_t^2 + 1/6*e$values[k]^2*delta_t^3
                    r3[k] <- delta_t +     e$values[k]*delta_t^2 + 2/3*e$values[k]^2*delta_t^3
                }      
            }
            R2 <- diag( r2 )
            R3 <- diag( r3 )

            mean <- Q %*% R1 %*% Qt %*% mean + Q %*% R2 %*% Qt %*% ai
            Sigma <- Q %*% R1 %*% Qt %*% Sigma %*% Q %*% R1 %*% Qt + Gammai %*% t(Gammai) %*% Q %*% R3 %*% Qt            
        }
           
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