setClass(
    Class = "PhenotypicOU",
    representation = representation(
        matrixCoalescenceTimes="matrix"
    ),
    contains="PhenotypicModel"
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicOU",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for an OU process ***\n")
            beginning <- Sys.time()
        }
        # params is the vector of parameters, but we do not know in which order. So we get back the parameter values with a, A, Gamma and the initial conditions.
        initialCondition <- object@initialCondition(params)
        aAGamma1 <- object@aAGamma(1, params)
        m0 <- initialCondition$mean[1]
        v0 <- initialCondition$var[1,1]
        psi <- aAGamma1$A[1,1]
        theta <- aAGamma1$a(0)[1]/psi
        sigma <- aAGamma1$Gamma(0)[1,1]

        mean <- theta + (m0 - theta)*exp(-psi*diag(object@matrixCoalescenceTimes))
        n = length(mean)
        Sigma <- diag(0, n)
        for(k in 1:n){
            for(l in k:n){
                valeur <- exp(-psi*(object@matrixCoalescenceTimes[k,k] + object@matrixCoalescenceTimes[l,l] - 2*object@matrixCoalescenceTimes[k,l])) * ( sigma**2/(2*psi) + exp(-2*psi*object@matrixCoalescenceTimes[k,l])*(v0 - sigma**2/(2*psi) ) )
                Sigma[k,l] <- valeur
                Sigma[l,k] <- valeur
            }
        }

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        mean <- matrix(data=mean, ncol=1)
        rownames(mean) <- object@tipLabels
        rownames(Sigma) <- object@tipLabels
        colnames(Sigma) <- object@tipLabels
        return(list(mean = mean, Sigma = Sigma))
    }
)
