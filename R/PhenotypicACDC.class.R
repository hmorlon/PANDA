setClass(
    Class = "PhenotypicACDC",
    representation = representation(
        matrixCoalescenceTimes="matrix"
    ),
    contains="PhenotypicModel"
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicACDC",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for an EB process ***\n")
            beginning <- Sys.time()
        }
        # params is the vector of parameters, but we do not know in which order. So we get back the parameter values with a, A, Gamma and the initial conditions.
        initialCondition <- object@initialCondition(params)
        aAGamma1 <- object@aAGamma(1, params)
        m0 <- initialCondition$mean[1]
        v0 <- initialCondition$var[1,1]
        sigma0 <- aAGamma1$Gamma(0)[1,1]
        r <- 2*log(sigma0 / aAGamma1$Gamma(1)[1,1])

        n = length(object@matrixCoalescenceTimes[1,])
        mean <- rep(m0, n)
        Sigma <- v0 + (sigma0**2/(2*r))*(exp(2*r*object@matrixCoalescenceTimes) - 1)

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
