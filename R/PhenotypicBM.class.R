setClass(
    Class = "PhenotypicBM",
    representation = representation(
        matrixCoalescenceTimes="matrix"
    ),
    contains="PhenotypicModel"
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicBM",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for a BM ***\n")
            beginning <- Sys.time()
        }
        # params is the vector of parameters, but we do not know in which order. So we get back the parameter values with a, A, Gamma and the initial conditions.
        initialCondition <- object@initialCondition(params)
        aAGamma1 <- object@aAGamma(1, params)
        m0 <- initialCondition$mean[1]
        v0 <- initialCondition$var[1,1]
        b <- aAGamma1$a(0)[1]
        sigma <- aAGamma1$Gamma(0)[1,1]

        mean <- m0 + b*diag(object@matrixCoalescenceTimes)
        Sigma <- v0 + sigma**2*object@matrixCoalescenceTimes

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        mean <- matrix(data=mean, ncol=1)
        return(list(mean = mean, Sigma = Sigma))
    }
)
