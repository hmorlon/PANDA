setClass(
    Class = "PhenotypicDD",
    representation = representation(
        matrixCoalescenceJ="matrix",
        nLivingLineages="numeric"
    ),
    contains="PhenotypicModel"
)


setMethod(
    f="getTipDistribution",
    signature="PhenotypicDD",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Computation of tip traits distribution through the analytical formula for a DD ***\n")
            beginning <- Sys.time()
        }
        # Warning : parameters must be given in the following order :
        m0 <- params[1]
        v0 <- params[2]
        r <- params[3]
        sigma0 <- params[4]
        N <- length(object@period)

        vectVariances <- rep(0, times=N)
        buffer <- 0
        for(j in 1:(N-1)){
            vectVariances[j] <- buffer
            buffer <- buffer + exp(2*r*object@nLivingLineages[j])*(object@period[j+1] - object@period[j])
        }
        vectVariances[N] <- buffer

        Sigma <- v0 + sigma0^2 * apply(object@matrixCoalescenceJ, c(1,2), function(x){ vectVariances[x] })
        mean <- rep(m0, times=length(Sigma[,1]))

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
