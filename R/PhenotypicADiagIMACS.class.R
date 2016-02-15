setClass(
    Class = "PhenotypicADiagIMACS",
    representation = representation(),
    contains="PhenotypicModel"
)

setMethod(
    f="getTipDistribution",
    signature="PhenotypicADiagIMACS",
    definition=function(object, params, v=FALSE){
        if(v){
            cat("*** Analytical computation of tip traits distribution ***\n(Method working for models with a constant, A = aI + bU, and Gamma constant)\n")
            beginning <- Sys.time()
        }
        # Initialisation of the distribution at the beginning of the process
        initialCondition <- object@initialCondition(params)
        mean <- initialCondition$mean
        Sigma <- initialCondition$var

        theta <-params[3]
        psi   <-params[4]
        S     <-params[5]
        sigma <-params[6] 

        # Initialisation of the eigenvalues    
        vp1 <--(psi+S)
        vp2 <--psi 

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
            vectorU <- aAGammai$v
            n = length(vectorU)  
           
            delta_t <- object@period[i+1]-object@period[i]
            vp1dt <- vp1*delta_t
            vp2dt <- vp2*delta_t
            exp_vp1dt <- exp(vp1dt)
            exp_vp2dt <- exp(vp2dt)
            exp_diff_vp <-(exp_vp2dt-exp_vp1dt)/n
            SumColSigma <- colSums(Sigma) 
            one <- rep(1,n)
            M  <- outer(SumColSigma,one)
            MT <- outer(one,SumColSigma)
            Sigma <- (exp_vp1dt^2)*Sigma + (exp_vp1dt*exp_diff_vp)*(M+MT) +  matrix(exp_diff_vp^2*sum(SumColSigma),n,n)  
            
            mean <- exp_vp1dt*mean + exp_diff_vp*sum(mean)*one

            if (sigma != 0) {
                 if (abs(vp1dt) > 1e-3) {
                    v3 <- (1-exp_vp1dt^2)/(2*vp1)
                }else{
                    v3 <- -delta_t*(1+vp1dt*(1+(2/3)*vp1dt))
                } 
                if (abs(vp2dt) > 1e-3) {
                    v4 <- (1-exp_vp2dt^2)/(2*vp2)
                }else{
                    v4 <- -delta_t*(1+vp2dt*(1+(2/3)*vp2dt))
                } 
                Sigma <- Sigma - diag(rep(sigma^2*v3,n)) + matrix(sigma^2*(v3-v4)/n,n,n)  
            }

           if ((theta*psi) != 0) {
                 if (abs(vp1dt) > 1e-3) {
                    v3 <- (1-exp_vp1dt)/vp1 
                }else{
                    v3 <- -delta_t*(1+vp1dt*(1/2+(1/6)*vp1dt))
                } 
                if (abs(vp2dt) > 1e-3) {
                    v4 <- (1-exp_vp2dt)/vp2
                }else{
                    v4 <- -delta_t*(1+vp2dt*(1/2+(1/6)*vp2dt))
                } 
                mean <- mean - (theta*psi*v4)*one
            }
        }

        if(v){
            end <- Sys.time()
            cat("Computation time :", format(end-beginning), "\n")
        }

        return(list(mean = matrix(mean,ncol=1), Sigma = Sigma))
    }
)

