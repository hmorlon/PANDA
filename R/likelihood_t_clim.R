################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 01-08-2015                                               ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################



likelihood_t_env<-function(phylo, data, par=NULL, model=c("ClimExp", "ClimLin")){

# Use mvMORPH for computing the log-likelihood
require(mvMORPH)
##---------------------Functions_used_internally--------------------------------##

    ## Function to scale the tree to parameters of the climatic model
    CLIMtransform<-function(phy, beta, mtot, times, funEnv, sigma=NULL, model, tips){
        
        # Not yet used (fixed at 0), depends on wether the tree have extant species
        maxdiff<-0
        res <- phy
        
        if(model=="ClimExp"){
            
        # because the curve start from the present to the past and provided values go from the past to the present
            f<-function(x){sigma*exp(beta*funEnv((mtot+maxdiff)-x))}
            
        }else if(model=="ClimLin"){
            # sigma is explicitely introduced here
            f<-function(x){sigma+(beta-sigma)*funEnv((mtot+maxdiff)-x)}
        }
        
        # Transforms the branch-lengths of the tree
        for (i in 1:length(phy$edge.length)) {
            bl <- phy$edge.length[i]
            age <- times[phy$edge[i, 1] - tips]
            res$edge.length[i] <- integrate(f,lower=age, upper=(age + bl), subdivisions=200,rel.tol = .Machine$double.eps^0.05)$value
        }
        phy<-res
        return(phy)
    }

##------------------------------------------------------------------------------##

## Parameterization

# Default model
    model<-model[1]
    
# Number of tips
    tips <- length(phylo$tip.label)
    
# Check if the climatic function is provided
    if(is.null(par[["fun"]])){
        stop("Please provide a time-function")
    }else if(!is.function(par$fun)){
        
        # env_data is a dataframe with two columns: 1) is time; 2) is datapoints
        if(is.null(par[["df"]])){
            par$df <- smooth.spline(par$fun[,1], par$fun[,2])$df
        }
        spline_result <- sm.spline(par$fun[,1],par$fun[,2], df=par$df)
        env_func <- function(t){predict(spline_result,t)}
        
        # if we don't provide a time step in par we take the time steps of the dataset?
        t<-unique(par$fun[,1])
        
        # If the ClimLin model is used we should standardize the environmental curve: i.e. sigma is the maximum rate value.
        # the user can choose by specifying it in the "par" list
        
        if(is.null(par[["scale"]])){
            # We build the interpolated smoothing spline function
            par$fun<-splinefun(t,env_func(t))
        }else{
            curve_int<-env_func(t)
            curve_scaled=scale(curve_int,min(curve_int,na.rm=T),max(curve_int, na.rm=T)-min(curve_int,na.rm=T))
            par$fun<-splinefun(t,curve_scaled)
        }
        
        # Otherwise we assume that the environnemental function is given
    }
    
# Check if the branching time is provided
    if(is.null(par[["times"]])){
        warning("The branching time for the \"phylo\" object was not provided by the user")
        par$times<-branching.times(phylo)
        
        # Set the root to zero
        par$times<-max(par$times)-par$times
    }

# Check if the root value is provided (could be used with an mcmc setting)
    if(is.null(par[["mu"]])){
        par$mu<-NULL
    }

# Check if the tree is in prunning-wise order
    if(is.null(par[["check"]])){
        par$check<-TRUE
    }

# Check if there is polytomies
    if (phylo$Nnode != tips - 1) {
        stop("You can't use this function with polytomies, transform the tree using \"multi2di\" function first")
    }
    
# Check if measurment error is provided
    if(is.null(par[["error"]])){
        is_error<-FALSE
    }else{
        is_error<-TRUE
    }

# root age
    mtot<-max(par$times)




## Transform the tree and return the log-likelihood

if(model=="ClimExp"){
    
    # Check the parameters
    if(is.null(par[["sig2"]]) | is.null(par[["beta"]]))  {
         stop("Please provide parameters values for \"sig2\" and \"beta\" ")
    }
     
    # Sigma is not provided but analytically computed instead
    phylo <- CLIMtransform(phylo, beta=par$beta, mtot=mtot, times=par$times, funEnv=par$fun, sigma=par$sig2, model=model, tips=tips)
   
    # Add measurement error
    if(is_error){
        phylo$edge.length[par$index_error]<-phylo$edge.length[par$index_error]+par$error^2 # assume the "se" are provided in the error vector
    }
   
 
    # Compute the log-likelihood
    LL<-mvLL(phylo,data,method="pic",param=list(estim=FALSE, check=par$check, mu=par$mu, sigma=1))$logl
   
   
   

}else if(model=="ClimLin"){
    
    # Check the parameters
    if(is.null(par[["sig2"]]) | is.null(par[["beta"]])) {
        stop("Please provide parameters values for \"sig2\" and \"beta\" ")
    }
    
    # Transform the tree
    phylo <- CLIMtransform(phylo, beta=par$beta, mtot=mtot, times=par$times, funEnv=par$fun, sigma=par$sig2, model=model, tips=tips)
    
    # Add measurement error
    if(is_error){
        phylo$edge.length[par$index_error]<-phylo$edge.length[par$index_error]+par$error^2 # assume the "se" are provided in the error vector
    }
    
    # Compute the log-likelihood
    LL<-mvLL(phylo,data,method="pic",param=list(estim=FALSE, check=par$check, mu=par$mu, sigma=1))$logl
}

if(is.na(LL) | is.infinite(LL)){
return(-1000000)
} # If we use Infinity some optimizer fails; e.g. L-BFGS-B
return(LL)
    
}