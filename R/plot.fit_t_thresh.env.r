################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 01-08-2015                                               ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################

plot.threshML<-function(x,steps=100,...){
    
    if(x$model!="Clim") stop("Plot option is only available for the climatic model ")
    # Rates through time function
    fun_temp<-function(x,temp,param){
            beta<-as.numeric(param[1])
            rate<-(exp(beta*temp(x)))
            return(rate)
        }
    
    # Times steps
    t <- seq(0,x$tot_time, length.out=steps)
    
    # Rates through time
    rates <- fun_temp( x=t, temp=x$env_func, param=x$param)
    
    plot(-t, rates, type='l', xlab="Times", ylab=bquote(paste("Evolutionary rates ", sigma)), ...)
    results<-list(time_steps=t, rates=rates)
    invisible(results)
}

# Allows drawing lines and superposing various results

lines.threshML<-function(x,steps=100,...){
    
    if(x$model!="Clim") stop("Plot option is only available for the climatic model ")
    # Rates through time function
    fun_temp<-function(x,temp,param){
            beta<-as.numeric(param[1])
            rate<-(exp(beta*temp(x)))
            return(rate)
        }
    
    # Times steps
    t <- seq(0,x$tot_time, length.out=steps)
    
    # Rates through time
    rates <- fun_temp( x=t, temp=x$env_func, param=x$param)
    
    lines(-t, rates, type='l', ...)
    results<-list(time_steps=t, rates=rates)
    invisible(results)
}
