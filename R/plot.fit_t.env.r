################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 01-08-2015                                               ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################

plot.fit_t.env<-function(x,...){
    
    # Rates through time function
    if(x$model=="ClimExp"){
        fun_temp<-function(x,temp,agetot,beta,sig){
            rate<-(sig*exp(beta*temp(x)))
            return(rate)
        }
    }else{
        fun_temp<-function(x,temp,agetot,beta,sig){
            rate<-sig+(beta-sig)*temp(x)
            return(rate)
        }
    }
    
    # Times steps
    t <- seq(0,x$tot_time, length.out=100)
    
    # Rates through time
    rates<-fun_temp( x=t, temp=x$env_func, agetot=x$tot_time, beta=x$b, sig=x$sig2)
    
    plot(t, rates, type='l', xlab="Times", ylab="Evolutionary rates", main="Evolutionary rate through time", ...)
}

# Allows drawing lines and superposing various results

lines.fit_t.env<-function(x,...){
    
    # Rates through time function
    if(x$model=="ClimExp"){
        fun_temp<-function(x,temp,agetot,beta,sig){
            rate<-(sig*exp(beta*temp(x)))
            return(rate)
        }
    }else{
        fun_temp<-function(x,temp,agetot,beta,sig){
            rate<-sig+(beta-sig)*temp(x)
            return(rate)
        }
    }
    
    # Times steps
    t <- seq(0,x$tot_time, length.out=100)
    
    # Rates through time
    rates<-fun_temp( x=t, temp=x$env_func, agetot=x$tot_time, beta=x$b, sig=x$sig2)
    
    lines(t, rates, type='l', ...)
}