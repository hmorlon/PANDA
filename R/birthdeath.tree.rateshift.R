library(geiger)
library(ape)
library(phytools)

birthdeath.tree.rateshift <- function (theta,lamb_par, mu_par,
                                       f.lamb=function(x,y){y}, f.mu=function(x,y){y}, lamb_shift=0,relative_death=0,
                                       new_lamb_law="uniform",new_mu_law="uniform",alpha=0, sigma=0.1,lamb_max=1,lamb_min=0,mu_min=0,mu_max=0, time.stop = 0, 
                                       sigma_mu=0,taxa.stop = Inf, return.all.extinct=TRUE, prune.extinct=TRUE, condition="time",nShiftMax=Inf)

# theta probability to have a shift at a speciation event (independant for the two new lineages)
# lamb_par parameters of the birth rate function, i.e. initial speciation rate if you keep f.lamb the default
# mu_par parameters of the birth rate function, i.e. initial extinction rate if you keep f.mu the default
# f.lamb , f.mu could theoretically be used to have rates variing through time but I have to modify the function first, so for now it must stay the default
# condition="time" the process stop after a certain time or number of taxa (if "taxa") is reached (or before if it goes extinct and return.all.extinct=TRUE)
# time.stop, taxa.stop time or number of tips before the process is stoped
  
# new_lamb_law the probability law from which new speciation rates are drawn. Can be "uniform", "normal", "lognormal", "normal*t" (normal with standard 
# deviation proportional to the length of the branch), "lognormal*shift", "normal+shift", "normal*shift" (lognormal or 
# normal with mode the previous speciation rate * (or +) lamb_shift)
  
# new_mu_law the probability law from which new death rates are drawn. Can be "uniform", "normal", "lognormal","normal*t" (normal with standard 
# deviation proportional to the length of the branch), "diversify" (constant diversification rate) or "turnover" (constant turnover rate)
# Yet I never used it with varying extinction rate, so I am not sure there would be no problem. To have constant extinction rate 
# set new_mu_law="uniform", mu_min=mu_max=mu_par

# relative_death For new_mu_law="diversify" or "turnover", diversification or turnover rate. Not used in other cases
# sigma=0.1 standard deviation for new_lamb_law = "normal", "lognormal", "normal*t", "lognormal*shift", "normal+shift", "normal*shift"
# sigma_mu standard deviation for new_mu_law = "normal", "lognormal", "normal*t"

# lamb_max,lamb_min,mu_min,mu_max, limits of the uniform laws
# return.all.extinct if FALSE the fuction is runned until we get a tree that don't become extinct before the stopping condition is reached 
# prune.extinct are extinct taxa removed from the phylogeny ?

# nShiftMax : if this number of shifts in the phylogeny is reached, theta is set to 0

# The function returns a list list with:
# list$tree the resulting tree
# list$rates the rate category for each branch
# list$lamb and list$mu the speciation and extinctio rate of each rate category (so that speciation rate of each branch is obtained 
# with list$lamb$par[list$rates])
# list$times and list$nblineages speciation times and corresponding number of lineages 

{
  lamb=list(fun=f.lamb,par=lamb_par)
  mu=list(fun=f.mu,par=mu_par)
  rates=c()
  nShift=0
  
  while (1) {
    
    nblineages<-c(1)
    times<-c(0)
    b<-lamb$fun(0,lamb$par)
    d<-mu$fun(0,mu$par)
    dt <- rexp(1,(b + d))
    #print(dt)
    t<-dt
    
    if ((t >= time.stop| 1>taxa.stop) & condition=="time") {
      t <- time.stop
    	alive<-1
    	rates<-c(1)
    	times<-c(times,t)
    	nblineages<-c(nblineages,1)
    	break
    }
    
    if (taxa.stop==1 & condition=="taxa") {
      alive<-1
      rates<-c(1)
      break
    }
            	
    r <- runif(1)
    
    if (r>b/(b + d)){
      #print("die")
      times<-c(times,dt)
      nblineages<-c(nblineages,0)
      alive<-rep(FALSE,1)
    }else{
      u=runif(2)
      for(i in 1:2){
        if(u[i]<theta & nShift<nShiftMax){
          nShift=nShift+1
          rates=c(rates,length(lamb$par)+1)
          lamb$fun=c(lamb$fun,function(x,y){exp(-alpha*x)*y})
          if(new_lamb_law=="uniform"){
            lamb$par=c(lamb$par,runif(1,min=lamb_min,max=lamb_max))
          }else if (new_lamb_law=="normal"){
              new_lambda=rnorm(1,mean=lamb$par[[1]],sd=sigma)
              while(new_lambda<0){
                new_lambda=rnorm(1,mean=lamb$par[[1]],sd=sigma)
              }
              lamb$par=c(lamb$par,new_lambda)
          }else if (new_lamb_law=="lognormal"){
              lamb$par=c(lamb$par,rlnorm(1, meanlog = log(lamb$par[[1]]),sdlog = sigma))
          }else if (new_lamb_law=="lognormal*shift"){
            lamb$par=c(lamb$par,rlnorm(1, meanlog = log(lamb$par[[1]]*lamb_shift),sdlog = sigma))
          }else if (new_lamb_law=="normal*t"){
            new_lambda=rnorm(1,mean=lamb$par[[i]],sd=sigma*t)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb$par[[i]],sd=sigma*t)
            }
            lamb$par=c(lamb$par,new_lambda)
          }else if (new_lamb_law=="normal+shift"){
            new_lambda=rnorm(1,mean=lamb$par[[1]]+lamb_shift,sd=sigma)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb$par[[1]]+lamb_shift,sd=sigma)
            }
            lamb$par=c(lamb$par,new_lambda)
          }else if (new_lamb_law=="normal*shift"){
            new_lambda=rnorm(1,mean=lamb$par[[1]]*lamb_shift,sd=sigma)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb$par[[1]]*lamb_shift,sd=sigma)
            }
            lamb$par=c(lamb$par,new_lambda)
          }
          mu$fun=c(mu$fun,function(x,y){y})
          if(new_mu_law=="uniform"){
            mu$par=c(mu$par,runif(1,min=mu_min,max=mu_max))
          }else if (new_mu_law=="normal"){
            new_mu=rnorm(1,mean=mu$par[[1]],sd=sigma)
            while(new_mu<0){
              new_mu=rnorm(1,mean=mu$par[[1]],sd=sigma)
            }
            mu$par=c(mu$par,new_mu)
          }else if (new_mu_law=="diversify"){
            new_mu=lamb$par[length(lamb$par)]-relative_death
            mu$par=c(mu$par,max(0,new_mu))
          }else if (new_mu_law=="turnover"){
            new_mu=lamb$par[length(lamb$par)]*relative_death
            mu$par=c(mu$par,new_mu)
          }else if (new_mu_law=="lognormal"){
            mu$par=c(mu$par,rlnorm(1, meanlog = log(mu$par[[1]]),sdlog = sigma_mu))
          }else if (new_mu_law=="normal*t"){
            new_mu=rnorm(1,mean=mu$par[[i]],sd=sigma_mu*t)
            while(new_mu<0){
              new_mu=rnorm(1,mean=mu$par[[i]],sd=sigma_mu*t)
            }
            mu$par=c(mu$par,new_mu)
          }
        }else{
          rates=c(rates,1)
        }
      }

    	edge <- rbind(c(1, 2), c(1, 3))
    	edge.length <- rep(NA, 2)
    	stem.depth <- rep(t, 2)
    	alive <- rep(TRUE, 2)
    	times<-c(times,dt)
    	nblineages<-c(nblineages,sum(alive))
    	next.node <- 4}
    	
    	repeat {
    	  if (sum(alive) == 0)	 break
    	  if (sum(alive)==taxa.stop){#} & condition=="taxa") {
    	    if(length(lamb$par)<=1){
    	      b=max(lamb$fun(t,lamb$par),.Machine$double.eps)
    	      d=mu$fun(t,mu$par)
    	      totalrate=sum(alive)*(b+d)
    	    }else{
    	      b<-sapply(1:length(lamb$par),function(x){max(lamb$fun[[x]](t,lamb$par[[x]]),.Machine$double.eps)})
    	      d<-sapply(1:length(mu$par),function(x){mu$fun[[x]](t,mu$par[[x]])})
    	      totalrate=sum(b[rates][alive])+sum(d[rates][alive])}
    	    
    	    dt <- rexp(1, totalrate)
    	    t <- t + dt
    	    
    	    break
    	  }else{
    		
    	  if(length(lamb$par)<=1){
    	    b=lamb$fun(t,lamb$par)
    	    d=mu$fun(t,mu$par)
    	    totalrate=sum(alive)*(b+d)
    	 }else{
    	   b<-sapply(1:length(lamb$par),function(x){max(lamb$fun[[x]](t,lamb$par[[x]]),.Machine$double.eps)})
    	   d<-sapply(1:length(mu$par),function(x){mu$fun[[x]](t,mu$par[[x]])})
         totalrate=sum(b[rates][alive])+sum(d[rates][alive])}
    		
#     	 print(totalrate)
#     	 if(totalrate==0) print(b)
        dt <- rexp(1, totalrate)
        # while(is.nan(dt)){print(b); print(d);print(totalrate);dt=rexp(1, totalrate)}
    		t <- t + dt
    		
    		if (t >= time.stop & condition=="time") {
    		  t <- time.stop
    			times<-c(times,t)
    			nblineages<-c(nblineages,sum(alive))
    			break
    		}
    		
    		r <- runif(1)
    		s=0
    		continue=T
    		i=1
    		while(i<=length(b) & continue){
    		  s=s+b[i]*sum((rates==i)[alive])/totalrate
    	
    		  if(r<=s){
    		    continue=F
    		    #print("speciation")
    		    if(length(which(rates==i & alive))>1){
    		      random_lineage = sample(which(rates==i & alive),1)
    		    }else{
    		      random_lineage = which(rates==i & alive)
    		    }
    		    parent <- edge[random_lineage, 2]
    		    alive[random_lineage] <- FALSE
    		    edge <- rbind(edge, c(parent, next.node), c(parent,next.node + 1))
    		    next.node <- next.node + 2
    		    alive <- c(alive, TRUE, TRUE)
    		    stem.depth <- c(stem.depth, t, t)
    		    #print(stem.depth)
    		    x <- which(edge[, 2] == parent)
    		    edge.length[x] <- t - stem.depth[x]
    		    edge.length <- c(edge.length, NA, NA)
    		    #print(edge.length)
    		    times<-c(times,t)
    		    nblineages<-c(nblineages,sum(alive))
    		    u=runif(2)
    		    for(j in 1:2){
    		      if(u[j]<theta & nShift<nShiftMax){
    		        nShift=nShift+1
    		        rates=c(rates,length(lamb$par)+1)
    		        lamb$fun=c(lamb$fun,function(x,y){exp(-alpha*x)*y})
    		        if(new_lamb_law=="uniform"){
    		          lamb$par=c(lamb$par,runif(1,min=lamb_min,max=lamb_max))
    		        }else if (new_lamb_law=="normal"){
    		          new_lambda=rnorm(1,mean=lamb$par[[i]],sd=sigma)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb$par[[i]],sd=sigma)
    		          }
    		          lamb$par=c(lamb$par,new_lambda)
    		        }else if (new_lamb_law=="lognormal"){
    		          lamb$par=c(lamb$par,rlnorm(1, meanlog = log(lamb$par[[i]]),sdlog = sigma))
    		        }else if (new_lamb_law=="lognormal*shift"){
    		          #mean lamb$par[[i]] iff lamb_shift==exp(-sigma^2/2)
    		          lamb$par=c(lamb$par,rlnorm(1, meanlog = log(lamb$par[[i]]*lamb_shift),sdlog = sigma))
    		        }else if (new_lamb_law=="normal*t"){
    		          new_lambda=rnorm(1,mean=lamb$par[[i]],sd=sigma*edge.length[x])
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb$par[[i]],sd=sigma*edge.length[x])
    		          }
    		          lamb$par=c(lamb$par,new_lambda)
    		        }else if (new_lamb_law=="normal+shift"){
    		          new_lambda=rnorm(1,mean=lamb$par[[i]]+lamb_shift,sd=sigma)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb$par[[i]]+lamb_shift,sd=sigma)
    		          }
    		          lamb$par=c(lamb$par,new_lambda)
    		        }else if (new_lamb_law=="normal*shift"){
    		          new_lambda=rnorm(1,mean=lamb$par[[i]]*lamb_shift,sd=sigma)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb$par[[i]]*lamb_shift,sd=sigma)
    		          }
    		          lamb$par=c(lamb$par,new_lambda)
    		        }
    		        mu$fun=c(mu$fun,function(x,y){y})
    		        if(new_mu_law=="uniform"){
    		          mu$par=c(mu$par,runif(1,min=mu_min,max=mu_max))
    		        }else if (new_mu_law=="normal"){
    		          new_mu=rnorm(1,mean=mu$par[[i]],sd=sigma_mu)
    		          while(new_mu<0){
    		            new_mu=rnorm(1,mean=mu$par[[i]],sd=sigma_mu)
    		          }
    		          mu$par=c(mu$par,new_mu)
    		        }else if (new_mu_law=="turnover"){
    		          new_mu=lamb$par[length(lamb$par)]*relative_death
    		          mu$par=c(mu$par,new_mu)
    		        }else if (new_mu_law=="diversify"){
    		          new_mu=lamb$par[length(lamb$par)]+relative_death
    		          mu$par=c(mu$par,new_mu)
    		        }else if (new_mu_law=="lognormal"){
    		          mu$par=c(mu$par,rlnorm(1, meanlog = log(mu$par[[i]]),sdlog = sigma_mu))
    		        }
    		        else if (new_mu_law=="normal*t"){
    		          new_mu=rnorm(1,mean=mu$par[[i]],sd=sigma_mu*edge.length[x])
    		          while(new_mu<0){
    		            new_mu=rnorm(1,mean=mu$par[[i]],sd=sigma_mu*edge.length[x])
    		          }
    		          mu$par=c(mu$par,new_mu)
    		        }
    		      }else{
    		        rates=c(rates,i)
    		      }}
    		    
    		    #print(times)
    		    #print(nblineages)
    		  }
    		  i=i+1
    		}
    		i=1
    		while(i<=length(b) & continue){
    		  s=s+d[i]*sum((rates==i)[alive])/totalrate
    		  
    		  if(r<=s){    							
    		    if(length(which(rates==i & alive))>1){
    		      random_lineage = sample(which(rates==i & alive),1)
    		    }else{
    		      random_lineage = which(rates==i & alive)
    		    }
    		    continue=F
    		    edge.length[random_lineage] <- t - stem.depth[random_lineage]
    		    alive[random_lineage] <- FALSE
    		    times<-c(times,t)
    		    nblineages<-c(nblineages,sum(alive))}
    		  i=i+1
    		}
    	}}
    	
    	
  if (return.all.extinct == TRUE | sum(alive) > 0) {
    #print("return.tree")
    break
    }
    }
    							
	if (sum(alive)==0) {obj<-NULL}
 	else if (sum(alive)==1) {obj<-1}
 	else {
 	  edge.length[alive] <- t - stem.depth[alive]
 	  n <- -1
 	  for (i in 1:max(edge)) {
 	  if (any(edge[, 1] == i)) {
 	  edge[which(edge[, 1] == i), 1] <- n
 	  edge[which(edge[, 2] == i), 2] <- n
 	  n <- n - 1
 	  }
 	}
    
 	edge[edge > 0] <- 1:sum(edge > 0)
  tip.label <- 1:sum(edge > 0)
  mode(edge) <- "character"
  mode(tip.label) <- "character"
  obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label)

  class(obj) <- "phylo"
  obj <- old2new.phylo(obj)
   rep=rigth.order(obj,rates)
   obj=rep$tree
   rates=rep$rates
  if (prune.extinct){
    rep=prune.extinct.with.rates(obj,rates)
    obj=rep$tree
    rates=rep$rates
    }
    #obj<-drop.extinct(obj)}
  return(list("tree"=obj,"times"=times,"nblineages"=nblineages,"rates"=rates,"lamb"=lamb,"mu"=mu))
 	}
  }

prune.extinct.with.rates=function(phy,rates)
  #used in the simulation function to remove extinct taxa withour loosing the rate information
  
{
  obj=list(tree=phy,rates=rates)
  extinct=is.extinct(obj$tree,tol=max(obj$tree$edge.length)/1000)
  nodes=extinct
  if(length(extinct)>0){
    for(i in 1:length(extinct)){
      edge=which(obj$tree$edge[,2]==which(obj$tree$tip.label==extinct[i]))
      if(obj$tree$edge[edge,1]==(obj$tree$Nnode+2)){
        edge=which(obj$tree$edge[,1]==(obj$tree$Nnode+2))
      }else{
        edge=c(which(obj$tree$edge[,1]==obj$tree$edge[edge,1]))
      }
      obj$rates=obj$rates[-edge]
      obj$tree=drop.tip(obj$tree,extinct[i])
    }}
  return(obj)
}

rigth.order=function(phy,rates){
  n=phy$Nnode+1
  root=which(sapply(1:(phy$Nnode*2+1),function(x){!(x %in% phy$edge[,2])}))
  next_node=c(root)
  order=rep(0,2*n-1)
  is.tip=c()
  i=1
  while(length(next_node)>0){
    order[next_node[1]]=i
    offspring=phy$edge[phy$edge[,1]==next_node[1],2]
    next_node=c(offspring,next_node[-1])
    if(length(offspring)==0) is.tip=c(is.tip,i)
    i=i+1
  }
  phy$edge[,1]=order[phy$edge[,1]]
  phy$edge[,2]=order[phy$edge[,2]]
  order=order(phy$edge[,2])
  phy$edge=phy$edge[order,]
  rates=rates[order]
  phy$edge.length=phy$edge.length[order]
  itip=1
  inode=n+1
  newnames=c()
  for(i in 1:(2*n-1)){
    if(i %in% is.tip) {
      newnames=c(newnames,itip)
      itip=itip+1
    }else{
        newnames=c(newnames,inode)
        inode=inode+1
      }
  }
  phy$edge=matrix(newnames[phy$edge],ncol=2)
  return(list(tree=phy,rates=rates))
}

mean.rate=function(tree,rate){
  nh=nodeHeights(tree)
  times=sort(unique(c(nh[,1],nh[,2])))
  meanRate=c()
  nb=0
  keep.i=c()
  for(i in 1:(length(times)-1)){
    t=times[i]
    which.edge=(nh[,1]<=t & nh[,2]>t)
    if(sum(which.edge)>nb){
      nb=sum(which.edge)
    meanRate=c(meanRate,sum(rate[which.edge])/sum(which.edge))
    keep.i=c(keep.i,i)
    }
  }
  keep.i=c(keep.i,length(times))
  return(rbind(times[keep.i],c(meanRate[1],meanRate)))
}
# 
# X=mean.rate(tree,true.rate)
# plot.default(t(X), xaxs = "r", yaxs = "r",  type = "S")
