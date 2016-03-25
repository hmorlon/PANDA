library(geiger)
library(ape)

birthdeath.tree.rateshift <- function (theta,lamb_par, mu_par,
                                       f.lamb=function(x,y){exp(-alpha*x)*y}, f.mu=function(x,y){y}, lamb_shift=0,
                                       new_lamb_law="uniform",new_mu_law="uniform",alpha=0, sigma=0.1,lamb_max=1,lamb_min=0,mu_min=0,mu_max=0, time.stop = 0, 
                                       sigma_mu=0,taxa.stop = 0, return.all.extinct=TRUE, prune.extinct=TRUE, condition="time")
# start with only 1 lineage, not 2.
# (not yet, for the moment constant rate only, but can easily be implemented) allow rate variation in time as specified by f.lamb, f.mu and the parameters lamb_par and mu_par, but 
# the parameters of the time variation are given from the past to the present
# extinct lineages can be pruned. Also return number of lineages through time.
# Note: the resulting tree has no root length. This root lenght is given by the first time when there is an event, i.e. $times[2] 
#At a speciation event, each new lineage get a new rate with probability theta. This new rate is uniformly 
#drawn in [0, lamb_max]

#new_lamb_law can be uniform (between 0 and lamb_max), normal or lognormal (parameter sigma)
#If we want the same model as what we have in the likelihood, we need theta=1, new_lamb_law="normal", mu_par=mu_max=0
{
  lamb=list(fun=f.lamb,par=lamb_par)
  mu=list(fun=f.mu,par=mu_par)
  rates=c()
  
  while (1) {
    
    nblineages<-c(1)
    times<-c(0)
    b<-lamb$fun(0,lamb$par)
    d<-mu$fun(0,mu$par)
    dt <- rexp(1,(b + d))
    #print(dt)
    t<-dt
    
    if (t >= time.stop & condition=="time") {
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
        if(u[i]<theta){
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
    	  if (sum(alive)==taxa.stop & condition=="taxa") {
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
    		      if(u[j]<theta){
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
    		        }
    		        else if (new_lamb_law=="normal*t"){
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


