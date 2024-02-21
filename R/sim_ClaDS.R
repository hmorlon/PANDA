
sim_ClaDS <- function (lambda_0, mu_0,
                       new_lamb_law="lognormal*shift",new_mu_law="turnover",
                       condition="time", time_stop = 0, taxa_stop = Inf,
                       sigma_lamb=0.1, alpha_lamb=1, lamb_max=1,lamb_min=0,
                       sigma_mu=0, alpha_mu=1, mu_min=mu_0,mu_max=mu_0, 
                       theta=1,nShiftMax=Inf,
                       return_all_extinct=FALSE,prune_extinct=TRUE,
                       maxRate=Inf)
{

  relative_death=mu_0
  if(new_mu_law=="turnover"){ mu_0=mu_0*lambda_0}
  
  while (1) {
    
    lamb=lambda_0
    mu=mu_0
    rates=c()
    nShift=0
    tooHigh=F
    
    nblineages<-c(1)
    times<-c(0)
    b<-lamb
    d<-mu
    dt <- rexp(1,(b + d))
    #print(dt)
    t<-dt
    
    if ((t >= time_stop| 1>taxa_stop) & condition=="time") {
      t <- time_stop
    	alive<-1
    	rates<-c(1)
    	times<-c(times,t)
    	nblineages<-c(nblineages,1)
    	break
    }
    
    if (taxa_stop==1 & condition=="taxa") {
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
          rates=c(rates,length(lamb)+1)
          if(new_lamb_law=="uniform"){
            lamb=c(lamb,runif(1,min=lamb_min,max=lamb_max))
          }else if (new_lamb_law=="normal"){
              new_lambda=rnorm(1,mean=lamb[[1]],sd=sigma_lamb)
              while(new_lambda<0){
                new_lambda=rnorm(1,mean=lamb[[1]],sd=sigma_lamb)
              }
              lamb=c(lamb,new_lambda)
          }else if (new_lamb_law=="lognormal"){
              lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]),sdlog = sigma_lamb))
          }else if (new_lamb_law=="lognormal*shift"){
            lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]*alpha_lamb),sdlog = sigma_lamb))
          }else if (new_lamb_law=="lognormal*t"){
            lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]),sdlog = sigma_lamb*t))
          }else if (new_lamb_law=="logbrownian"){
            lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]),sdlog = sigma_lamb*sqrt(t)))
          }else if (new_lamb_law=="normal*t"){
            new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma_lamb*t)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma_lamb*t)
            }
            lamb=c(lamb,new_lambda)
          }else if (new_lamb_law=="normal+shift"){
            new_lambda=rnorm(1,mean=lamb[[1]]+alpha_lamb,sd=sigma_lamb)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb[[1]]+alpha_lamb,sd=sigma_lamb)
            }
            lamb=c(lamb,new_lambda)
          }else if (new_lamb_law=="normal*shift"){
            new_lambda=rnorm(1,mean=lamb[[1]]*alpha_lamb,sd=sigma_lamb)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb[[1]]*alpha_lamb,sd=sigma_lamb)
            }
            lamb=c(lamb,new_lambda)
          }
          if(new_mu_law=="uniform"){
            mu=c(mu,runif(1,min=mu_min,max=mu_max))
          }else if (new_mu_law=="normal"){
            new_mu=rnorm(1,mean=mu[[1]],sd=sigma_mu)
            while(new_mu<0){
              new_mu=rnorm(1,mean=mu[[1]],sd=sigma_mu)
            }
            mu=c(mu,new_mu)
          }else if (new_mu_law=="diversify"){
            new_mu=lamb[length(lamb)]-relative_death
            mu=c(mu,max(0,new_mu))
          }else if (new_mu_law=="turnover"){
            
            new_mu=lamb[length(lamb)]*relative_death
            mu=c(mu,new_mu)
          }else if (new_mu_law=="lognormal"){
            mu=c(mu,rlnorm(1, meanlog = log(mu[[1]]),sdlog = sigma_mu))
          }else if (new_mu_law=="lognormal*shift"){
            mu=c(mu,rlnorm(1, meanlog = log(mu[[1]]*alpha_mu),sdlog = sigma_mu))
          }else if (new_mu_law=="normal*t"){
            new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*t)
            while(new_mu<0){
              new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*t)
            }
            mu=c(mu,new_mu)
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
    	  if (sum(alive)>=taxa_stop){#} & condition=="taxa") {
    	    if(length(lamb)<=1){
    	      b=max(lamb,.Machine$double.eps)
    	      d=mu
    	      totalrate=sum(alive)*(b+d)
    	    }else{
    	      b<-sapply(1:length(lamb),function(x){max(lamb[[x]],.Machine$double.eps)})
    	      d<-sapply(1:length(mu),function(x){mu[[x]]})
    	      totalrate=sum(b[rates][alive])+sum(d[rates][alive])}
    	    
    	    dt <- rexp(1, totalrate)
    	    t <- t + dt
    	    
    	    break
    	  }else{
    		
    	  if(length(lamb)<=1){
    	    b=lamb
    	    d=mu
    	    totalrate=sum(alive)*(b+d)
    	    live_rates=rates[alive]
    	 }else{
    	   live_rates=rates[alive]
    	   b=lamb
    	   d=mu
         totalrate=sum(b[live_rates])+sum(d[live_rates])}

        dt <- rexp(1, totalrate)
    		t <- t + dt
    		
    		if (t >= time_stop & condition=="time") {
    		  t <- time_stop
    			times<-c(times,t)
    			nblineages<-c(nblineages,sum(alive))
    			break
    		}
    		
    		if (any(lamb>maxRate)) {
    		  tooHigh=T
    		  break
    		}
    		
    		r <- runif(1)
    		s=0
    		continue=T
    		if(theta<1) {
    		  Walive=unique(live_rates)
    		}else{
    		    Walive=live_rates[order(b[live_rates],decreasing = T)]
    		  }
    		k=1
    		sumBirth=sum(b[live_rates])/totalrate
    		i=Walive[k]
    		if(r<=sumBirth){
    		while(i<=length(b) & continue){
    		  if(theta<1){
    		    s=s+b[i]*sum((live_rates==i))/totalrate
    		  }else{
    		    s=s+b[i]/totalrate
    		    }
    		  if(r<=s){
    		    continue=F
    		    #print("speciation")
    		    if(theta<1){if(length(which(rates==i & alive))>1){
    		      random_lineage = sample(which(rates==i & alive),1)
    		    }else{
    		      random_lineage = which(rates==i & alive)
    		    }}else{
    		      random_lineage = which(rates==i)
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
    		        rates=c(rates,length(lamb)+1)
    		        if(new_lamb_law=="uniform"){
    		          lamb=c(lamb,runif(1,min=lamb_min,max=lamb_max))
    		        }else if (new_lamb_law=="normal"){
    		          new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma_lamb)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma_lamb)
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }else if (new_lamb_law=="lognormal"){
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]),sdlog = sigma_lamb))
    		        }else if (new_lamb_law=="lognormal*shift"){
    		          #mean lamb[[i]] iff alpha_lamb==exp(-sigma_lamb^2/2)
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]*alpha_lamb),sdlog = sigma_lamb))
    		        }else if (new_lamb_law=="lognormal*t"){
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]),sdlog = sigma_lamb*edge.length[x]))
    		        }else if (new_lamb_law=="logbrownian"){
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]),sdlog = sigma_lamb*sqrt(edge.length[x])))
    		        }else if (new_lamb_law=="normal*t"){
    		          new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma_lamb*edge.length[x])
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma_lamb*edge.length[x])
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }else if (new_lamb_law=="normal+shift"){
    		          new_lambda=rnorm(1,mean=lamb[[i]]+alpha_lamb,sd=sigma_lamb)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]]+alpha_lamb,sd=sigma_lamb)
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }else if (new_lamb_law=="normal*shift"){
    		          new_lambda=rnorm(1,mean=lamb[[i]]*alpha_lamb,sd=sigma_lamb)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]]*alpha_lamb,sd=sigma_lamb)
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }
    		        if(new_mu_law=="uniform"){
    		          mu=c(mu,runif(1,min=mu_min,max=mu_max))
    		        }else if (new_mu_law=="normal"){
    		          new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu)
    		          while(new_mu<0){
    		            new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu)
    		          }
    		          mu=c(mu,new_mu)
    		        }else if (new_mu_law=="turnover"){
    		          new_mu=lamb[length(lamb)]*relative_death
    		          mu=c(mu,new_mu)
    		        }else if (new_mu_law=="diversify"){
    		          new_mu=lamb[length(lamb)]+relative_death
    		          mu=c(mu,new_mu)
    		        }else if (new_mu_law=="lognormal"){
    		          mu=c(mu,rlnorm(1, meanlog = log(mu[[i]]),sdlog = sigma_mu))
    		        }else if (new_mu_law=="lognormal*shift"){
    		          mu=c(mu,rlnorm(1, meanlog = log(mu[[i]]*alpha_mu),sdlog = sigma_mu))
    		        }
    		        else if (new_mu_law=="normal*t"){
    		          new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*edge.length[x])
    		          while(new_mu<0){
    		            new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*edge.length[x])
    		          }
    		          mu=c(mu,new_mu)
    		        }
    		      }else{
    		        rates=c(rates,i)
    		      }}
    		    
    		    #print(times)
    		    #print(nblineages)
    		  }
    		  k=k+1
    		  i=Walive[k]
    		}}else{
    		  s=sumBirth
    		}
    		k=1
    		i=Walive[k]
    		while(i<=length(b) & continue){
    		  if(theta<1){
    		    s=s+d[i]*sum((live_rates==i))/totalrate
    		  }else{
    		    s=s+d[i]/totalrate
    		    }
    		  
    		  if(r<=s){    							
    		    if(theta<1){if(length(which(rates==i & alive))>1){
    		      random_lineage = sample(which(rates==i & alive),1)
    		    }else{
    		      random_lineage = which(rates==i & alive)
    		    }}else{
    		      random_lineage = which(rates==i)
    		    }
    		    continue=F
    		    edge.length[random_lineage] <- t - stem.depth[random_lineage]
    		    alive[random_lineage] <- FALSE
    		    times<-c(times,t)
    		    nblineages<-c(nblineages,sum(alive))}
    		  k=k+1
    		  i=Walive[k]
    		}
    	}}
    	
    	
  if (return_all_extinct == TRUE | sum(alive) > 0) {
    #print("return.tree")
    break
    }
    }
    							
	if ((sum(alive)==0 & prune_extinct) | (length(nblineages)==2 & !prune_extinct)) {obj<-NULL; root_length=t} #
 	else if (sum(alive)==1 & prune_extinct) {obj<-list(nbTaxa=1,"maxRate"=tooHigh); root_length=t}
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
  if (prune_extinct){
    rep=prune_extinct.with.rates(obj,rates)
    obj=rep$tree
    rates=rep$rates
  }
   root_length=t-max(node.depth.edgelength(obj))

 	}
  return(list("tree"=obj,"times"=times,"nblineages"=nblineages,"rates"=rates,"lamb"=lamb,"mu"=mu,"maxRate"=tooHigh,"root_length"=root_length))
  
  }

prune_extinct.with.rates=function(phy,rates,extinct=NULL)
  #used in the simulation function to remove extinct taxa withour loosing the rate information
  
{
  obj=list(tree=phy,rates=rates)
  if(is.null(extinct)) extinct=is.extinct(obj$tree,tol=max(obj$tree$edge.length)/100000)
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
  root=which(sapply(1:max(phy$edge),function(x){!(x %in% phy$edge[,2])}))
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
  if(inherits(rates,"list")){
    rates=lapply(rates,function(x){x[order]})
  }else{
    rates=rates[order]
  }
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

mean_rate=function(tree,rate){
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
# X=mean_rate(tree,true.rate)
# plot.default(t(X), xaxs = "r", yaxs = "r",  type = "S")
