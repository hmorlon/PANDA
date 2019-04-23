createLikelihood_ClaDS0 <- function(phylo, root_depth=0, relative=T){
  
  nbtips = Ntip(phylo)
  edge = phylo$edge
  nedges=(2*(nbtips - 1))
  edge.length = phylo$edge.length
  
  type = rep(NA, Ntip(phylo))
  parents = rep(NA, Ntip(phylo))
  
  for (i in 1:nedges) {
    parent=which(phylo$edge[,2]==phylo$edge[i,1])
    if(length(parent)==0){
      type[i] = 1
      parents[i] = 0
      if ((phylo$edge[i,2]<=nbtips)){type[i]=2}
    }else if (!(phylo$edge[i,2]<=nbtips)){
      type[i] = 3    
      parents[i] = parent
    }else{
      type[i] = 4       
      parents[i] = parent
    }
  }
  
  roots = which(type <3)
  internal = which(type ==3)
  parentsInternal = parents[type == 3]
  internalAndRoots = which(type <4)
  internalAndTerminal = which(type >2)
  parentsInternalAndTerminal = parents[type >2]
  
  ancestors=list()
  for(i in 1:nedges){
    if(type[i]<3){
      ancestors[[i]]=-1
    }else{
      ancestors[[i]]=c(parents[i],ancestors[[parents[i]]][ancestors[[parents[i]]]>0])
    }
  }
  
    #lambda[1] is the basal speciation rate of the tree, and lambda[i] (i>1) is speciationRate[i]-speciationRate[parents[i]]
  ll<-function(lambda, sigma,alpha,lamb_shift=0){
    lambda2=lambda
    if(relative){
      lambda2=.Call("relToAbs", lambda=as.numeric(lambda2), parents=as.integer(parents), length=as.integer(length(lambda2)),
                    PACKAGE = "RPANDA") 
    }else{
      lambda=.Call("absToRel", lambda=as.numeric(lambda2), parents=as.integer(parents), length=as.integer(length(lambda2)),
                   PACKAGE = "RPANDA") 
    }
    # speciation events - notSpecation contribution 
    log.Likelihood = .Call("loglik", lambda=as.numeric(lambda), lambda2=as.numeric(lambda2), 
                           sigma=as.numeric(sigma),alpha=as.numeric(alpha), times=as.numeric(edge.length), 
                           internalAndRoots=as.integer(internalAndRoots), nNodes=as.integer(nbtips-1), root_depth=as.numeric(root_depth),
                           PACKAGE = "RPANDA")
    
    return(log.Likelihood)
  }
  
  relToAbs=function(lambda){
    rep=.Call("relToAbs",lambda=as.numeric(lambda), parents=as.integer(parents), length=as.integer(length(lambda)),
              PACKAGE = "RPANDA")
    return(rep)
  }
  
  absToRel=function(lambda){
    rep=.Call("absToRel",lambda=as.numeric(lambda), parents=as.integer(parents), length=as.integer(length(lambda)),
              PACKAGE = "RPANDA")
    return(rep)
  }
  res = list(ll = ll,relToAbs=relToAbs,absToRel=absToRel)
  
  
  return(res)
  
}

get_rates <- function(phylo, lambda, log_scale = TRUE){
  
  Ancestors=NULL
  nbtips = Ntip(phylo)
  edge = phylo$edge
  nedges=(2*(nbtips - 1))
  edge.length = phylo$edge.length
  
  if(is.null(Ancestors)){
    type = rep(NA, Ntip(phylo))
    parents = rep(NA, Ntip(phylo))
    
    for (i in 1:nedges) {
    parent=which(phylo$edge[,2]==phylo$edge[i,1])
    if(length(parent)==0){
      type[i] = 1
      parents[i] = NA
    }else if (!(phylo$edge[i,2]<=nbtips)){
      type[i] = 2    
      parents[i] = parent
    }else{
      type[i] = 3       
      parents[i] = parent
    }
  }
  
  roots = which(type ==1)
  internal = which(type ==2)
  parentsInternal = parents[type == 2]
  internalAndRoots = which(type <3)
  internalAndTerminal = which(type >1)
  parentsInternalAndTerminal = parents[type >1]
  
  ancestors=list()
  for(i in 1:nedges){
    if(type[i]==1){
      ancestors[[i]]=-1
    }else{
      ancestors[[i]]=c(parents[i],ancestors[[parents[i]]][ancestors[[parents[i]]]>0])
    }
  }
  }else{ancestors=Ancestors}
  
  if(log_scale){
    lambda2=c(lambda[1],sapply(1:(nedges),function(i){lambda[i+1]+lambda[1]+sum(lambda[ancestors[[i]][ancestors[[i]]>0]+1])}))
  }else{
    lambda2=c(lambda[1],sapply(1:(nedges),function(i){lambda[i+1]*lambda[1]*prod(lambda[ancestors[[i]][ancestors[[i]]>0]+1])}))
  }
  return(lambda2)
  
}

get_ancestors <- function(phylo){
  
  
  nbtips = Ntip(phylo)
  edge = phylo$edge
  nedges=(2*(nbtips - 1))
  edge.length = phylo$edge.length
  
  type = rep(NA, Ntip(phylo))
  parents = rep(NA, Ntip(phylo))
  
  for (i in 1:nedges) {
    parent=which(phylo$edge[,2]==phylo$edge[i,1])
    if(length(parent)==0){
      type[i] = 1
      parents[i] = NA
    }else if (!(phylo$edge[i,2]<=nbtips)){
      type[i] = 2    
      parents[i] = parent
    }else{
      type[i] = 3       
      parents[i] = parent
    }
  }
    
    ancestors=list()
    for(i in 1:nedges){
      if(type[i]==1){
        ancestors[[i]]=-1
      }else{
        ancestors[[i]]=c(parents[i],ancestors[[parents[i]]][ancestors[[parents[i]]]>0])
      }
    }
  return(ancestors)
  
}
