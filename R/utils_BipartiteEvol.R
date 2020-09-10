######## Make a well conformed tree ################################

# phy a phylo object
# trait a list of vector of trait corresponding to each edge


rigth.order.BE=function(phy,trait){
  
  root=unique(phy$edge[,1])[which(sapply(unique(phy$edge[,1]),function(x){!(x %in% phy$edge[,2])}))]  # root of the phylogeny
  nextNode=c(root)                                  # next node
  order=rep(0,phy$Nnode+length(phy$tip.label))      # rank of each node in the tree traversal
  isTip=c()                                         # rank of the tips in the tree traversal
  i=1                                               # current rank
  label=rep(0,length(phy$tip.label))                # tip labels
  
  # fill order, isTip and label by traversing the tree (Deep First Search ?)
  while(length(nextNode)>0){
    order[nextNode[1]]=i
    offspring=phy$edge[phy$edge[,1]==nextNode[1],2]
    
    if(length(offspring)==0) {
      # then nextNode[1] is a tip
      isTip=c(isTip,i)
      label[length(isTip)]=phy$tip.label[nextNode[1]]
    }
    
    nextNode=c(offspring,nextNode[-1])            # update next node
    i=i+1
  }
  
  # rename the nodes
  phy$edge[,1]=order[phy$edge[,1]]
  phy$edge[,2]=order[phy$edge[,2]]
  phy$tip.label=label[1:length(isTip)]   # if the length was not correct in the first place or tips where badly labeled
  
  # reorder the edges
  order=order(phy$edge[,2])
  phy$edge=phy$edge[order,]
  
  # reorder trait values
  if(! is.null(trait)){
    for(i in 1:length(trait)){
      trait[[i]]=(trait[[i]])[order]
    }
  }
  
  # reorder edge lengths
  phy$edge.length=phy$edge.length[order]
  
  # rename the nodes according to their nature (tip <= ntip, internal nodes > ntip)
  iTip=1                               # current tip name
  iNode=length(phy$tip.label)+1        # current node name
  newNames=c()                         # vector of new names
  
  for(i in 1:(phy$Nnode+length(phy$tip.label))){
    if(i %in% isTip) {
      newNames=c(newNames,iTip)
      iTip=iTip+1
    }else{
      newNames=c(newNames,iNode)
      iNode=iNode+1
    }
  }
  
  phy$edge=matrix(newNames[phy$edge],ncol=2)   # change the names in edge
  
  return(list(tree=phy,trait=trait))
}



######## Compute number of mutation between tips ###################

# gen an object obtained with make.gen


compute.dist=function(gen, verbose=T){
  
  N=length(gen$tip.label)               # number of tips
  distance=matrix(0,nrow=N,ncol=N)      # initiate distance matrix
  
  #main loop
  for(i in 1:nrow(gen$edge)){
    
    if(gen$nMut[i]>0){
      #add the number of mutation between the individuals separated by this edge
      if (verbose) cat("\r",i,"\r")
      if(gen$edge[i,2]<=N){
        tip=as.integer(gen$tip.label[gen$edge[i,2]])
        distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
        distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
      }else{
        tip=as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)
        distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
        distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
      }
    }
  }  # end main loop
  
  distance=distance[as.integer(gen$tip.label),as.integer(gen$tip.label)]
  return(distance)
}



######## Remove taxa when there are traits #########################

# phy a phylo object
# nodes the labels of the tips to supress 
# traits a list of vector of traits corresopnding to the edges of the tree
# additiveTrait vector of indice of extensive traits (eg branch length, nmut...)


prune.with.traits=function(phy,nodes,traits,extensiveTrait=c()){
  
  D=length(traits)                   # dimension of trait space
  obj=list(tree=phy,trait=traits)    # returned object
  edges=c()                          # removed edges
  if(length(extensiveTrait)>0){
    intensiveTrait=(1:D)[-extensiveTrait]
  }else{
    intensiveTrait=1:D
  }
  
  if(length(nodes)>0){
    # main loop
    for(i in 1:length(nodes)){
      edge=which(obj$tree$edge[,2]==which(obj$tree$tip.label==nodes[i]))   # edge leading to nodes[i]
      
      while(edge %in% edges){
        # ie edge has already been supressed (all the sister clades of nodes[i] have ben supressed)
        # we have to supress the parent edge
        edge=which(obj$tree$edge[,2]==obj$tree$edge[edge,1])
      }
      
      sisterEdge=which(obj$tree$edge[,1]==obj$tree$edge[edge,1])  # sister edges of edge (including edge)
      sisterEdge=sisterEdge[which(!(sisterEdge %in% edges))]      # remaining sister edges of edge
      
      if(length(sisterEdge)<3){
        # in that case we supress edge and its sister edge
        edges=c(edges,sisterEdge)
        parent=which(obj$tree$edge[,2]==obj$tree$edge[edge,1])
        
        if(length(extensiveTrait)>0){
          # for an extensive trait we have to sum the trait of the two combined edges
          for(j in extensiveTrait){
            obj$trait[[j]][parent]=obj$trait[[j]][sisterEdge[sisterEdge!=edge]]+obj$trait[[j]][parent]
          }
        }
        
      }else{
        # in that case we only have to supress edge, no edges are combined
        edges=c(edge,edges)
      }
    }  # end main loop
    
    for(j in 1:D){
      obj$trait[[j]]=obj$trait[[j]][-edges]
    }
    
    obj$tree=drop.tip(obj$tree,which(obj$tree$tip.label %in% nodes))
  }
  return(obj)
}

