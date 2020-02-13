tree_change_geo_dependent <-
function(name, index,host_tree,ksi,geo_signal,stochastic_map){
  maxlen <- max(node.depth.edgelength(host_tree))
  stochastic_map_slice <- RPANDA::CreateClassObject(stochastic_map,rnd=9)
  
  #### compute the density 
  nb_slice <- length(stochastic_map_slice$spans)
  step_slice <- maxlen/10001
  step_density <- seq(step_slice,maxlen,step_slice)[1:10000]
  density_switches <- c()
  for (slice in 1:nb_slice){
    list_pop <- table(stochastic_map_slice$class.object[[slice]][,2])
    
    proba_slice <- sum(list_pop[which((list_pop-1)>0)]) +   geo_signal*sum(list_pop[which((list_pop-1)==0)])
    
    density_switches[intersect(which(step_density<=stochastic_map_slice$times[slice+1]), which(step_density>=stochastic_map_slice$times[slice]))] <- proba_slice
  }
  density_switches[which(density_switches=="NaN")] <- 0
  density_switch_position <- density_switches/sum(density_switches)
  pdf(paste("simulations/expected_distribution_switches_position_geo_dependent_",name,"_",index,".pdf",sep=""),width=6,height = 4)
  plot(density_switch_position,pch=19,col="#145a32",xlab="Position",ylab="Density")
  dev.off()
  write.table(density_switch_position,paste("simulations/expected_distribution_switches_position_geo_dependent_",name,"_",index,".txt",sep=""),row.names = F,col.names = F)
  
  # simulate switches
  loop=0
  while (loop<ksi){
    loop=0
    switches <- matrix(0,nrow=3,ncol=ksi)
    z <- sort(sample(step_density,size=ksi,prob=density_switch_position,replace=F))
    switches[3,] <- z 
    output <- host_tree
    for (ind_ksi in 1:ksi) { 
      i <- z[ind_ksi]
      node_length <- node.depth.edgelength(output)
      i <- i - (maxlen - max(node_length)) 
      branches <- intersect(which(node_length[output$edge[,1]]<i),which(node_length[output$edge[,2]]>i))
      
      geo_branches <- c() # geography on each branch
      for (branch in branches){
        geo_branches <- c(geo_branches, names(stochastic_map$maps[[branch]])[min(which(cumsum(stochastic_map$maps[[branch]])>z[ind_ksi]-node.depth.edgelength(host_tree)[host_tree$edge[branch,1]]))])
      }
      names(geo_branches) <- branches
      
      proba_branches <- c() # proba switch on each branch
      for (branch in 1:length(branches)){
        proba_branches[branch] <- 2*(length(which(geo_branches==geo_branches[branch]))-1)
        proba_branches[branch] <- proba_branches[branch] + 2*geo_signal*length(which(geo_branches!=geo_branches[branch]))
      }
      
      b <- sample(x=branches, size=1, prob=proba_branches)
      
      proba_bs <- c() # proba switch on each branch
      for (branch in 1:length(branches)){
        proba_bs[branch] <- 0
        if (names(geo_branches)[branch]!=b){
          if (geo_branches[branch]==geo_branches[which(names(geo_branches)==b)]){ proba_bs[branch] <- 1
          } else { proba_bs[branch] <- geo_signal } 
        }}
      names(proba_bs) <- branches
      bs <- sample(x=branches, size=1, prob=proba_bs)
      
      switches[1,ind_ksi] <- b
      switches[2,ind_ksi] <- bs 
      output <- tree_switch(output,i,b,bs,node_length)
      loop <- loop +1 }}
  return(list(output,switches))}
