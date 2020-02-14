tree_change_host_dependent <-
function(name, index,host_tree,ksi,maxlen,host_signal){
  maxlen <- max(node.depth.edgelength(host_tree))
  step_slice <- maxlen/10001
  density_switch_position <- c()
  position_cood <- seq(step_slice,maxlen,step_slice)[1:10000]
  for (position_branch in position_cood){
    sliced_tree <- host_tree
    if (position_branch<maxlen){
      #sliced_sub_trees <- phytools::treeSlice(sliced_tree,slice=position_branch, trivial=TRUE)
      sliced_sub_trees <- tree_slice(sliced_tree,slice=position_branch)
      for (i in 1:length(sliced_sub_trees)){if (Ntip(sliced_sub_trees[[i]])>1){
        sliced_tree <- drop.tip(sliced_tree,tip=sliced_sub_trees[[i]]$tip.label[2:Ntip(sliced_sub_trees[[i]])])}}
      for (i in which(node.depth.edgelength(sliced_tree)>position_branch)){sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)] <- sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)]-(maxlen-position_branch)}}
    sum_dist_tips <- 0
    dist_nodes <- cophenetic.phylo(sliced_tree)
    for (i in 1:Ntip(sliced_tree)){
      ##### EXPONENTIAL (if host_signal=0 ; then it is uniform)
      for(j in setdiff(1:Ntip(sliced_tree),i)){ sum_dist_tips <-  sum_dist_tips + exp(-host_signal*dist_nodes[i,j])}
    }
    sum_dist_tips <- sum_dist_tips/(Ntip(sliced_tree)-1)
    
    density_switch_position <- c(density_switch_position,sum_dist_tips)
  }
  density_switch_position <- density_switch_position/sum(density_switch_position)
  pdf(paste("simulations/expected_distribution_switches_position_host_dependent_",name,"_",index,".pdf",sep=""),width=6,height = 4)
  plot(density_switch_position,pch=19,col="#145a32",xlab="Position",ylab="Density")
  dev.off()
  write.table(density_switch_position,paste("simulations/expected_distribution_switches_position_host_dependent_",name,"_",index,".txt",sep=""),row.names = F,col.names = F)
  
  loop=0
  while (loop<ksi){
    loop=0
    switches <- matrix(0,nrow=3,ncol=ksi)
    z <- sort(sample(position_cood,size=ksi,prob=density_switch_position,replace=F))
    switches[3,] <- z 
    output <- host_tree
    for (ind_ksi in 1:ksi) { 
      i <- z[ind_ksi]
      node_length <- node.depth.edgelength(output)
      i <- i - (maxlen - max(node_length)) 
      branches <- intersect(which(node_length[output$edge[,1]]<i),which(node_length[output$edge[,2]]>i))
      
      tip_label_branch <- c()
      for (branch in branches){
        tip_label <-  get_tip_label(tree=output,branch=branch) 
        tip_label_branch <- rbind(tip_label_branch,cbind(rep(branch,length(tip_label)),tip_label))}
      tip_label_branch <- data.frame(tip_label_branch)
      colnames(tip_label_branch) <- c("branch","tip_label")
      
      sliced_tree <- output
      write.tree(sliced_tree,paste0("tree",name,index,".tre"))
      sliced_tree <- read.tree(paste0("tree",name,index,".tre"))
      file.remove(paste0("tree",name,index,".tre"))
      
      #sliced_sub_trees <- phytools::treeSlice(sliced_tree,slice=i, trivial=TRUE)
      sliced_sub_trees <- tree_slice(sliced_tree,slice=i)
      for (j in 1:length(sliced_sub_trees)){if (Ntip(sliced_sub_trees[[j]])>1){
        for (dropped_tip in sliced_sub_trees[[j]]$tip.label[2:Ntip(sliced_sub_trees[[j]])]){sliced_tree <- drop.tip(sliced_tree,tip=dropped_tip)}}}
      for (j in which(node.depth.edgelength(sliced_tree)>i)){sliced_tree$edge.length[which(sliced_tree$edge[,2]==j)] <- sliced_tree$edge.length[which(sliced_tree$edge[,2]==j)]-(max(node_length)-i)  }
      
      dist_nodes <- cophenetic.phylo(sliced_tree)
      dist_nodes <- exp(-host_signal*dist_nodes) 
      proba_switch <- c()
      for (j in 1:(ncol(dist_nodes)-1)){for (k in (j+1):(ncol(dist_nodes))){proba_switch <- rbind(proba_switch, c(rownames(dist_nodes)[k],colnames(dist_nodes)[j],dist_nodes[k,j]) )}}
      proba_switch <- data.frame(proba_switch)
      colnames(proba_switch) <- c("B1","B2","Proba")
      proba_switch$Proba <- as.numeric(as.character(proba_switch$Proba))
      proba_switch$Proba <- proba_switch$Proba/sum(proba_switch$Proba)
      
      combi_branches <- sample(1:nrow(proba_switch),size=1,prob = proba_switch$Proba)
      
      b <- as.character(proba_switch[combi_branches,sample(1:2,1)])
      bs <- setdiff(c(as.character(proba_switch[combi_branches,1]),as.character(proba_switch[combi_branches,2])),b)
      b <- as.numeric(as.character(tip_label_branch$branch[which(tip_label_branch$tip_label==b)]))
      bs <- as.numeric(as.character(tip_label_branch$branch[which(tip_label_branch$tip_label==bs)]))
      
      if (length(bs)==1){
        switches[1,ind_ksi] <- b
        switches[2,ind_ksi] <- bs 
        output <- tree_switch(output,i,b,bs,node_length)
        loop <- loop +1 }}}
  return(list(output,switches))}
