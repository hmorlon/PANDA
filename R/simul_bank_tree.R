simul_bank_tree <-
function(nb_ksi,name,nb_tree,lambda,seed=1){
  if (!file.exists(paste("host_tree_",name,".tre",sep=""))) stop("Please provide the host tree (format .tre) in the working directory")
  tree <- read.tree(paste("host_tree_",name,".tre",sep=""))
  tree$edge.length <- tree$edge.length/sum(tree$edge.length)
  set.seed(seed+1)
  ksi <- lambda[nb_ksi]
  maxlen <- max(node.depth.edgelength(tree))
  print(noquote(paste("Ksi: ",ksi,sep="")))
  list_tree <- vector("list", nb_tree)
  list_switches <- vector("list", nb_tree)
  for (i in (1:nb_tree)){out_trees <- tree_change(tree,ksi,maxlen)
  list_tree[[i]] <- out_trees[[1]] 
  list_switches[[i]] <- out_trees[[2]]
  }
  save(list_tree,file=paste("simulated_trees/simulated_trees_",name,"_",ksi,".RData",sep=""))
  save(list_switches,file=paste("simulated_trees/simulated_switches_",name,"_",ksi,".RData",sep=""))
  rm(list_tree)
  rm(list_switches)
}
