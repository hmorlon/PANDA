simul_bank_tree <-
function(nb_ksi,name, provided_tree=NULL,nb_tree=10000,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),seed=1){
  
  if (is.null(provided_tree)){
    if (!file.exists(paste("host_tree_",name,".tre",sep=""))){stop(print("Please provide the host tree (format .tre) in the working directory"))}
    provided_tree <- read.tree(paste("host_tree_",name,".tre",sep=""))
  }
  tree <- provided_tree
  if (!is.binary(tree)) stop(print("Please provide a binary host tree"))
  if (!is.rooted(tree)) stop(print("Please provide a rooted host tree"))
  if (!is.ultrametric(tree)) stop(print("Please provide an ultrametric host tree"))
  
  
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
