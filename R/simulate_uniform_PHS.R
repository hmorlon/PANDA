simulate_uniform_PHS <-
function(name,est_ksi,nb_tree, host_tree, geo=FALSE, host=FALSE, stochastic_map=NULL, nb_cores=1){
  
  if(!exists("path")) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  dir.create(file.path(path, "PHS/"), showWarnings = FALSE)
  
  print(est_ksi)
  
  ####  Simulate uniform switches
  simul_uniform_distribution_tree <-function(tree,est_ksi,nb_tree=10000,seed=1){
    tree$edge.length <- tree$edge.length/sum(tree$edge.length)
    set.seed(seed+1)
    maxlen <- max(node.depth.edgelength(tree))
    uniform_switches <- vector("list", nb_tree)
    for (i in (1:nb_tree)){out_trees <- tree_change(tree,est_ksi,maxlen)
    uniform_switches[[i]] <- out_trees[[2]]}
    
    resampled_trees <- sample(seq(1:nb_tree), size=nb_tree, replace = T)
    uniform_switches <- uniform_switches[resampled_trees]
    
    uniform_switches <- data.frame(matrix(unlist(uniform_switches),nrow=3))
    return(uniform_switches)
  }
  
  simulate_tree_uniform_PHS <- function(j,name,est_ksi,nb_tree, host_tree, geo, host, stochastic_map,simul_uniform_distribution_tree){
    
    # simulate switches
    table_simul_switches <- simul_uniform_distribution_tree(tree=host_tree,est_ksi=est_ksi,nb_tree=nb_tree, seed=j)
    table_simul_switches <- data.frame(t(table_simul_switches[c(1,2,3),]))
    colnames(table_simul_switches) <- c("branch_depature","branch_arrival","position")
    
    # infer ratio geography
    if (geo==TRUE){
      
      table_simul_switches$geography <- NA
      
      for (i in 1:nrow(table_simul_switches)){
        
        branch <- table_simul_switches$branch_depature[i]
        b <- names(stochastic_map$maps[[branch]])[min(which(cumsum(stochastic_map$maps[[branch]])>table_simul_switches$position[i] - node.depth.edgelength(host_tree)[host_tree$edge[branch,1]] ))]  
        
        branch <- table_simul_switches$branch_arrival[i]
        bs <- names(stochastic_map$maps[[branch]])[min(which(cumsum(stochastic_map$maps[[branch]])>table_simul_switches$position[i] - node.depth.edgelength(host_tree)[host_tree$edge[branch,1]] ))]
        
        if (b==bs){table_simul_switches$geography[i] <- "intra"}else{table_simul_switches$geography[i] <- "inter"}
      }
      
      ratio_simul_geo <- length(which(table_simul_switches$geography=="intra"))/length(table_simul_switches$geography) 
    }
    
    # infer mean host-switch phylogenetic distances
    if (host==TRUE){
      table_simul_switches$host_relatedness <- NA
      
      maxlen <- max(node.depth.edgelength(host_tree))
      branches <- 1:(2*Ntip(host_tree)-2)
      tip_label_branch <- c()
      for (branch in branches){
        tip_label <-  get_tip_label(tree=host_tree,branch=branch) 
        tip_label_branch <- rbind(tip_label_branch,cbind(rep(branch,length(tip_label)),tip_label))}
      tip_label_branch <- data.frame(tip_label_branch)
      colnames(tip_label_branch) <- c("branch","tip_label")
      
      compute_host_relatedness_uniform <- function(switch){
        sliced_tree <- host_tree
        sliced_sub_trees <- phytools::treeSlice(sliced_tree,slice=table_simul_switches$position[switch], trivial=TRUE)
        for (i in 1:length(sliced_sub_trees)){if (Ntip(sliced_sub_trees[[i]])>1){
          sliced_tree <- drop.tip(sliced_tree,tip=sliced_sub_trees[[i]]$tip.label[2:Ntip(sliced_sub_trees[[i]])])}}
        for (i in which(node.depth.edgelength(sliced_tree)>table_simul_switches$position[switch])){sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)] <- sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)]-(maxlen-table_simul_switches$position[switch])}
        
        dist_nodes <- cophenetic.phylo(sliced_tree)
        rep1 <- intersect(tip_label_branch$tip_label[which(tip_label_branch$branch==table_simul_switches$branch_depature[switch])],sliced_tree$tip.label)
        rep2 <- intersect(tip_label_branch$tip_label[which(tip_label_branch$branch==table_simul_switches$branch_arrival[switch])],sliced_tree$tip.label)
        
        return(dist_nodes[rep1,rep2])
      }
      table_simul_switches$host_relatedness <- sapply(1:nrow(table_simul_switches), function(switch) compute_host_relatedness_uniform(switch))
      
      mean_simul_host <-  mean(c(table_simul_switches$host_relatedness))
    }
    return(c(ratio_simul_geo,mean_simul_host))
  }
  output <- mclapply(1:100, simulate_tree_uniform_PHS, mc.cores=nb_cores,name=name,est_ksi=est_ksi,nb_tree=nb_tree, host_tree=host_tree, geo=geo, host=host, stochastic_map=stochastic_map,simul_uniform_distribution_tree=simul_uniform_distribution_tree)
  
  ratio_simul_geo <- sapply(output,"[[",1) 
  mean_simul_host <- sapply(output,"[[",2) 
  save(ratio_simul_geo,mean_simul_host,file=paste0("PHS/uniform_switches_simulations_",name,"_",est_ksi,".Rdata"))
}
