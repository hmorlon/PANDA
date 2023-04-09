sim_microbiota <-
function(name, name_index, simul, mu=1, n=20, provided_tree=NULL, N=300, proportion_variant=0.1, model="uniform", mean=0.5, sd=0.01, host_signal=10, geo_signal=0, stochastic_map=NULL, path=NULL, seed=1, nb_cores=1, scale=TRUE, delta=0, figure=FALSE, ...){
  if(!exists("name")) stop(print("Please provide the name of the simulations"))
  if(!exists("simul")) {simul  <- c(rep(0,5),rep(1,5),rep(3,5),rep(5,5),rep(10,5),rep("indep",5))}
  if(!exists("name_index")) {name_index  <- sapply(1:length(simul), function(i) paste("S",i,sep=""))}
  if(is.null(path)) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  setwd(path)
  
  dir.create(file.path(path, "simulations/"), showWarnings = FALSE)
  dir.create(file.path(path, "data/"), showWarnings = FALSE)
  dir.create(file.path(path, "figures/"), showWarnings = FALSE)
  
  set.seed(seed)
  
  if (!class(provided_tree)=="phylo"){
    print("Simulation of an host tree")
    host_tree <- pbtree(n=n)
    host_tree$tip.label <- sample(paste(rep("H",n),1:n,sep=""))}else{
      if (!is.binary(provided_tree)) stop(print("Please provide a binary host tree"))
      if (!is.rooted(provided_tree)) stop(print("Please provide a rooted host tree"))
      if (!is.ultrametric(provided_tree)) stop(print("Please provide an ultrametric host tree"))
      host_tree <- provided_tree
      write.tree(host_tree,file=paste("host_tree_",name,".tre",sep=""))
    }
  
  host_tree <- ladderize(host_tree)
  
  if (scale==TRUE) {host_tree$edge.length <- host_tree$edge.length/sum(host_tree$edge.length)} #host tree scaled with total branch length=1
  write.tree(host_tree,file=paste("host_tree_",name,".tre",sep=""))
  host_tree <- read.tree(file=paste("host_tree_",name,".tre",sep=""))
  n <- Ntip(host_tree)
  output <- mclapply(1:length(name_index),simulate_OTU_alignment,mc.cores=nb_cores,host_tree=host_tree,name=name,seed=seed,name_index=name_index,mu=mu,n=n,N=N,proportion_variant=proportion_variant,simul=simul,model=model,mean=mean,sd=sd,host_signal=host_signal,geo_signal=geo_signal,stochastic_map=stochastic_map,path=path,delta=delta,figure=figure)
}
