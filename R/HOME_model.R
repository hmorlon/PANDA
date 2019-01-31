HOME_model <-
function(name,name_index,nb_cores=1,seed=3,nb_tree=5000,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),raref=FALSE,empirical=TRUE,randomize=TRUE,nb_random=10,...){
  
  if(!exists("name")) stop("Please provide the name of the dataset ")
  if(!exists("name_index")) stop("Please provide the name of the different OTU alignments ")
  if(!exists("path")) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  setwd(path)
  
  if (!exists("path_alignment")){ path_alignment <- path}
  if (exists("provided_tree")){
    host_tree <- provided_tree
    if (!file.exists(paste("host_tree_",name,".tre",sep=""))){ write.tree(host_tree,file=paste("host_tree_",name,".tre",sep=""))}}
  if (!file.exists(paste("host_tree_",name,".tre",sep=""))) stop("Please provide the host tree (format .tre) in the working directory")
  if (!is.binary(read.tree(paste("host_tree_",name,".tre",sep="")))) stop("Please provide a binary host tree")
  if (!is.rooted(read.tree(paste("host_tree_",name,".tre",sep="")))) stop("Please provide a rooted host tree")
  if (!is.ultrametric(read.tree(paste("host_tree_",name,".tre",sep="")))) stop("Please provide an ultrametric host tree")
  
  print("Data preparation:")
  output <- mclapply(1:length(name_index), prepare_data_HOME, mc.cores=nb_cores,name=name,name_index=name_index,path=path,path_alignment=path_alignment)
  
  print("Tree bank simulations:")
  output <- mclapply(1:length(lambda),simul_bank_tree,mc.cores=nb_cores,name=name,nb_tree=nb_tree,lambda=lambda,seed=seed)
  
  print("Global inference:")
  for (index in name_index){output <- fit_HOME(index=index,name=name,nb_tree=nb_tree,lambda=lambda,nb_cores=nb_cores,raref=raref)}
  
  print("Initial output:")
  output <- mclapply(1:length(name_index), output_results_HOME, mc.cores=nb_cores,name=name,name_index=name_index,lambda=lambda,nb_tree=nb_tree,empirical=empirical,randomize=F,raref=raref)
  
  print("Model selection:")
  for (index in name_index){output <- model_selection_HOME(index=index,name=name,nb_tree=nb_tree,lambda=lambda,nb_cores=nb_cores,seed=seed,nb_random=nb_random)}
  
  print("Output:")
  output <- mclapply(1:length(name_index), output_results_HOME, mc.cores=nb_cores,name=name,name_index=name_index,lambda=lambda,nb_tree=nb_tree,empirical=empirical,randomize=T,raref=raref)
  
}
