model_selection_HOME <-
function(index,name,nb_tree,lambda,nb_cores,seed,nb_random=10,...){
  print(noquote(paste("Index: ",index,sep="")))
  if (!file.exists("data/data_model_",name,"_",index,".RData",sep="")) stop("Please start by running the previous steps of HOME (fit_HOME...)")
  load(paste("data/data_model_",name,"_",index,".RData",sep=""))
  if (N_variant>0){
    output <- selection_vertical_transmission(name,index)
    if (file.exists(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""))){file.remove(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""))}
    for (replicate in 1:nb_random){
      print(noquote(paste("Replicate: ",replicate,sep="")))
      output <- independent_evolution(replicate,name,index,seed,nb_tree,lambda,nb_cores)}
    table <- read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T)
  }}
