model_selection_HOME <-
function(index,name,nb_tree=10000,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),nb_cores=1,seed=1,nb_random=10,overwrite=TRUE,...){
  
  print(noquote(paste("Index: ",index,sep="")))
  if (!file.exists(paste("data/data_model_",name,"_",index,".RData",sep=""))) stop(print("Please start by running the previous steps of HOME (fit_HOME...)"))
  
  if (file.exists(paste("data/alignment_variant_",name,"_",index,".fas",sep=""))){
    load(paste("data/data_model_",name,"_",index,".RData",sep=""))
    if (N_variant>0){
      output <- selection_vertical_transmission(name,index)
      run <- 1
      if (file.exists(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""))){
        if (overwrite==TRUE){ file.remove(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""))
        }else{ if (nrow(read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T))==nb_random) {run <- 0 }}}
      if (run == 1){
        for (replicate in 1:nb_random){
          print(noquote(paste("Replicate: ",replicate,sep="")))
          output <- independent_evolution(replicate,name,index,seed,nb_tree,lambda,nb_cores)} 
        table <- read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T)
      }else{print(noquote("overwrite==FALSE: model selection is not performed again"))}
    }}}
