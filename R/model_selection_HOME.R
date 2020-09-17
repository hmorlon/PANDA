model_selection_HOME <-
function(index,name,nb_tree=10000,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),nb_cores=1,seed=1,nb_random=10,overwrite=TRUE,...){
  
  print(noquote(paste("Index: ",index,sep="")))
  if (!file.exists(paste("data/data_model_",name,"_",index,".RData",sep=""))) stop(print("Please start by running the previous steps of HOME (fit_HOME...)"))
  
  if (file.exists(paste("data/alignment_variant_",name,"_",index,".fas",sep=""))){
    N_variant <- NULL
    load(paste("data/data_model_",name,"_",index,".RData",sep=""))
    if (N_variant>0){
      
      # strict vertical transmission
      output <- selection_vertical_transmission(name,index)
      
      # independent evolutions
      results <- read.table(paste("results/results_",name,"_",index,".txt",sep=""),header=T)
      results <- results[order(results$ksi),]
      results <- results[which(is.finite(results$minloglik)),]
      est_ksi <- results$ksi[which.min(results$minloglik)]
      est_mu <- results$mu[which.min(results$minloglik)]
      
      run <- 1
      nb_run <- 0
      if (file.exists(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""))){
        if (overwrite==TRUE){ file.remove(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""))
        }else{ if (nrow(read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T))==nb_random) {
          run <- 0 }else{nb_run <- nrow(read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T))}}}
      if (run == 1){
        
        # goes up to nb_random unless there are already more than 0.05*nb_run non-rejection of the null hypothesis
        replicate <- nb_run+1
        pvalue <- 0
        
        while ((replicate<=nb_random)&(pvalue<=0.05)){
          print(noquote(paste("Replicate: ",replicate,sep="")))
          output <- independent_evolution(replicate,name,index,seed,nb_tree,lambda,nb_cores)
        replicate <- replicate + 1
        table <- read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T)
        pvalue <- table$mu
        
        if (replicate>10) {pvalue <- max(c(length(which(table$ksi<=est_ksi))/nb_random, length(which(round(table$mu,digits=10)<=round(est_mu,digits=10)))/nb_random))}
        }
      }else{print(noquote("overwrite==FALSE: model selection is not performed again"))}
    }}}
