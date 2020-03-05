independent_evolution <-
function(replicate,name,index,seed,nb_tree,lambda,nb_cores){
  set.seed(seed+replicate)
  load(paste("data/data_model_",name,"_",index,".RData",sep=""))
  index_randomize <- sample(rownames(variant_sequences))
  index_randomize <- rbind(rownames(variant_sequences),index_randomize)
  rownames(index_randomize)[1] <- "original_index"
  rownames(variant_sequences) <- index_randomize[2,]
  write.table(index_randomize, paste("model_selection/randomize_indexes_",name,"_",index,"_replicate_",replicate,".txt",sep=""), col.names=F,quote=F)
  
  output <- optimize(f=inference_vertical_transmission,lower=0.0001,upper=50,tol=0.05,name=name,index=index,sequences=variant_sequences)
  
  results <- cbind(0,output$minimum,output$objective)
  colnames(results) <- c("ksi","mu","minloglik")
  write.table(results, paste("model_selection/results_randomize_",name,"_",index,"_replicate_",replicate,".txt",sep=""), append=FALSE, row.names=F,quote = F,sep="\t")
  
  inference_ksi_indep <- function(ksi){
    load(file=paste("simulated_trees/simulated_trees_",name,"_",ksi,".RData",sep=""))
    if (n!=Ntip(list_tree[[1]])) {missing_symbiont <- setdiff(list_tree[[1]]$tip.label,rownames(variant_sequences))
    for(i in 1:nb_tree){for (missing in missing_symbiont){list_tree[[i]] <- drop.tip(list_tree[[i]], tip=missing)}}}
    
    output <- optimize(f=inference_switches, lower=0.0001, upper=50,tol=0.05,randomize=T,ksi=ksi,name=name,index=index,sequences=variant_sequences,nb_tree=nb_tree,list_tree=list_tree,eig_val=eig_val, eig_vect=eig_vect, ivp=ivp, propinv=propinv)
    
    results <- cbind(ksi,output$minimum,output$objective)
    colnames(results) <- c("ksi","mu","minloglik")
    write.table(results, paste("model_selection/results_randomize_",name,"_",index,"_replicate_",replicate,".txt",sep=""), append=TRUE, col.names=FALSE, row.names=F,quote = F,sep="\t")
  }
  
  output <- mclapply(lambda,inference_ksi_indep,mc.cores=nb_cores)
  
  results <- read.table(paste("model_selection/results_randomize_",name,"_",index,"_replicate_",replicate,".txt",sep=""), header=TRUE,sep="\t")
  results <- cbind(c(replicate),results[which.min(results$minloglik),])
  colnames(results) <- c("replicate","ksi","mu","minloglik")
  options(warn=-1)
  write.table(results, paste("results/model_selection_independent_",name,"_",index,".txt",sep=""), row.names=F,quote = F,sep="\t",append=T,col.names=!file.exists(paste("results/model_selection_independent_",name,"_",index,".txt",sep="")))
  options(warn=0)
  
}
