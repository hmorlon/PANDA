inference_switches <-
function(mu,ksi,randomize,nb_cores,name,index,sequences,nb_tree,list_tree){
  load(file=paste("data/data_model_",name,"_",index,".RData",sep=""))
  #output <- unlist(mclapply(1:nb_tree,LL_tree,mc.cores=nb_cores,mu=mu,list_tree=list_tree,N_variant=N_variant,n=n,eigQ=eigQ,ivp=ivp,sequences=sequences,PI=PI))
  output <- unlist(lapply(1:nb_tree, function (i) LL_tree (i,mu=mu,list_tree=list_tree,N_variant=N_variant,n=n,eigQ=eigQ,ivp=ivp,sequences=sequences,PI=PI)))
  
  rm(list_tree)
  if (all(is.finite(output))) {
    max_likelihood <- max(-unlist(output))
    likelihood <- -(log(sum(exp(-unlist(output)-max_likelihood)))-log(nb_tree)+max_likelihood)
    if (randomize==F){write.table(data.frame(cbind(1:nb_tree,unlist(output))),file=paste("results/optim_ll_",name,"_",index,"_",ksi,".txt",sep=""),col.names = F,row.names = F,quote=F)}
  } else {likelihood <- Inf}
  return(likelihood)
}
