inference_vertical_transmission <-
function(mu,name,index,sequences){
  
  load(paste("data/data_model_",name,"_",index,".RData",sep=""))
  
  n=nrow(sequences)
  N=ncol(sequences)
  
  fit <- LL_tree(mu,host_tree,sequences,n,N,eig_val, eig_vect, ivp, propinv)
  
  return(fit)
}
