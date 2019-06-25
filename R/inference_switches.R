inference_switches <-
function(mu,ksi,randomize,name,index,sequences,nb_tree,list_tree,eig_val, eig_vect, ivp, propinv){

  N=ncol(sequences)
  n=length(sequences)
  
  output <- unlist(lapply(1:nb_tree, function (i) LL_tree(mu,symbiont_tree=list_tree[[i]],sequences,n,N,eig_val, eig_vect, ivp, propinv)))
  
  rm(list_tree)
  if (all(is.finite(output))) {
    max_likelihood <- max(-unlist(output))
    likelihood <- -(log(sum(exp(-unlist(output)-max_likelihood)))-log(nb_tree)+max_likelihood)
    if (randomize==F){write.table(data.frame(cbind(1:nb_tree,unlist(output))),file=paste("results/optim_ll_",name,"_",index,"_",ksi,".txt",sep=""),col.names = F,row.names = F,quote=F)}
  } else {likelihood <- Inf}
  return(likelihood)
}
