inference_vertical_transmission <-
function(mu,name,index,sequences){
  
  load(paste("data/data_model_",name,"_",index,".RData",sep=""))
  sequences <- rbind(sequences[host_tree$tip.label,,drop=F], sequences[nrow(sequences),])
  nodes <- n+order(node.depth.edgelength(host_tree)[(n+1):(2*n-1)],decreasing =T)
  likelihood <- LL(mu,host_tree,nodes,Nd=ncol(sequences),n,eigQ,ivp,sequences,PI)
  likelihood <- likelihood + sum(as.numeric(sequences[nrow(sequences),]))*log(1-exp(-mu*sum(host_tree$edge.length)))
  return(likelihood)
}
