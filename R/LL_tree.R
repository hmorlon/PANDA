LL_tree <-
function(j,mu,list_tree,N_variant,n,eigQ,ivp,sequences,PI){
  symbiont_tree <- list_tree[[j]]
  sequences <- rbind(sequences[symbiont_tree$tip.label,,drop=F], sequences[nrow(sequences),])
  nodes <- n+order(node.depth.edgelength(symbiont_tree)[(n+1):(2*n-1)],decreasing =T)
  fit <- 0
  fit <- min(Inf, LL(mu,symbiont_tree,nodes,Nd=ncol(sequences),n,eigQ,ivp,sequences,PI))
  fit <- fit + sum(as.numeric(sequences[nrow(sequences),]))*log(1-exp(-mu*sum(symbiont_tree$edge.length)))
  return(fit)
}
