LL_tree <-
function(mu,symbiont_tree,sequences,n,N,eig_val, eig_vect, ivp, propinv){
  
  symbiont_tree <- reorder(symbiont_tree, "postorder")
  sequences <- sequences[symbiont_tree$tip.label,,drop=F]
  symbiont_tree$edge.length <- mu*symbiont_tree$edge.length
  
  nodes <- unique(symbiont_tree$edge[, 1]) - 1
  edges <- symbiont_tree$edge[, 2] -1
  el <- symbiont_tree$edge.length
  
  fit <- LL(mu,symbiont_tree,sequences, n, nodes, edges, el, eig_val, eig_vect, ivp, propinv, N)
  fit <- fit + N*log(1-exp(-sum(symbiont_tree$edge.length)))
  return(fit)
}
