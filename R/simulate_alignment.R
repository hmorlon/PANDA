simulate_alignment <-
function(tree, mu=0.5, N=300, proportion_variant=0.1){
  host_tree <- tree
  n <- Ntip(host_tree)
  
  simulated_mu <- mu
  maxlen <- max(node.depth.edgelength(host_tree))
  tree <- host_tree
  
  propinv <- c(0.25,0.25,0.25,0.25)
  maxlen <- max(node.depth.edgelength(tree))
  a <- 1
  b <- 4
  c <- 1
  d <- 1
  e <- 4
  f <- 1
  Q <-  (as.matrix(rbind(c(0,a,b,c)*t(propinv),c(a,0,d,e)*t(propinv),c(b,d,0,f)*t(propinv),c(c,e,f,0)*t(propinv))))
  diag(Q) <- - apply(Q,1,sum)
  Q <- -Q/Q[4,4]
  eigQ <- eigen(Q)
  eig_val <- eigQ$values
  eig_val_mu <- mu*eig_val
  eig_vect <- eigQ$vectors
  ivp <- solve(eig_vect)
  nodes<-order(node.depth.edgelength(tree),decreasing=F)[-1]
  original_sequences <-  matrix(0,nrow=n,ncol=N)
  list_nuc <- c("a","c","g","t")
  list_possible_nuc <- rep(list_nuc,n)
  for (nu in 1:N){
    if (runif(1,0,1)<proportion_variant){
      L <-  matrix(0,nrow=(2*n-1),ncol=4)
      L[n+1,sample(1:4,size=1)] <- 1
      for(i in nodes) {
        ind <- which(tree$edge[,2]==i)
        v <- tree$edge[ind,1]
        t <- as.numeric(tree$edge.length[ind])
        proba_nu <-  L[v,]%*%eig_vect%*%diag(exp(t*eig_val_mu))%*%ivp
        L[i,sample(1:4,size=1,prob=proba_nu)] <- 1}
      n_seq <- list_possible_nuc[as.logical(t(L[1:n,]))]
      original_sequences[,nu]<- t(t(n_seq))
    } else {original_sequences[1:n,nu] <- sample(list_nuc,size=1)}}
  row.names(original_sequences) <- tree$tip.label
  
  return(original_sequences)
}
