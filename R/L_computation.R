L_computation <-
function(k,mu,n_tree,nodes,n,eigQ,ivp,sequences,PI) {
  for (i in names(which(sequences[,k]=='-'))){n_tree <- drop.tip(n_tree,tip=i)}
  n_n<-Ntip(n_tree)
  n_seq<-sequences[which(sequences[,k]!='-'),k]
  if(n!=n_n){nodes<-n_n+order(node.depth.edgelength(n_tree)[(n_n+1):(2*n_n-1)],decreasing=T)}
  L = matrix(0,nrow=(2*n_n-1),ncol=4)
  L[which(n_seq=='a'),1]<-1
  L[which(n_seq=='c'),2]<-1
  L[which(n_seq=='g'),3]<-1
  L[which(n_seq=='t'),4]<-1
  L[which(n_seq=='n'),1:4]<-0.25
  for(i in nodes) {
    v <- n_tree$edge[which(n_tree$edge[,1]==i),2]
    t <- as.numeric(n_tree$edge.length[which(n_tree$edge[,1]==i)])
    L[i,] <-  (L[v[1],]%*%eigQ$vectors%*%diag(exp(t[1]*mu*eigQ$values))%*%ivp) * (L[v[2],]%*%eigQ$vectors%*%diag(exp(t[2]*mu*eigQ$values))%*%ivp)}
  K <- as.numeric(n_seq[n_n+1])*log(PI%*%L[n_n+1,])
  return(K)
}
