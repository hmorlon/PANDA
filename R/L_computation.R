L_computation <-
function(k,mu,n_tree,sequences,n, nodes, edges, el, eig_val, eig_vect, ivp, propinv) {
  n_seq<-sequences[,k,drop=F]
  if ("-" %in% n_seq){
    gaps <- which(n_seq=='-')
    n_tree <- drop.tip(n_tree,tip=rownames(n_seq)[gaps])
    n_seq <- n_seq[-gaps,]
    n <- Ntip(n_tree)}
  L = matrix(0,nrow=(2*n-1),ncol=4)
  L[n_seq=='a',1]<-1
  L[n_seq=='c',2]<-1
  L[n_seq=='g',3]<-1
  L[n_seq=='t',4]<-1
  return(llpruning(nodes,edges, el, L, eig_val, eig_vect, ivp, propinv))
}
