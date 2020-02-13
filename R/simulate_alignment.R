simulate_alignment <-
function(iter,name_index,name,host_tree,mu,n,seed,N,proportion_variant,simul,model,path=getwd(),mean=0.5,sd=0.01,host_signal=10,geo_signal=0,stochastic_map=NULL,...){
  if(!exists("path")) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  setwd(path)
  
  set.seed(seed+iter)
  index <- name_index[iter]
  
  if (simul[iter]=="indep"){ # independent simulations
    simulated_traceable_ksi <- simul[iter]
    ksi <- 0
    indep_tree <- pbtree(n=n)
    indep_tree$tip.label <- sample(host_tree$tip.label)
    indep_tree$edge.length <- indep_tree$edge.length/sum(indep_tree$edge.length)  #host tree scaled with total branch length=1
    write.tree(indep_tree,file=paste("independent_tree_",name,"_",index,".tre",sep=""))
    tree <- host_tree <- indep_tree
  } else {ksi <- as.numeric(simul[iter])} # normal simulations
  
  simulated_mu <- mu
  maxlen <- max(node.depth.edgelength(host_tree))
  if ((simul[iter]!="indep")&(ksi==0)){
    tree <- host_tree
    simulated_traceable_ksi <- ksi}
  
  if ((simul[iter]!="indep")&(ksi>0)){
    if (model=="uniform") {output <- tree_change(host_tree,ksi,maxlen)}
    if (model=="temporal") {output <- tree_change_temporal(host_tree,ksi,maxlen,mean=maxlen*mean,sd)}
    if (model=="host_dependent") {output <- tree_change_host_dependent(name,index,host_tree,ksi,maxlen,host_signal)}
    if (model=="geo_dependent") {output <- tree_change_geo_dependent(name,index,host_tree,ksi,geo_signal,stochastic_map)}
    
    tree <- output[[1]]
    switches <- output[[2]]
    write.tree(tree,file=paste("simulations/symbiont_tree_",name,"_",index,".tre",sep=""))
    
    tree <- read.tree(file=paste("simulations/symbiont_tree_",name,"_",index,".tre",sep=""))
    invisible(capture.output(plot_simulated_switches(n=n,host_tree=host_tree,name=name,index=index,switches=switches)))
    
    switches[1,] <- as.character(switches[1,])
    switches[2,] <- as.character(switches[2,])
    write.table(switches, paste("simulations/simulated_switches_",name,"_",index,".txt",sep=""),col.names=F,row.names = F) 
  }
  
  maxlen <- max(node.depth.edgelength(tree))
  a <- 1
  b <- 4
  c <- 1
  d <- 1
  e <- 4
  f <- 1
  propinv <- c(0.25,0.25,0.25,0.25)
  Q <-  (as.matrix(rbind(c(0,a,b,c)*t(propinv),c(a,0,d,e)*t(propinv),c(b,d,0,f)*t(propinv),c(c,e,f,0)*t(propinv))))
  diag(Q) <- - apply(Q,1,sum)
  Q <- -Q/Q[4,4]
  eigQ <- eigen(Q)
  eig_val <- eigen(Q)$values
  eig_vect <- eigen(Q)$vectors
  ivp <- solve(eig_vect)
  nodes<-order(node.depth.edgelength(tree),decreasing=F)[-1]
  original_sequences <-  matrix(0,nrow=n,ncol=N)
  for (nu in 1:N){
    if (runif(1,0,1)<proportion_variant){
      L <-  matrix(0,nrow=(2*n-1),ncol=4)
      L[n+1,sample(1:4,size=1)] <- 1
      for(i in nodes) {
        v <- tree$edge[which(tree$edge[,2]==i),1]
        t <- as.numeric(tree$edge.length[which(tree$edge[,2]==i)])
        proba_nu <-  L[v,]%*%eigQ$vectors%*%diag(exp(t*mu*eigQ$values))%*%ivp
        L[i,sample(1:4,size=1,prob=proba_nu)] <- 1}
      n_seq <- rep(c("a","c","g","t"),n)[as.logical(t(L[1:n,]))]
      original_sequences[,nu]<- t(t(n_seq))
    } else {original_sequences[1:n,nu] <- sample(c("a","c","g","t"),size=1)}}
  row.names(original_sequences) <- tree$tip.label
  write.dna(original_sequences,paste("alignment_",name,"_",index,".fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=N)
  
  ####  compute simulated likelihood 
  variant_sequences <-  read.dna(paste("alignment_",name,"_",index,".fas",sep=""),format="fasta",as.character=T)
  rownames(variant_sequences) <- gsub(" ","",rownames(variant_sequences))
  for (missing in setdiff(tree$tip.label,rownames(variant_sequences))){tree <-drop.tip(tree,missing)}
  n <- nrow(variant_sequences)
  for (i in N:1) {if (length(unique(variant_sequences[,i]))==1){variant_sequences<-variant_sequences[,-i,drop=F]}}
  N_invariant <- N - ncol(variant_sequences)
  N_variant <- ncol(variant_sequences)
  variant_sequences <-  variant_sequences[,colSums(variant_sequences=='-') < n-1,drop=F]
  variant_sequences <- variant_sequences[tree$tip.label,,drop=F]
  
  if (N_variant>0){
    simulated_likelihood <- LL_tree(mu,tree,variant_sequences,n=nrow(variant_sequences),N=ncol(variant_sequences),eig_val, eig_vect, ivp, propinv)
    save(simulated_likelihood,simulated_mu,seed,simul,proportion_variant,file=paste("data/simulation_data_",name,"_",index,".RData",sep=""))
  }else{save(simulated_mu,seed,simul,proportion_variant,file=paste("data/simulation_data_",name,"_",index,".RData",sep=""))}}
