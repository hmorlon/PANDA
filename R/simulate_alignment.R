simulate_alignment <-
function(iter,name_index,name,host_tree,mu,n,seed,N,proportion_variant,simul,model,...){
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
    tree <- output[[1]]
    switches <- output[[2]]
    write.tree(tree,file=paste("simulations/symbiont_tree_",name,"_",index,".tre",sep=""))
    tree$edge.length <- tree$edge.length/sum(tree$edge.length) # scaled only for the traceable switches research 
    tips_alignment <- host_tree$tip.label
    tree <- read.tree(file=paste("simulations/symbiont_tree_",name,"_",index,".tre",sep=""))
    invisible(capture.output(plot_simulated_switches(n=n,host_tree=host_tree,name=name,index=index,switches=switches,path=path)))
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
  PI <- c(0.25,0.25,0.25,0.25)
  Q <-  (as.matrix(rbind(c(0,a,b,c)*t(PI),c(a,0,d,e)*t(PI),c(b,d,0,f)*t(PI),c(c,e,f,0)*t(PI))))
  diag(Q) <- - apply(Q,1,sum)
  Q <- -Q/Q[4,4] 
  eigQ <- eigen(Q)
  ivp <- solve(eigQ$vectors)
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
  row.names(variant_sequences) <- gsub(" ","",row.names(variant_sequences))
  for (missing in setdiff(tree$tip.label,row.names(variant_sequences))){tree <-drop.tip(tree,missing)}
  n <- nrow(variant_sequences)
  for (i in N:1) {if (length(unique(variant_sequences[,i]))==1){variant_sequences<-variant_sequences[,-i,drop=F]}}
  N_invariant <- N - ncol(variant_sequences)
  N_variant <- ncol(variant_sequences)
  variant_sequences <-  variant_sequences[,colSums(variant_sequences=='-') < n-1,drop=F]
  variant_sequences <- variant_sequences[tree$tip.label,,drop=F]
  variant_sequences <- rbind(variant_sequences, rep(0, ncol(variant_sequences)))
  duplicated <- duplicated(variant_sequences, MARGIN = 2)
  for (i in which(!duplicated)) {variant_sequences[n+1,i] <- sum(apply(variant_sequences, 2, identical, variant_sequences[,i]))}
  variant_sequences <- variant_sequences[,!duplicated,drop=F]
  variant_sequences <- rbind(variant_sequences[tree$tip.label,,drop=F], variant_sequences[n+1,])
  nodes <- n+order(node.depth.edgelength(tree)[(n+1):(2*n-1)],decreasing =T)
  simulated_likelihood <- LL(mu=mu,symbiont_tree=tree,nodes=nodes,Nd=ncol(variant_sequences),n=n,eigQ=eigQ,ivp=ivp,sequences=variant_sequences,PI=PI)
  simulated_likelihood <- simulated_likelihood + sum(as.numeric(variant_sequences[nrow(variant_sequences),]))*log(1-exp(-mu*sum(tree$edge.length)))
  save(simulated_likelihood,simulated_mu,seed,simul,proportion_variant,file=paste("data/simulation_data_",name,"_",index,".RData",sep=""))
}
