prepare_data_HOME <-
function(iter,name,name_index,provided_tree=NULL,...){ 
  if(!exists("path")) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  dir.create(file.path(path, "results/"), showWarnings = FALSE)
  dir.create(file.path(path, "data/"), showWarnings = FALSE)
  dir.create(file.path(path, "simulated_trees/"), showWarnings = FALSE)
  dir.create(file.path(path, "figures/"), showWarnings = FALSE)
  dir.create(file.path(path, "rarefactions/"), showWarnings = FALSE)
  dir.create(file.path(path, "model_selection/"), showWarnings = FALSE)
  index <- name_index[iter]
  print(noquote(paste("Index: ",index,sep="")))
  
  setwd(path)
  set.seed(1)
  
  if (!exists("path_alignment")){ path_alignment <- path}
  
  #### Step 1 : Load the host tree ####
  if (!is.null(provided_tree)){
    host_tree <- provided_tree
    if (!file.exists(paste("host_tree_",name,".tre",sep=""))){ write.tree(host_tree,file=paste("host_tree_",name,".tre",sep=""))}}
  if (!file.exists(paste("host_tree_",name,".tre",sep=""))) stop(print("Please provide the host tree (format .tre) in the working directory"))
  if (!is.binary(read.tree(paste("host_tree_",name,".tre",sep="")))) stop(print("Please provide a binary host tree"))
  if (!is.rooted(read.tree(paste("host_tree_",name,".tre",sep="")))) stop(print("Please provide a rooted host tree"))
  if (!is.ultrametric(read.tree(paste("host_tree_",name,".tre",sep="")))) stop(print("Please provide an ultrametric host tree"))
  
  host_tree <- read.tree(paste("host_tree_",name,".tre",sep=""))
  host_tree <- ladderize(host_tree)
  host_tree$edge.length <- host_tree$edge.length/sum(host_tree$edge.length)  #host tree scaled with total branch length=1
  
  if (length(which(host_tree$edge.length==0))>0) stop(print("Please provide an host tree with positive branch lengths"))
  
  #### Step 2 : Load the symbiont sequences ####
  if (!file.exists(paste(path_alignment,"/alignment_",name,"_",index,".fas",sep=""))) stop(print(paste("Please provide an nucleotidic alignment (format .fas) in the working directory (or in path_alignment/) for the index ",index,sep="")))
  variant_sequences <-  read.dna(paste(path_alignment,"/alignment_",name,"_",index,".fas",sep=""),format="fasta",as.character=T)
  if (length(which(!rownames(variant_sequences) %in% host_tree$tip.label))>0) stop(print(paste("Please provide an nucleotidic alignment with names of sequences matching the names of the tips of the host tree for the index ",index,sep="")))
  
  # replace "n" nucleotide by gaps
  variant_sequences[which(!variant_sequences %in% c("a","t","g","c","-"))] <- "-"
  
  n <- nrow(variant_sequences)
  if (n>=5){
    N <- ncol(variant_sequences)
    for (i in N:1) {if (length(unique(variant_sequences[,i]))==1){variant_sequences<-variant_sequences[,-i,drop=F]}}
    variant_sequences <-  variant_sequences[,colSums(variant_sequences=='-') < n-1,drop=F]
    
    if(ncol(variant_sequences)>0){
      for (k in ncol(variant_sequences):1){if(length(unique(variant_sequences[which(variant_sequences[,k]!='-'),k]))==1){variant_sequences <-  variant_sequences[,-(k),drop=F]}}
      N_invariant <- N - ncol(variant_sequences)
      N_variant <- ncol(variant_sequences)
      if(N_variant>0){
        write.dna(variant_sequences,paste("data/alignment_variant_",name,"_",index,".fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=N)
        
        #### Step 3 : Order the host tree ####
        for (missing in setdiff(host_tree$tip.label,rownames(variant_sequences))){host_tree <-drop.tip(host_tree,missing)}  #row.names(variant_sequences)[1:n]
        r <- n+1 # root
        
        #### Step 4 : Substitution model ####
        sequences_model <- read.phyDat(paste("data/alignment_variant_",name,"_",index,".fas",sep=""),format="fasta")
        model_test <- modelTest(sequences_model,model=c("K80","F81","HKY"),G=F,I=F) #JC 
        model_evo <- eval(get(model_test$Model[which.min(model_test$BIC)], attr(model_test, "env")), env=attr(model_test, "env"))
        selected_model <- model_test$Model[which.min(model_test$BIC)]
        propinv <- t(as.matrix(model_evo$bf))
        Q <-  matrix(rbind(t(c(0,model_evo$Q[1],model_evo$Q[2],model_evo$Q[3])*t(propinv)),t(c(model_evo$Q[1],0,model_evo$Q[4],model_evo$Q[5])*t(propinv)),t(c(model_evo$Q[2],model_evo$Q[4],0,model_evo$Q[6])*t(propinv)),t(c(model_evo$Q[3],model_evo$Q[5],model_evo$Q[6],0)*t(propinv))),nrow = 4,ncol=4)
        diag(Q) <- - apply(Q,1,sum)
        Q <- -Q/Q[4,4] 
        eig_val <- eigen(Q)$values
        eig_vect <- eigen(Q)$vectors
        ivp <- solve(eig_vect)
        
        save(host_tree,variant_sequences,selected_model,n,r,N,N_invariant,N_variant,propinv,Q,eig_val,eig_vect,ivp,path,file=paste("data/data_model_",name,"_",index,".RData",sep=""))
        
      }else{save(host_tree,n,N,N_invariant,N_variant,path,file=paste("data/data_model_",name,"_",index,".RData",sep=""))
        write.table(c("NO VARIATION WITHIN ALIGNMENT",N_variant),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)
        write.table(c("NO VARIATION WITHIN ALIGNMENT",N_variant),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)}
    }else{N_invariant <- N
    N_variant <- 0
    save(host_tree,n,N,N_invariant,N_variant,path,file=paste("data/data_model_",name,"_",index,".RData",sep=""))
    write.table(c("NO VARIATION WITHIN ALIGNMENT",0),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)
    write.table(c("NO VARIATION WITHIN ALIGNMENT",0),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)}
  }else{#save(host_tree,n,path,file=paste("data/data_model_",name,"_",index,".RData",sep=""))
    N <- NA
    N_invariant <- NA
    N_variant <- NA
    save(host_tree,n,N,N_invariant,N_variant,path,file=paste("data/data_model_",name,"_",index,".RData",sep=""))
    write.table(c("NOT ENOUGHT HOSTS",n),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)
    write.table(c("NOT ENOUGHT HOSTS",n),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)}
}
