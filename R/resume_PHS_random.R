resume_PHS_random <-
function(name,index,est_ksi,est_mu, nb_tree, host_tree, geo=FALSE, host=FALSE, stochastic_map=NULL, nb_cores=1, seed){
  if(est_ksi>0){
    transparent_theme_no <- theme(panel.grid = element_blank(),
                                  axis.line = element_line("black"),
                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                  plot.background = element_rect(fill = "transparent",colour = NA),
                                  legend.key=element_blank(),legend.background=element_blank())
    
    
    #### Compute likelihood random switches  (100 replicates)
    
    ratio_simul_geo <- c()
    mean_simul_host <- c()
    
    eig_val <- NULL
    eig_vect <- NULL
    eig_val <- NULL
    ivp <- NULL
    propinv <- NULL
    list_switches <- NULL
    
    for (replicate in 1:100){
      set.seed(seed+replicate)
      load(paste("data/data_model_",name,"_",index,".RData",sep=""))
      index_randomize <- sample(rownames(variant_sequences))#[1:n])
      index_randomize <- rbind(rownames(variant_sequences),index_randomize)
      rownames(variant_sequences) <- index_randomize[2,]
      
      # compute likelihood randomizations
      if (!file.exists(paste("infer_PHS/optim_ll_random_",name,"_",index,"_",replicate,".txt",sep=""))){
        load(file=paste("simulated_trees/simulated_trees_",name,"_",est_ksi,".RData",sep=""))
        if (n!=Ntip(list_tree[[1]])) {missing_symbiont <- setdiff(list_tree[[1]]$tip.label,rownames(variant_sequences))
        for(i in 1:nb_tree){for (missing in missing_symbiont){list_tree[[i]] <- drop.tip(list_tree[[i]], tip=missing)}}}
        N=ncol(variant_sequences)
        n=length(variant_sequences)
        
        output <- unlist(lapply(1:nb_tree, function (i) LL_tree(est_mu,symbiont_tree=list_tree[[i]],variant_sequences,n,N,eig_val, eig_vect, ivp, propinv)))
        write.table(data.frame(cbind(round(unlist(output),digits = 5))),file=paste("infer_PHS/optim_ll_random_",name,"_",index,"_",replicate,".txt",sep=""),col.names = F,row.names = F,quote=F)
      }
      
      trees_ll <- read.table(paste("infer_PHS/optim_ll_random_",name,"_",index,"_",replicate,".txt",sep=""))
      trees_ll <- cbind(1:length(trees_ll$V1),trees_ll$V1)
      colnames(trees_ll) <- c("index","minloglik")
      load(file=paste("simulated_trees/simulated_switches_",name,"_",est_ksi,".RData",sep=""))
      
      max_likelihood <- max(-trees_ll$minloglik)
      proba_trees <- exp(-trees_ll$minloglik-max_likelihood)*exp(max_likelihood)
      if (all(proba_trees==0)){proba_trees <- exp(-trees_ll$minloglik-max_likelihood)}
      
      list_tree <- sample(x = trees_ll$index, prob=proba_trees, size = nb_tree, replace = T)
      
      list_switches_selected <- list_switches
      for (i in 1:length(list_tree)){
        index_tree <- list_tree[i]
        list_switches_selected[[i]] <- rbind(rbind(list_switches[[index_tree]],rep(trees_ll$minloglik[index_tree],ncol(list_switches[[index_tree]]))),rep(index_tree,ncol(list_switches[[index_tree]])))}
      selected_switches <- matrix(unlist(list_switches_selected),nrow=5)    
      
      table_inferred_switches <- data.frame(t(selected_switches[c(1,2,3,4,5),]))
      colnames(table_inferred_switches) <- c("branch_depature","branch_arrival","position","likelihood","tree")
      
      table_inferred_switches <- aggregate(list(numdup=rep(1,nrow(table_inferred_switches))), table_inferred_switches, length)
      
      if (geo==TRUE){
        compute_geo_switches <- function(i)  {
          branch <- table_inferred_switches$branch_depature[i]
          
          b <- names(stochastic_map$maps[[branch]])[min(which(cumsum(stochastic_map$maps[[branch]])>table_inferred_switches$position[i] - node.depth.edgelength(host_tree)[host_tree$edge[branch,1]] ))]  #z[ind_ksi]-node.depth.edgelength(host_tree)[host_tree$edge[branch,1]]))])
          
          branch <- table_inferred_switches$branch_arrival[i]
          bs <- names(stochastic_map$maps[[branch]])[min(which(cumsum(stochastic_map$maps[[branch]])>table_inferred_switches$position[i] - node.depth.edgelength(host_tree)[host_tree$edge[branch,1]] ))]
          
          if (b==bs){return("intra")}else{return("inter")}
        }
        table_inferred_switches$geography <- sapply(1:nrow(table_inferred_switches), function(i) compute_geo_switches(i))
        
        # number of duplicates
        ratio_simul_geo <- c(ratio_simul_geo, sum(c(table_inferred_switches$numdup[which(table_inferred_switches$geography=="intra")]))/sum(c(table_inferred_switches$numdup)) )
        
      }
      
      if (host==TRUE){    
        ###  Compute the host relatedness:
        maxlen <- max(node.depth.edgelength(host_tree))
        branches <- 1:(2*Ntip(host_tree)-2)
        tip_label_branch <- c()
        for (branch in branches){
          tip_label <-  get_tip_label(tree=host_tree,branch=branch)
          tip_label_branch <- rbind(tip_label_branch,cbind(rep(branch,length(tip_label)),tip_label))}
        tip_label_branch <- data.frame(tip_label_branch)
        colnames(tip_label_branch) <- c("branch","tip_label")
        
        compute_host_relatedness <- function(switch){
          sliced_tree <- host_tree
          sliced_sub_trees <- tree_slice(sliced_tree,slice=table_inferred_switches$position[switch])
          for (i in 1:length(sliced_sub_trees)){if (Ntip(sliced_sub_trees[[i]])>1){
            sliced_tree <- drop.tip(sliced_tree,tip=sliced_sub_trees[[i]]$tip.label[2:Ntip(sliced_sub_trees[[i]])])}}
          for (i in which(node.depth.edgelength(sliced_tree)>table_inferred_switches$position[switch])){sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)] <- sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)]-(maxlen-table_inferred_switches$position[switch])}
          
          dist_nodes <- cophenetic.phylo(sliced_tree)
          rep1 <- intersect(tip_label_branch$tip_label[which(tip_label_branch$branch==table_inferred_switches$branch_depature[switch])],sliced_tree$tip.label)
          rep2 <- intersect(tip_label_branch$tip_label[which(tip_label_branch$branch==table_inferred_switches$branch_arrival[switch])],sliced_tree$tip.label)
          
          return(dist_nodes[rep1,rep2]) 
        }
        table_inferred_switches$host_relatedness <- sapply(1:nrow(table_inferred_switches), function(switch) compute_host_relatedness(switch))
        
        # number of duplicates
        mean_simul_host <- c(mean_simul_host, sum(c(table_inferred_switches$numdup*table_inferred_switches$host_relatedness))/sum(c(table_inferred_switches$numdup)) )
      }
    }
    
    write.table(ratio_simul_geo,paste0("infer_PHS/simul_test_geo_signal_",name,"_",index,".txt"),row.names = F,col.names = F,quote=F)
    write.table(mean_simul_host,paste0("infer_PHS/simul_test_host_signal_",name,"_",index,".txt"),row.names = F,col.names = F,quote=F)
    
    #### Resume estimated switches
    trees_ll <- read.table(paste("results/optim_ll_",name,"_",index,"_",est_ksi,".txt",sep=""))
    trees_ll <- cbind(1:length(trees_ll$V1),trees_ll$V1)
    colnames(trees_ll) <- c("index","minloglik")
    load(file=paste("simulated_trees/simulated_switches_",name,"_",est_ksi,".RData",sep=""))
    
    max_likelihood <- max(-trees_ll$minloglik)
    proba_trees <- exp(-trees_ll$minloglik-max_likelihood)*exp(max_likelihood)
    if (all(proba_trees==0)){proba_trees <- exp(-trees_ll$minloglik-max_likelihood)}
    
    list_tree <- sample(x = trees_ll$index, prob=proba_trees, size = nb_tree, replace = T)
    
    list_switches_selected <- list_switches
    for (i in 1:length(list_tree)){
      index_tree <- list_tree[i]
      list_switches_selected[[i]] <- rbind(rbind(list_switches[[index_tree]],rep(trees_ll$minloglik[index_tree],ncol(list_switches[[index_tree]]))),rep(index_tree,ncol(list_switches[[index_tree]])))}
    selected_switches <- matrix(unlist(list_switches_selected),nrow=5)    
    
    table_inferred_switches <- data.frame(t(selected_switches[c(1,2,3,4,5),]))
    colnames(table_inferred_switches) <- c("branch_depature","branch_arrival","position","likelihood","tree")
    
    table_inferred_switches <- aggregate(list(numdup=rep(1,nrow(table_inferred_switches))), table_inferred_switches, length)
    
    pvalue_geo <- NA
    pvalue_host <- NA
    
    #### Infer influence of geography 
    if (geo==TRUE){
      compute_geo_switches <- function(i)  {
        branch <- table_inferred_switches$branch_depature[i]
        
        b <- names(stochastic_map$maps[[branch]])[min(which(cumsum(stochastic_map$maps[[branch]])>table_inferred_switches$position[i] - node.depth.edgelength(host_tree)[host_tree$edge[branch,1]] ))]  #z[ind_ksi]-node.depth.edgelength(host_tree)[host_tree$edge[branch,1]]))])
        
        branch <- table_inferred_switches$branch_arrival[i]
        bs <- names(stochastic_map$maps[[branch]])[min(which(cumsum(stochastic_map$maps[[branch]])>table_inferred_switches$position[i] - node.depth.edgelength(host_tree)[host_tree$edge[branch,1]] ))]
        
        if (b==bs){return("intra")}else{return("inter")}
      }
      table_inferred_switches$geography <- sapply(1:nrow(table_inferred_switches), function(i) compute_geo_switches(i))
      
      # number of duplicates
      ratio_inferred <- sum(c(table_inferred_switches$numdup[which(table_inferred_switches$geography=="intra")]))/sum(c(table_inferred_switches$numdup)) 
      
      pdf(paste("infer_PHS/distribution_test_geo_signal_",name,"_",index,".pdf",sep=""),width=6,height=5) 
      print(ggplot(data.frame(ratio_simul_geo))+geom_histogram(aes(x=ratio_simul_geo),fill="#1b4f72")+
              geom_vline(xintercept = ratio_inferred,col=c("#d35400"),size=1.5) + transparent_theme_no + xlab("Ratio intra/inter populations") )
      dev.off()
      
      pvalue_geo=length(which(ratio_simul_geo>=ratio_inferred))/length(ratio_simul_geo) #p-value
      
      write.table(pvalue_geo,paste0("infer_PHS/p-value_test_geo_signal_",name,"_",index,".txt"),row.names = F,col.names = F,quote=F)
    }
    
    #### Infer infuence host relatedness
    if (host==TRUE){    
      
      ###  Compute the host relatedness:
      maxlen <- max(node.depth.edgelength(host_tree))
      branches <- 1:(2*Ntip(host_tree)-2)
      tip_label_branch <- c()
      for (branch in branches){
        tip_label <-  get_tip_label(tree=host_tree,branch=branch)
        tip_label_branch <- rbind(tip_label_branch,cbind(rep(branch,length(tip_label)),tip_label))}
      tip_label_branch <- data.frame(tip_label_branch)
      colnames(tip_label_branch) <- c("branch","tip_label")
      
      compute_host_relatedness <- function(switch){
        sliced_tree <- host_tree
        sliced_sub_trees <- tree_slice(sliced_tree,slice=table_inferred_switches$position[switch])
        for (i in 1:length(sliced_sub_trees)){if (Ntip(sliced_sub_trees[[i]])>1){
          sliced_tree <- drop.tip(sliced_tree,tip=sliced_sub_trees[[i]]$tip.label[2:Ntip(sliced_sub_trees[[i]])])}}
        for (i in which(node.depth.edgelength(sliced_tree)>table_inferred_switches$position[switch])){sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)] <- sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)]-(maxlen-table_inferred_switches$position[switch])}
        
        dist_nodes <- cophenetic.phylo(sliced_tree)
        rep1 <- intersect(tip_label_branch$tip_label[which(tip_label_branch$branch==table_inferred_switches$branch_depature[switch])],sliced_tree$tip.label)
        rep2 <- intersect(tip_label_branch$tip_label[which(tip_label_branch$branch==table_inferred_switches$branch_arrival[switch])],sliced_tree$tip.label)
        
        return(dist_nodes[rep1,rep2]) 
      }
      table_inferred_switches$host_relatedness <- sapply(1:nrow(table_inferred_switches), function(switch) compute_host_relatedness(switch))
      
      # number of duplicates
      mean_infer_host <-  sum(c(table_inferred_switches$numdup*table_inferred_switches$host_relatedness))/sum(c(table_inferred_switches$numdup))
      
      pdf(paste("infer_PHS/distribution_test_host_signal_",name,"_",index,".pdf",sep=""),width=6,height=5) 
      print(ggplot(data.frame(mean_simul_host))+geom_histogram(aes(x=mean_simul_host),fill="#1b4f72")+
              geom_vline(xintercept = mean_infer_host,col=c("#d35400"),size=1.5) + transparent_theme_no + xlab("Mean phylogenetic relatedness") )
      dev.off()
      
      pvalue_host=length(which(mean_simul_host<=mean_infer_host))/length(mean_simul_host) #p-value
      write.table(pvalue_host,paste0("infer_PHS/p-value_test_host_signal_",name,"_",index,".txt"),row.names = F,col.names = F,quote=F)
    }
    return(c(pvalue_host, pvalue_host))
  }else{print("No host-switches inferred")
    return(c(NA, NA))}
}
