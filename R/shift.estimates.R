shift.estimates <- function(phylo, data, sampling.fractions, comb.shift,
                            models = c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR"),
                            backbone.option = "backbone2", multi.backbone = F, 
                            np.sub = 4, n.max = NULL, rate.max = NULL, Ncores = 1){
  env.func <- environment()
  options(echo = T)
  
  #### argument check ####
  if(!inherits(data, "data.frame")){
    stop("object \"data\" is not of class \"data.frame\"")
  }
  
  if(!inherits(phylo, "phylo")){
    stop("object \"phylo\" is not of class \"phylo\"")
  } else {
    phylo$node.label <- c(c(Ntip(phylo)+1):c(Ntip(phylo)+Nnode(phylo)))
  }
  
  if(!inherits(sampling.fractions, "data.frame")){
    stop("object \"sampling.fractions\" is not of class \"data.frame\"")
  }
  
  if(any("Species" %in% colnames(data)) == F){
    stop("No column named \"Species\" in the database.
         \nPlease rename the corresponding column with the name \"Species\".")
  }
  
  if(is(comb.shift)[1] != "list"){
    stop("Object \"comb.shift\" is incorrect.")
  } else {
    ALL_bck_comb <- comb.shift
  }
  
  if(any(!phylo$tip.label %in% data$Species)){
    cat("The following tips are not in the database.\n \n")
    cat(phylo$tip.label[!phylo$tip.label %in% data$Species], "\n")
    stop()
  }
  
  if(all(sapply(comb.shift, is.null)) & multi.backbone == T){
    multi.backbone <- F
    cat("\nWarnings: There is no combination with multiple backbones.\n")
  }
  
  if(!np.sub %in% 1:4 & np.sub != "no_extinction"){
    cat("\nArgument \"np.sub\" is incorrect.")
    stop()
  }
  
  if(!backbone.option %in% c("backbone1","backbone2")){
    stop("\"backbone.option\" argument is incorrect." )
  }
  
  if(!is.null(n.max) & !is.null(rate.max)){
    cat("\nArguments \"rate.max\" and \"n.max\" cannot be used together.")
    stop()
  }
  
  # Final list to return
  phylo <- ladderize(phylo, F) # Mandatory to match node.labels (to check)
  
  all_res <- rep(list(NULL),4)
  
  totalsp <- nrow(data) # number of species in datanomic database
  data_phylo <- data[data$Species %in% phylo$tip.label,]
  
  ####_____________ ####
  #### WHOLE TREE ####
  
  cat("\n Estimating shifts of diversification from phylogeny.\n")
  cat("-----------------------------------------------------------")
  if(!is.null(n.max)){
    cat("\nA constrain will be applied: maximum diversity value is fixed at", n.max, "\n")
  }
  
  if(!is.null(rate.max)){
    cat("\nA constrain will be applied: maximum rate value is fixed at", rate.max, "\n")
  }

  f1 <- Ntip(phylo)/totalsp
  cat("\n--- WHOLE TREE ---\n \n")
  
  cat("\n","Sampling fraction =",  paste0(Ntip(phylo), "/", nrow(data), " (",round(f1,3)*100," %)"), "\n") 
  
  res_phylo <- div.models(phylo = phylo, tot_time = max(node.age(phylo)$ages), f = f1,
                        cond = "crown", models = models, n.max = n.max, rate.max = rate.max)
  
  res_phylo[,-1] <- apply(res_phylo[,-1], 2, as.numeric)
  
  for(i in 1:nrow(res_phylo)){res_phylo$Parameters[i] <- 4-sum(is.na(res_phylo[i,]))}
  
  res_phylo <- res_phylo[,c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta")]
  
  all_res[[1]] <- res_phylo ; names(all_res)[[1]] <- "whole_tree"
  
  all_tested_nodes <- unique(c(unique(unlist(strsplit(names(comb.shift), "[.]"))), unique(unlist(comb.shift))))
  
  all_lineages <- unique(unlist(strsplit(names(ALL_bck_comb), "[.]")))
  
  ALL_clade_names <- rep(list(NULL), length(unique(unlist(all_tested_nodes))))
  
  for(pot_names in 1:length(ALL_clade_names)){
    ALL_clade_names[pot_names] <- list(phylo$tip.label[unlist(Descendants(phylo, as.numeric(all_tested_nodes[pot_names])))])
  }
  names(ALL_clade_names) <- all_tested_nodes
  
  ALL_nodes_ages <- as.data.frame(apply(data.frame(nodesID=names(branching.times(phylo)),ages=branching.times(phylo)), 2, as.numeric))
  
  ALL_clade_names1 <- ALL_clade_names[all_lineages]
  
  #### SUBCLADES ####
  cat("\n--- SUBCLADES ---\n \n")
  ## Subclades trees
  # Check whether all groups are monophyletic 
  subclade_check <- NULL
  for(subclade in 1:length(ALL_clade_names1)){
    if(!is.monophyletic(phy = phylo, tips = ALL_clade_names1[[subclade]])){
      subclade_check <- c(subclade_check, names(ALL_clade_names1[subclade]))
    }
  }
  if(!is.null(subclade_check)){
    cat("The following subclades are not monophyletic:")
    cat("",subclade_check,sep="\n  - ")
    stop()
  }
  
  phylo2 <- list()
  for(cla in 1:length(ALL_clade_names1)){
    phylo2[[cla]] <- subtree(phylo, unlist(ALL_clade_names1[names(ALL_clade_names1)[cla]]))
  }
  
  names(phylo2) <- names(ALL_clade_names1)
  
  node_order <- names(branching.times(phylo)[order(branching.times(phylo))])
  node_order <- node_order[node_order %in% names(ALL_clade_names1)]
  ALL_clade_names2 <- ALL_clade_names1[match(node_order, names(ALL_clade_names1))]
  
  for(cla in 1:length(ALL_clade_names2)){
    
    if(is.na(sampling.fractions$sp_tt[as.numeric(names(ALL_clade_names2)[cla])])){
      ancest_cla <- Ancestors(phylo, as.numeric(names(ALL_clade_names2)[cla]))
      cla_no_NA <- sampling.fractions$nodes[sampling.fractions$nodes %in% ancest_cla &
                                              !is.na(sampling.fractions$sp_in)]
      ancest_cla <- ancest_cla[ancest_cla %in% cla_no_NA]
      ancest_cla <- ancest_cla[ancest_cla %in% as.numeric(names(ALL_clade_names2))]
      
      ancest_cla <- ancest_cla[1]
      
      sampling.fractions$sp_tt[as.numeric(names(ALL_clade_names2)[cla])] <- sampling.fractions$sp_tt[ancest_cla]
      sampling.fractions$sp_in[as.numeric(names(ALL_clade_names2)[cla])] <- sampling.fractions$sp_in[ancest_cla]
    }
    
  }
  
  totalsp2 <- as.list(sampling.fractions$sp_tt[as.numeric(names(ALL_clade_names1))])
  phylosp2 <- as.list(sampling.fractions$sp_in[as.numeric(names(ALL_clade_names1))])
  names(totalsp2) <- sampling.fractions$nodes[as.numeric(names(ALL_clade_names1))]
  names(phylosp2) <- sampling.fractions$nodes[as.numeric(names(ALL_clade_names1))]
  
  f2 <- as.list(unlist(phylosp2)/unlist(totalsp2))
  # totalsp2 <- totalsp2[names(totalsp2) %in% "others" == F]
  
  final2<-list()
  # all parameters
  
  backbone <- rep(list(F),length(phylo2))
  spec_times <- unlist(list(rep(list(NULL),length(phylo2))),recursive = F)
  branch_times <- unlist(list(rep(list(NULL),length(phylo2))),recursive = F)
  tot_time2 <- unlist(list(rep(list(NULL),length(phylo2))),recursive = F)
  cond <- unlist(list(rep(list(NULL),length(phylo2))),recursive = F)
  
  # backbone selection
  if(backbone.option == "backbone1"){
    cond <- rep(list("stem"),length(phylo2))
    for (i in 1:length(phylo2)){
      tot_time2[[i]] <- max(node.age(phylo2[[i]])$ages) + phylo2[[i]]$root.edge
    }
  }
  if(backbone.option == "backbone2"){
    cond <- rep(list("crown"),length(phylo2))
    for (i in 1:length(phylo2)){
      tot_time2[[i]] <- max(node.age(phylo2[[i]])$ages)
    }
  }
  
  names(tot_time2) <- names(phylo2)
  
  names_phylo2 <- sampling.fractions[sampling.fractions$nodes %in% as.numeric(names(phylo2)),]
  names_phylo2 <- names_phylo2$data[match(names_phylo2$nodes, as.numeric(names(phylo2)))]
  
  names_phylo2 <- ifelse(is.na(names_phylo2), paste("at node", names(phylo2)), paste("for", names_phylo2))
  
  models.sub <- models
  if(np.sub == 1){
    models.sub <- models.sub[models.sub %in% "BCST"]
  }
  
  if(np.sub == 2){
    models.sub <- models.sub[models.sub %in% c("BCST", "BCST_DCST", "BVAR")]
  }
  
  if(np.sub == 3){
    models.sub <- models.sub[models.sub %in% c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR")]
  }
  
  if(np.sub == "no_extinction"){
    models.sub <- models.sub[models.sub %in% c("BCST", "BVAR")]
  }
  
  # subclades to parallelize
  #cl <- makeCluster(Ncores, type="SOCK")
  #clusterExport(cl, list("phylo2", "tot_time2", "backbone", "spec_times", "branch_times", "f2", "models", "div.models",
  #                       "fitLikelihood", "suppressWarnings", "try","getLikelihood", ".Phi", "integrate", ".Psi", ".Integrate"), envir = env.func)
  
  final2 <- lapply(seq_along(phylo2), function(i){
    cat(i, "/", length(phylo2))
    n_tree <- sampling.fractions$sp_in[sampling.fractions$nodes == as.numeric(names(phylo2)[i])]
    n_tot <- sampling.fractions$sp_tt[sampling.fractions$nodes == as.numeric(names(phylo2)[i])]
    cat("\n","Sampling fraction", names_phylo2[i], "=",  paste0(n_tree, "/", n_tot, " (",round(f2[[i]],3)*100," %)"), "\n") 
    
    results <- div.models(phylo = phylo2[[i]], tot_time = tot_time2[[i]], f = f2[[i]],
                                 cond = F, models = models.sub, n.max = n.max, rate.max = rate.max)
    results1 <- div.models(phylo2[[i]], tot_time2[[i]], f = f2[[i]],
                                  cond = cond[[i]], models = models.sub, n.max = n.max, rate.max = rate.max, verbose = F)
    
    results <- merge(results[, c("Models", "Parameters", "logL", "AICc")], results1[,c("Models", "Lambda", "Alpha", "Mu", "Beta")], by = "Models")
    
    # adding a parameter for the location of the shift (to modify for the printing)
    results$AICc <- 2 * -results$logL + 2 * (results$Parameters+1) + (2 * (results$Parameters+1) * ((results$Parameters+1) + 1))/(n_tree - (results$Parameters+1) - 1)
    
    
    results[,-1] <- apply(results[,-1], 2, as.numeric)
    results$Parameters <- results$Parameters+1
    
    results$delta_AICc <- results$AICc - min(results$AICc)
    results <- results[order(results$delta_AICc),]
    
    #final2[[i]] <- results
    return(results)
    
  })
  
  names(final2) <- names(phylo2)
  names_final2 <- paste(names(phylo2),lapply(final2,function(x) paste(x$Models[x$delta_AICc < 2], collapse = "/")), sep = "/")
  
  all_res[2] <- list(final2)
  names(all_res)[2] <- "subclades"
  
  # Best model by clades
  best_subclades <- lapply(final2, function(x) x[x$AIC == min(x$AICc),])
  best_subclades_df <- data.frame(matrix(unlist(best_subclades), nrow = length(best_subclades), byrow = T))
  names(best_subclades_df) <- names(best_subclades[[1]])
  best_subclades_df$Clades <- names(best_subclades)
  best_subclades_df <- best_subclades_df[,c("Clades", names(best_subclades[[1]]))]
  best_subclades_df <- best_subclades_df[,-length(best_subclades_df)]
  best_subclades_df[,-c(1,2)] <- apply(best_subclades_df[,-c(1,2)], 2, as.numeric)
  
  ####_____________ ####
  #### BACKBONES ####
  cat("\n--- BACKBONES ---\n")
  
  if(length(ALL_bck_comb) == length(ALL_bck_comb[sapply(ALL_bck_comb, is.null)])){
    multi.backbone <- F
  }
  
  if(multi.backbone == T){
    
    cat("\n MULTI-BACKBONE OPTION ACTIVATED \n \n", length(unlist(ALL_bck_comb, recursive = F))+length(ALL_bck_comb[sapply(ALL_bck_comb, is.null)]), "combinations will be compared. \n")
  } else {
    cat("\n", length(ALL_bck_comb),"combinations will be compared. \n")
  }
  
  ALL_backbones <- rep(list(NULL),length(ALL_bck_comb))
  best_backbones <- NULL
  ALL_final3 <- rep(list(NULL),length(ALL_bck_comb))
  names(ALL_final3) <- sapply(ALL_bck_comb, function(x) paste(paste(x, collapse = "."),"_bck",sep = ""))
  ALL_TOTAL <- NULL
  clades <- names(totalsp2)
  
  # Individual branching times ####
  # calculated for backbone2 option and then transform for backbone1
  ALL_branch_times_clades <- rep(list(NULL),length(all_tested_nodes))
  names(ALL_branch_times_clades) <- all_tested_nodes
  
  for(clade in 1:length(all_tested_nodes)){
    
    parental_node <- Ancestors(phylo, as.numeric(all_tested_nodes[clade]), type = "parent")
    if(parental_node == Ntip(phylo)+1){ # this clade is at the root
      root_clade <- 1
    } else {
      root_clade <- 0
    }
    branch_times_clade <- unlist(list(rep(list(NULL),1 + root_clade)),recursive = F)
    
    bt_cl <- as.numeric(c(all_tested_nodes[clade], parental_node))
    # if(is.na(bt_cl[2])){}
    branch_times_clade[1] <- list(bt_cl)
    
    if(root_clade == 1){
      branch_times_clade[1 + root_clade] <- list(c(Siblings(phylo, max(bt_cl)),min(bt_cl)))
    }  
    ALL_branch_times_clades[[clade]] <- branch_times_clade
  }
  
  # LOOP ON COMBINATIONS ####
  # Parallelization of backbone models
  cat("\nDiversification models are running: \n")
  
  cl <- parallel::makeCluster(Ncores, type="SOCK")
  # Check what comes form RPANDA
  clusterExport(cl, list("all_comb_models", "subtree", "multi.backbone", "ALL_branch_times_clades", "phylo", "ALL_bck_comb","Descendants", "Ancestors", "Siblings", "getMRCA",
                         "get.node.ages", "drop.tip", "all_tested_nodes", "ALL_backbones", "node.age", "totalsp", "totalsp2", "Ntip", 
                         "div.models", "fit_bd_backbone","fit_bd_backbone_c", "likelihood_bd_backbone", "n.max", ".Phi", "integrate", ".Psi", ".Integrate",
                         "branching.times", "ALL_clade_names", "sampling.fractions", "backbone.option", "models", "ALL_final3", "rate.max",
                         "Children", "extract.clade.ln", "expand.grid", "get.branching.nodes"), envir = env.func)
  
  ALL_final3 <- ParallelLogger::clusterApply(cl, seq_along(ALL_bck_comb), all_comb_models, progressBar = T, stopOnError = T)
  stopCluster(cl)
  
  names(ALL_final3) <- names(ALL_bck_comb)
  if(multi.backbone == T){
    cat("\n\n--- Comparison(s) of the", length(unlist(ALL_bck_comb, recursive = F))+length(ALL_bck_comb[sapply(ALL_bck_comb, is.null)]), "combinations ---\n")
  } else {
    cat("\n\n--- Comparison(s) of the", length(ALL_bck_comb), "combinations ---\n")
  }
  
  best_ALL_final3 <- ALL_final3
  
  cat("\n")
  for(to in 1:length(ALL_final3)){
    
    for(sub_to in 1:length(ALL_final3[[to]])){
      for(ss_to in 1:length(ALL_final3[[to]][[sub_to]])){
        
        ALL_final3[[to]][[sub_to]][[ss_to]]$delta_AICc <- ALL_final3[[to]][[sub_to]][[ss_to]]$AICc - min(ALL_final3[[to]][[sub_to]][[ss_to]]$AICc)
        ALL_final3[[to]][[sub_to]][[ss_to]] <- ALL_final3[[to]][[sub_to]][[ss_to]][order(ALL_final3[[to]][[sub_to]][[ss_to]]$AICc),]
        ALL_final3[[to]][[sub_to]][[ss_to]][,-1] <- apply(ALL_final3[[to]][[sub_to]][[ss_to]][,-1], 2, as.numeric)
        best_ALL_final3[[to]][[sub_to]][[ss_to]] <- ALL_final3[[to]][[sub_to]][[ss_to]][1,]
      }
    }
  }
  
  # __ SELECTION #####
  if(multi.backbone == T){
    ALL_TOTAL <- as.data.frame(matrix(ncol = 4, nrow = length(unlist(ALL_bck_comb, recursive = F))+length(ALL_bck_comb[sapply(ALL_bck_comb, is.null)])))
  } else {
    ALL_TOTAL <- as.data.frame(matrix(ncol = 4, nrow = length(ALL_bck_comb)))
  }
  
  names(ALL_TOTAL) <- c("Combination", "Parameters","logL","AICc")
  
  n_to <- 1
  for(to in 1:length(best_ALL_final3)){
    names_subclades <- gsub("_bck","",unlist(strsplit(names(best_ALL_final3)[[to]], "[.]")))
    
    for(sub_to in 1:length(best_ALL_final3[[to]])){
      
      ALL_TOTAL[n_to,-1] <- apply(do.call(rbind.data.frame, lapply(best_ALL_final3[[to]][[sub_to]], function(x) x[, c("Parameters", "logL", "AICc")])), 2, sum)
      
      for(n_sub in 1:length(names_subclades)){
        ALL_TOTAL[n_to,-1] <- ALL_TOTAL[n_to,-1] + best_subclades_df[best_subclades_df$Clades %in% names_subclades[n_sub], c("Parameters","logL","AICc")]
      }
      
      ALL_TOTAL$Combination[n_to] <- paste0(paste(names_subclades, collapse = "."), "/", paste(names(best_ALL_final3[[to]][sub_to]), collapse = "."))
      n_to <- n_to + 1 
    }
  }
  
  res_phylo1 <- cbind("whole_tree", res_phylo[res_phylo$AICc == min(res_phylo$AICc), c("Parameters","logL","AICc")])
  names(res_phylo1) <- names(ALL_TOTAL)
  
  ALL_TOTAL <- rbind(ALL_TOTAL, res_phylo1)
  
  ALL_TOTAL$delta_AICc <- ALL_TOTAL$AICc - min(ALL_TOTAL$AICc)
  ALL_TOTAL <- ALL_TOTAL[order(ALL_TOTAL$delta_AICc),]
  row.names(ALL_TOTAL) <- NULL
  if(multi.backbone == F){
    ALL_TOTAL$Combination <- gsub("/","", ALL_TOTAL$Combination)
  }
  
  best_ALL_TOTAL <- ALL_TOTAL[ALL_TOTAL$delta_AICc < 2,]
  
  all_res[3] <- list(ALL_final3)
  all_res[4] <- list(ALL_TOTAL)
  names(all_res)[c(3,4)] <- c("backbones", "total")
  
  if(!"whole_tree" %in% best_ALL_TOTAL$Combination){
    cat("\n A total of", nrow(best_ALL_TOTAL), "combination(s) got the best fit(s) (delta AICc < 2). \n")
  } else{
    if(best_ALL_TOTAL$Combination[1] == "whole_tree"){
      cat("\n No shift has been detected. \n")
    } else {
      if(which(best_ALL_TOTAL$Combination == "whole_tree")-1 == 1){
        cat("\n",which(best_ALL_TOTAL$Combination == "whole_tree")-1,"combination has been detected. \n")
      } else {
        cat("\n",which(best_ALL_TOTAL$Combination == "whole_tree")-1,"combinations have been detected. \n")
      }
    }
  }
  return(all_res)
}