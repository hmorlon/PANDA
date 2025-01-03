shift.estimates <- function(phylo, data, sampling.fractions, comb.shift,
                            models = c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR"),
                            backbone.option = "crown.shift", multi.backbone = FALSE, 
                            np.sub = 4, rate.max = NULL, n.max = NULL, Ncores = 1){
  
  # reset global options when exiting the function
  old_options <- options()
  on.exit(options(old_options))
  
  env.func <- environment()
  options(echo = TRUE)
  
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
  
  if(any("Species" %in% colnames(data)) == FALSE){
    stop("No column named \"Species\" in the database.
         \nPlease rename the corresponding column with the name \"Species\".")
  }
  
  if(any(!phylo$tip.label %in% data$Species)){
    message("The following tips are not in the database.\n \n")
    message(phylo$tip.label[!phylo$tip.label %in% data$Species], "\n")
    stop()
  }
  
  if(!np.sub %in% 1:4 & np.sub != "no_extinction"){
    stop("\nArgument \"np.sub\" is incorrect.")
  }
  
  if(!backbone.option %in% c("stem.shift","crown.shift")){
    stop("\"backbone.option\" argument is incorrect." )
  }
  
  if(!is.null(n.max) & !is.null(rate.max)){
    stop("\nArguments \"rate.max\" and \"n.max\" cannot be used together.")
  }
  
  # Final list to return
  phylo <- ladderize(phylo, F) # Mandatory to match node.labels (to check)
  
  all_res <- rep(list(NULL),4)
  
  totalsp <- nrow(data) # number of species in datanomic database
  data_phylo <- data[data$Species %in% phylo$tip.label,]
  
  # sorting comb.shift
  # single backbone
  comb.shift1 <- comb.shift[sapply(strsplit(comb.shift, "/"), length) == 1]
  
  # multiple backbone
  comb.shift2 <- comb.shift[sapply(strsplit(comb.shift, "/"), length) == 2]
  
  ####_____________ ####
  #### WHOLE TREE ####
  
  message("\n Estimating shifts of diversification from phylogeny.\n")
  message("-----------------------------------------------------------")
  if(!is.null(n.max)){
    message("\nA constrain will be applied: maximum diversity value is fixed at ", n.max, "\n")
  }
  
  if(!is.null(rate.max)){
    message("\nA constrain will be applied: maximum rate value is fixed at ", rate.max, "\n")
  }
  
  f1 <- Ntip(phylo)/totalsp
  message("\n--- WHOLE TREE ---\n \n")
  
  message("\n","Sampling fraction = ",  paste0(Ntip(phylo), "/", nrow(data), " (",round(f1,3)*100,"%)"), "\n") 
  
  res_phylo <- div.models(phylo = phylo, tot_time = max(branching.times(phylo)), f = f1,
                          cond = "crown", models = models, n.max = n.max, rate.max = rate.max)
  
  res_phylo[,-1] <- apply(res_phylo[,-1], 2, as.numeric)
  
  for(i in 1:nrow(res_phylo)){res_phylo$Parameters[i] <- 4-sum(is.na(res_phylo[i,]))}
  
  res_phylo <- res_phylo[,c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta")]
  all_res[[1]] <- res_phylo
  names(all_res)[[1]] <- "whole_tree"
  
  # all nodes to be tested
  all_tested_nodes <- unique(unlist(strsplit(gsub("/", "", comb.shift1), "[.]")))
  
  ALL_clade_names <- rep(list(NULL), length(unique(unlist(all_tested_nodes))))
  
  for(pot_names in 1:length(ALL_clade_names)){
    ALL_clade_names[pot_names] <- list(phylo$tip.label[unlist(Descendants(phylo, as.numeric(all_tested_nodes[pot_names])))])
  }
  names(ALL_clade_names) <- all_tested_nodes
  ALL_nodes_ages <- as.data.frame(apply(data.frame(nodesID=names(branching.times(phylo)),ages=branching.times(phylo)), 2, as.numeric))
  
  #### SUBCLADES ####
  message("\n--- SUBCLADES ---\n \n")
  ## Subclades trees
  # Check whether all groups are monophyletic 
  subclade_check <- NULL
  for(subclade in 1:length(ALL_clade_names)){
    if(!is.monophyletic(phy = phylo, tips = ALL_clade_names[[subclade]])){
      subclade_check <- c(subclade_check, names(ALL_clade_names[subclade]))
    }
  }
  if(!is.null(subclade_check)){
    message("The following subclades are not monophyletic:")
    message("",subclade_check,sep="\n  - ")
    stop()
  }
  
  phylo2 <- list()
  for(cla in 1:length(ALL_clade_names)){
    phylo2[[cla]] <- subtree(phylo, ALL_clade_names[[names(ALL_clade_names)[cla]]])
  }
  names(phylo2) <- names(ALL_clade_names)
  
  node_order <- names(branching.times(phylo)[order(branching.times(phylo))])
  node_order <- node_order[node_order %in% names(ALL_clade_names)]
  ALL_clade_names2 <- ALL_clade_names[match(node_order, names(ALL_clade_names))]
  
  #### ??? ####
  #for(cla in 1:length(ALL_clade_names2)){
    
  #  if(is.na(sampling.fractions$sp_tt[as.numeric(names(ALL_clade_names2)[cla])])){
  #    ancest_cla <- Ancestors(phylo, as.numeric(names(ALL_clade_names2)[cla]))
  #    cla_no_NA <- sampling.fractions$nodes[sampling.fractions$nodes %in% ancest_cla &
  #                                            !is.na(sampling.fractions$sp_in)]
  #    ancest_cla <- ancest_cla[ancest_cla %in% cla_no_NA]
  #    ancest_cla <- ancest_cla[ancest_cla %in% as.numeric(names(ALL_clade_names2))]
      
  #    ancest_cla <- ancest_cla[1]
      
  #    sampling.fractions$sp_tt[as.numeric(names(ALL_clade_names2)[cla])] <- sampling.fractions$sp_tt[ancest_cla]
  #    sampling.fractions$sp_in[as.numeric(names(ALL_clade_names2)[cla])] <- sampling.fractions$sp_in[ancest_cla]
  #  }
  #}
  
  totalsp2 <- as.list(sampling.fractions$sp_tt[as.numeric(names(ALL_clade_names))])
  phylosp2 <- as.list(sampling.fractions$sp_in[as.numeric(names(ALL_clade_names))])
  names(totalsp2) <- sampling.fractions$nodes[as.numeric(names(ALL_clade_names))]
  names(phylosp2) <- sampling.fractions$nodes[as.numeric(names(ALL_clade_names))]
  
  f2 <- as.list(unlist(phylosp2)/unlist(totalsp2))
  # totalsp2 <- totalsp2[names(totalsp2) %in% "others" == FALSE]
  
  final2<-list()
  # all parameters
  
  backbone <- rep(list(F),length(phylo2))
  spec_times <- unlist(list(rep(list(NULL),length(phylo2))),recursive = FALSE)
  branch_times <- unlist(list(rep(list(NULL),length(phylo2))),recursive = FALSE)
  tot_time2 <- unlist(list(rep(list(NULL),length(phylo2))),recursive = FALSE)
  cond <- unlist(list(rep(list(NULL),length(phylo2))),recursive = FALSE)
  
  # backbone selection
  if(backbone.option == "stem.shift"){
    cond <- rep(list("stem"),length(phylo2))
    for (i in 1:length(phylo2)){
      tot_time2[[i]] <- max(node.age(phylo2[[i]])$ages) + phylo2[[i]]$root.edge
    }
  }
  if(backbone.option == "crown.shift"){
    cond <- rep(list("crown"),length(phylo2))
    for (i in 1:length(phylo2)){
      tot_time2[[i]] <- max(node.age(phylo2[[i]])$ages)
    }
  }
  
  names(tot_time2) <- names(phylo2)
  
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
  
  names_phylo2 <- sampling.fractions[sampling.fractions$nodes %in% as.numeric(names(phylo2)),]
  names_phylo2 <- names_phylo2$data[match(as.numeric(names(phylo2)), names_phylo2$nodes)]
  names_phylo2 <- ifelse(is.na(names_phylo2), paste("at node", names(phylo2)), paste("for", names_phylo2))
  
  # subclades to parallelize
  #cl <- makeCluster(Ncores, type="SOCK")
  #clusterExport(cl, list("phylo2", "tot_time2", "backbone", "spec_times", "branch_times", "f2", "models", "div.models",
  #                       "fitLikelihood", "suppressWarnings", "try","getLikelihood", ".Phi", "integrate", ".Psi", ".Integrate"), envir = env.func)
  
  final2 <- lapply(seq_along(phylo2), function(i){
    message(i, "/", length(phylo2))
    n_tree <- sampling.fractions$sp_in[sampling.fractions$nodes == as.numeric(names(phylo2)[i])]
    n_tot <- sampling.fractions$sp_tt[sampling.fractions$nodes == as.numeric(names(phylo2)[i])]
    message("\n","Sampling fraction of ", names_phylo2[i], " = ",  paste0(n_tree, "/", n_tot, " (",round(f2[[i]],3)*100," %)"), "\n") 
    
    results <- div.models(phylo = phylo2[[i]], tot_time = tot_time2[[i]], f = f2[[i]],
                          cond = FALSE, models = models.sub, n.max = n.max, rate.max = rate.max)
    results1 <- div.models(phylo2[[i]], tot_time2[[i]], f = f2[[i]],
                           cond = cond[[i]], models = models.sub, n.max = n.max, rate.max = rate.max, verbose = FALSE)
    
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
  best_subclades_df <- data.frame(matrix(unlist(best_subclades), nrow = length(best_subclades), byrow = TRUE))
  names(best_subclades_df) <- names(best_subclades[[1]])
  best_subclades_df$Clades <- names(best_subclades)
  best_subclades_df <- best_subclades_df[,c("Clades", names(best_subclades[[1]]))]
  best_subclades_df <- best_subclades_df[,-length(best_subclades_df)]
  best_subclades_df[,-c(1,2)] <- apply(best_subclades_df[,-c(1,2)], 2, as.numeric)
  
  ####_____________ ####
  #### BACKBONES ####
  message("\n--- BACKBONES ---\n")
  
  # multi.backbone change
  if(multi.backbone == FALSE){
    # should be comb.shift1
    message("\n", "The ", length(comb.shift1)," combinations with simple backbones will be compared.\n")
    comb.shift <- comb.shift1
  }
  
  if(multi.backbone == TRUE){
    # should be comb.shift2
    message("\n", "The ", length(comb.shift2)," combinations with multiple backbones will be compared.\n")
    comb.shift <- comb.shift2
  }
  
  if(multi.backbone == "all"){
    # should be comb.shift
    message("\n", "All the ", length(comb.shift), " combinations will be compared.\n")
  }
  
  ALL_backbones <- rep(list(NULL),length(comb.shift))
  best_backbones <- NULL
  ALL_final3 <- rep(list(NULL),length(comb.shift))
  
  ALL_TOTAL <- NULL
  clades <- names(totalsp2)
  
  # Individual branching times ####
  # calculated for crown.shift option and then transform for stem.shift
  ALL_branch_times_clades <- rep(list(NULL),length(all_tested_nodes))
  names(ALL_branch_times_clades) <- all_tested_nodes
  
  for(clade in 1:length(all_tested_nodes)){
    
    parental_node <- Ancestors(phylo, as.numeric(all_tested_nodes[clade]), type = "parent")
    branch_times_clade <- unlist(list(rep(list(NULL),1)),recursive = FALSE)
    
    bt_cl <- as.numeric(c(all_tested_nodes[clade], parental_node))
    branch_times_clade[1] <- list(bt_cl)
    
    ALL_branch_times_clades[[clade]] <- branch_times_clade
  }
  
  # LOOP ON COMBINATIONS ####
  # Parallelization of backbone models
  message("\nDiversification models are running: \n")
  
  cl <- parallel::makeCluster(Ncores, type="SOCK")
  
  clusterExport(cl, list("all_comb_models","comb.shift", "subtree", "multi.backbone", "ALL_branch_times_clades", "phylo","Descendants", "Ancestors", "Siblings", "getMRCA",
                         "get.node.ages", "drop.tip", "all_tested_nodes", "ALL_backbones", "node.age", "totalsp", "totalsp2", "Ntip", 
                         "div.models", "fit_bd_backbone","fit_bd_backbone_c", "likelihood_bd_backbone", "n.max", ".Phi", "integrate", ".Psi", ".Integrate",
                         "branching.times", "ALL_clade_names", "sampling.fractions", "backbone.option", "models", "ALL_final3", "rate.max",
                         "Children", "extract.clade.ln", "expand.grid", "get.branching.nodes"), envir = env.func)
  
  ALL_final3 <- ParallelLogger::clusterApply(cl, seq_along(comb.shift), all_comb_models, progressBar = TRUE, stopOnError = TRUE)
  stopCluster(cl)
  
  #ALL_final3 <- ParallelLogger::clusterApply(cl, 1, all_comb_models, progressBar = TRUE, stopOnError = TRUE)
  
  #ALL_final3 <- lapply(1, all_comb_models)
  
  names(ALL_final3) <- comb.shift
  #ALL_final3 <- lapply(ALL_final3, unlist, recursive = FALSE)
  
  message("\n\n--- Comparison(s) of the ", length(ALL_final3), " combinations ---\n\n")
  
  best_ALL_final3 <- ALL_final3
  
  for(to in 1:length(ALL_final3)){
    for(sub_to in 1:length(ALL_final3[[to]])){
      ALL_final3[[to]][[sub_to]]$delta_AICc <- ALL_final3[[to]][[sub_to]]$AICc - min(ALL_final3[[to]][[sub_to]]$AICc)
      ALL_final3[[to]][[sub_to]] <- ALL_final3[[to]][[sub_to]][order(ALL_final3[[to]][[sub_to]]$AICc),]
      ALL_final3[[to]][[sub_to]][,-1] <- apply(ALL_final3[[to]][[sub_to]][,-1], 2, as.numeric)
      best_ALL_final3[[to]][[sub_to]] <- ALL_final3[[to]][[sub_to]][1,]
    }
  }
  
  # __ SELECTION #####
  ALL_TOTAL <- as.data.frame(matrix(ncol = 4, nrow = length(comb.shift)))
  names(ALL_TOTAL) <- c("Combination", "Parameters","logL","AICc")
  
  for(to in 1:length(best_ALL_final3)){
    names_subclades <- unlist(strsplit(sapply(strsplit(names(best_ALL_final3)[[to]], "/"), "[[", 1), split = "[.]"))
    
    for(sub_to in 1:length(best_ALL_final3[[to]])){
      
      ALL_TOTAL[to,-1] <- apply(do.call(rbind.data.frame, lapply(best_ALL_final3[[to]], function(x) x[, c("Parameters", "logL", "AICc")])), 2, sum)
      
      for(n_sub in 1:length(names_subclades)){
        ALL_TOTAL[to,-1] <- ALL_TOTAL[to,-1] + best_subclades_df[best_subclades_df$Clades %in% names_subclades[n_sub], c("Parameters","logL","AICc")]
      }
    }
  }
  
  ALL_TOTAL$Combination <- names(best_ALL_final3)
  
  res_phylo1 <- cbind("whole_tree", res_phylo[res_phylo$AICc == min(res_phylo$AICc), c("Parameters","logL","AICc")])
  names(res_phylo1) <- names(ALL_TOTAL)
  
  ALL_TOTAL <- rbind(ALL_TOTAL, res_phylo1)
  
  ALL_TOTAL$delta_AICc <- ALL_TOTAL$AICc - min(ALL_TOTAL$AICc)
  ALL_TOTAL <- ALL_TOTAL[order(ALL_TOTAL$delta_AICc),]
  row.names(ALL_TOTAL) <- NULL
  best_ALL_TOTAL <- ALL_TOTAL[ALL_TOTAL$delta_AICc < 2,]
  
  all_res[3] <- list(ALL_final3)
  all_res[4] <- list(ALL_TOTAL)
  names(all_res)[c(3,4)] <- c("backbones", "total")
  
  if(!"whole_tree" %in% best_ALL_TOTAL$Combination){
    message("\n A total of ", nrow(best_ALL_TOTAL), " combination(s) got the best fit(s) (delta AICc < 2).\n")
  } else{
    if(best_ALL_TOTAL$Combination[1] == "whole_tree"){
      message("\n No shift has been detected.\n")
    } else {
      if(which(best_ALL_TOTAL$Combination == "whole_tree")-1 == 1){
        message("\n",which(best_ALL_TOTAL$Combination == "whole_tree")-1," combination is better than the homogeneous model (but non-significant).\n")
      } else {
        message("\n",which(best_ALL_TOTAL$Combination == "whole_tree")-1,"combinations are better than the homogeneous model (but non-significant).\n")
      }
    }
  }
  return(all_res)
}

