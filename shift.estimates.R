# This function tests all combinations of subclades for which
# sampling fraction is available thanks to the taxonomy.
# 
# Update from April 03, 2020
#
#
# Update from August 11, 2020
# Ncores argument allows to specifies number of cores used for parallelization
# 
# Upgrade from November 09, 2020
# Version 0.1.0
# - Multiple backbone has been successfully implemented
# - MODELc constrains diversification rate estimates to a limit maximum set by n.max
# 
# Upgrade from February 10, 2021
# Version 0.1.1
# - Coded to be used with hypothesis testing.
#
# Upgrade from April 07, 2021
# Version 0.1.2
# - Better selection of initial values for parameters
# - new argument np.sub to apply a reduced set of models to subclades.

# Upgrade from April 30, 2021
# Version 0.1.3
# - new argument rate.max to constrain maximum rate value (cannot be used with n.max)

shift.estimates <- function(phy, data, sampling.fractions, comb.shift,
                            models = c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR"),
                            np.sub = 4, multi.backbone = F, backbone.option = "backbone2",
                            n.max = NULL, rate.max = NULL, 
                            Ncores = 1){
  env.func <- environment()
  options(echo = T)
  #### Loading packages ####
  source <- getwd()
  pkgs <- c("ape", "phytools", "phangorn","picante","stats", "parallel", "snow", "geiger", "ParallelLogger")
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new.pkgs)){install.packages(new.pkgs)}
  lapply(pkgs, require, character.only = T)
  
  # RPANDA codes and functions ####
  # to check why these functions are not in RPANDA package
  
  # additional functions
  #source("./tools/functions.for.shift.estimates.R")
  #source("./get.ncomb.shift.R")
  
  setwd(source)
  #### argument check ####
  if(!inherits(data, "data.frame")){
    stop("object \"data\" is not of class \"data.frame\"")
  }
  
  if(!inherits(phy, "phylo")){
    stop("object \"phy\" is not of class \"phylo\"")
  } else {
    phy$node.label <- c(c(Ntip(phy)+1):c(Ntip(phy)+Nnode(phy)))
  }
  
  if(!inherits(sampling.fractions, "data.frame")){
    stop("object \"sampling.fractions\" is not of class \"data.frame\"")
  } else {
    if(any(names(sampling.fractions) == "taxo")){
      names(sampling.fractions)[names(sampling.fractions) == "taxo"] <- "data"
    }
  }
  
  if(any("Species" %in% colnames(data)) == F){
    stop("No column named \"Species\" in the database.
         \nPlease rename the corresponding column with the name \"Species\".")
  }
  
  if(is(ncomb.shift)[1] != "list"){
    stop("Object \"ncomb.shift\" is incorrect.")
  } else {
    ALL_bck_comb <- ncomb.shift
  }
  
  if(any(!phy$tip.label %in% data$Species)){
    cat("The following tips are not in the database.\n \n")
    cat(phy$tip.label[!phy$tip.label %in% data$Species], "\n")
    stop()
  }
  
  if(all(sapply(ncomb.shift, is.null)) & multi.backbone == T){
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
  
  # Final list to  to return
  phylo <- ladderize(phy, F) # Mandatory to match node.labels
  
  all_res <- rep(list(NULL),4)
  
  totalsp <- nrow(data) # number of species in datanomic database
  data_phylo <- data[data$Species %in% phy$tip.label,]
  
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
  
  cat("\n","Sampling fraction =",  paste0(Ntip(phy), "/", nrow(data), " (",round(f1,3)*100," %)"), "\n") 
  
  res_phy <- div.models_RPANDA(phylo = phy, tot_time = max(node.age(phylo)$ages), f = f1, 
                               cond = "crown", models = models, n.max = n.max, rate.max = rate.max)
  
  res_phy[,-1] <- apply(res_phy[,-1], 2, as.numeric)
  
  for(i in 1:nrow(res_phy)){res_phy$Parameters[i] <- 4-sum(is.na(res_phy[i,]))}
  
  res_phy <- res_phy[,c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta")]
  
  all_res[[1]] <- res_phy ; names(all_res)[[1]] <- "whole_tree"
  
  all_tested_nodes <- unique(c(unique(unlist(strsplit(names(ncomb.shift), "[.]"))), unique(unlist(ncomb.shift))))
  
  all_lineages <- unique(unlist(strsplit(names(ALL_bck_comb), "[.]")))
  
  ALL_clade_names <- rep(list(NULL), length(unique(unlist(all_tested_nodes))))
  
  for(pot_names in 1:length(ALL_clade_names)){
    ALL_clade_names[pot_names] <- list(phy$tip.label[unlist(Descendants(phy, as.numeric(all_tested_nodes[pot_names])))])
  }
  names(ALL_clade_names) <- all_tested_nodes
  
  ALL_nodes_ages <- as.data.frame(apply(data.frame(nodesID=names(branching.times(phy)),ages=branching.times(phy)), 2, as.numeric))
  
  ALL_clade_names1 <- ALL_clade_names[all_lineages]
  
  #### SUBCLADES ####
  cat("\n--- SUBCLADES ---\n \n")
  ## Subclades trees
  # Check whether all groups are monophyletic 
  subclade_check <- NULL
  for(subclade in 1:length(ALL_clade_names1)){
    if(!is.monophyletic(phy = phy, tips = ALL_clade_names1[[subclade]])){
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
    phylo2[[cla]] <- subtree(phy, unlist(ALL_clade_names1[names(ALL_clade_names1)[cla]]))
  }
  
  names(phylo2) <- names(ALL_clade_names1)
  
  node_order <- names(branching.times(phy)[order(branching.times(phy))])
  node_order <- node_order[node_order %in% names(ALL_clade_names1)]
  ALL_clade_names2 <- ALL_clade_names1[match(node_order, names(ALL_clade_names1))]
  
  for(cla in 1:length(ALL_clade_names2)){
    
    if(is.na(sampling.fractions$sp_tt[as.numeric(names(ALL_clade_names2)[cla])])){
      ancest_cla <- Ancestors(phy, as.numeric(names(ALL_clade_names2)[cla]))
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
  
  #cl <- makeCluster(Ncores, type="SOCK")
  #clusterExport(cl, list("phylo2", "tot_time2", "backbone", "spec_times", "branch_times", "f2", "models", "div.models",
  #                       "fitLikelihood", "suppressWarnings", "try","getLikelihood", ".Phi", "integrate", ".Psi", ".Integrate"), envir = env.func)
  
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
  
  final2 <- lapply(seq_along(phylo2), function(i){
    cat(i, "/", length(phylo2))
    n_tree <- sampling.fractions$sp_in[sampling.fractions$nodes == as.numeric(names(phylo2)[i])]
    n_tot <- sampling.fractions$sp_tt[sampling.fractions$nodes == as.numeric(names(phylo2)[i])]
    cat("\n","Sampling fraction", names_phylo2[i], "=",  paste0(n_tree, "/", n_tot, " (",round(f2[[i]],3)*100," %)"), "\n") 
    
    results <- div.models_RPANDA(phylo2[[i]], tot_time2[[i]], f = f2[[i]],
                                 cond = F, models = models.sub, n.max = n.max, rate.max = rate.max)
    results1 <- div.models_RPANDA(phylo2[[i]], tot_time2[[i]], f = f2[[i]],
                                  cond = cond[[i]], models = models.sub, n.max = n.max, rate.max = rate.max, verbose = F)
    
  
    results <- merge(results[, c("Models", "Parameters", "logL", "AICc")], results1[,c("Models", "Lambda", "Alpha", "Mu", "Beta")], by = "Models")
    
    results[,-1] <- apply(results[,-1], 2, as.numeric)
    
    results$delta_AICc <- results$AICc - min(results$AICc)
    results <- results[order(results$delta_AICc),]
    
    #final2[[i]] <- results
    return(results)
    
  })
  #, progressBar = T, stopOnError = T
  #stopCluster(cl)
  
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
    
    parental_node <- Ancestors(phy, as.numeric(all_tested_nodes[clade]), type = "parent")
    if(parental_node == Ntip(phy)+1){ # this clade is at the root
      root_clade <- 1
    } else {
      root_clade <- 0
    }
    branch_times_clade <- unlist(list(rep(list(NULL),1 + root_clade)),recursive = F)
    
    bt_cl <- as.numeric(c(all_tested_nodes[clade], parental_node))
    # if(is.na(bt_cl[2])){}
    branch_times_clade[1] <- list(bt_cl)
    
    if(root_clade == 1){
      branch_times_clade[1 + root_clade] <- list(c(Siblings(phy, max(bt_cl)),min(bt_cl)))
    }  
    ALL_branch_times_clades[[clade]] <- branch_times_clade
  }
  
  # LOOP ON COMBINATIONS ####
  
  # Parallelization of backbone models
  cat("\nDiversification models are running: \n")
  
  cl <- parallel::makeCluster(Ncores, type="SOCK")
  # Check what comes form RPANDA
  clusterExport(cl, list("all_comb_models", "subtree", "multi.backbone", "ALL_branch_times_clades", "phy", "ALL_bck_comb","Descendants", "Ancestors", "Siblings", "getMRCA",
                         "get.node.ages", "drop.tip", "all_tested_nodes", "ALL_backbones", "node.age", "totalsp", "totalsp2", "Ntip", 
                         "div.models_RPANDA", "fit_bd_backbone", "likelihood_bd_backbone", "n.max", ".Phi", "integrate", ".Psi", ".Integrate",
                         "branching.times", "ALL_clade_names", "sampling.fractions", "backbone.option", "models", "ALL_final3", "rate.max",
                         "Children", "extract.clade.ln", "expand.grid", "get.branching.nodes"), envir = env.func)
  
  ALL_final3 <- ParallelLogger::clusterApply(cl, seq_along(ALL_bck_comb), all_comb_models, progressBar = T, stopOnError = T)
  
  # for debugging
  # ALL_final3 <- lapply(1, all_comb_models)
  
  stopCluster(cl)
  
  sapply(ALL_final3, function(x) x[[1]][[1]]$Lambda < 0)
  
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
  
  #ALL_TOTAL_df <- ALL_TOTAL_df[-ALL_TOTAL_df$AICc < 0,]
  res_phy1 <- cbind("whole_tree", res_phy[res_phy$AICc == min(res_phy$AICc), c("Parameters","logL","AICc")])
  names(res_phy1) <- names(ALL_TOTAL)
  
  ALL_TOTAL <- rbind(ALL_TOTAL, res_phy1)
  
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
  
  setwd(source)
  
  cat("\n A total of", nrow(best_ALL_TOTAL), "combination(s) got the best fit(s) (delta AICc < 2). \n")
  
  
  return(all_res)
  
}