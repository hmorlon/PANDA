get.comb.shift <- function(phylo, data, sampling.fractions, clade.size = 5, Ncores = 1){
  
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
  }
  
  phylo$node.label <- c(Ntip(phylo) + 1):c(Ntip(phylo) + Nnode(phylo))
  
  if(any("Species" %in% colnames(data)) == FALSE){
    stop("No column named \"Species\" in the database
         \nPlease rename the corresponding column with the name \"Species\".")
  }
  
  if(any(!phylo$tip.label %in% data$Species)){
    message("The following tips are not in database.\n \n")
    message(phylo$tip.label[!phylo$tip.label %in% data$Species], "\n")
    stop()
  }
  
  if(phylo$Nnode + Ntip(phylo) != nrow(sampling.fractions) | is(sampling.fractions)[1]!="data.frame"){
    stop("object \"sampling.fractions is not of class \"data.frame\" or it does not correspond to the provided phylogeny")
  }

    # START
  #### SUBCLADES ####
  # Steps #
  # 1) list all potential subclades (named by nodes)

  potential_clades <- sampling.fractions[!is.na(sampling.fractions$to_test),]
  potential_clades$data <- ifelse(is.na(potential_clades$data), paste("f", sampling.fractions[!is.na(sampling.fractions$to_test), c("nodes")], sep="."), potential_clades$data)
  
  # 2) create coalescence of subclades to exclude combinaisons that are not possible
  
  coal_potential_clades <- Ancestors(phylo, potential_clades$nodes, "all")
  #coal_potential_clades <- lapply(coal_potential_clades, function(x) x[!x %in% 1:Ntip(phylo)])
  #coal_potential_clades <- lapply(coal_potential_clades, rev)
  max_length_coal <- max(sapply(coal_potential_clades, length))
  df_coal <- as.data.frame(matrix(nrow = length(coal_potential_clades), ncol = max_length_coal+1))
  for(n_coal in 1:length(coal_potential_clades)){
    df_coal[n_coal,] <- c(rep(NA, max_length_coal - length(coal_potential_clades[[n_coal]])), potential_clades$nodes[n_coal], coal_potential_clades[[n_coal]])
  }
  
  to_remove <- NULL
  df_coal_pat <- apply(df_coal, 1, function(x) paste(x[!is.na(x)], collapse = "."))
  for(coal_p in 1:length(df_coal_pat)){
    pattern <- df_coal_pat[coal_p]
    df_coal_pat_not_pattern <- df_coal_pat[df_coal_pat != pattern]
    if(any(grepl(pattern, df_coal_pat_not_pattern))){
      to_remove <- c(to_remove, pattern)
    }
  }
  
  diff_lineages <- df_coal_pat[!df_coal_pat %in% to_remove]
  diff_lineages <- strsplit(diff_lineages,"[.]")
  diff_lineages <- lapply(diff_lineages, function(x) x[!x %in% as.character(Ntip(phylo)+1)]) # remove the root
  diff_lineages <- lapply(diff_lineages, function(x) x[!x %in% sampling.fractions$nodes[is.na(sampling.fractions$f)]]) # remove node without f
  
  all_lineages <- unique(unlist(diff_lineages))
  # 3) for each group of nested subclades, find the best one (maybe test from present to past)
  
  #cat(length(diff_lineages), "lineages detected")
  #cat("\n--- Calculating all combinations of potential shifts from 1 to", length(diff_lineages),"---\n")
  
  n_clade_comb <- NULL
  n_c = 0
  for(x in 1:length(diff_lineages)){
    n_c = n_c + 1
    n_clade_comb[n_c] <- list(t(combn(1:length(diff_lineages),x)))
  }
  
  names(n_clade_comb) <- 1:length(diff_lineages)
  get.all.comb <- function(n){
    ALL_comb <- NULL
    poor_nodes <- sapply(unlist(diff_lineages), function(x) length(unlist(Descendants(phylo, as.numeric(x), type = "tips"))))
    poor_nodes <- names(poor_nodes[poor_nodes < clade.size])
    
    for(n_row in 1:nrow(n_clade_comb[[n]])){
      diff_lineages_n <- diff_lineages[n_clade_comb[[n]][n_row,]]
      shared_nodes <- sapply(unlist(diff_lineages_n), function(x) length(unlist(diff_lineages_n)[unlist(diff_lineages_n) %in% x]))
      shared_nodes <- names(shared_nodes[shared_nodes >= 2 ])
      diff_lineages_n <- lapply(diff_lineages_n, function(x) x[!x %in% c(shared_nodes, poor_nodes)])
      ALL_comb <- c(ALL_comb, apply(expand.grid(diff_lineages_n), 1, function(x) paste(x, collapse = ".")))
    }
    ALL_comb <- unique(ALL_comb)
    print("next")
    return(ALL_comb)
    
  }
  message("\n SIMPLE BACKBONES:\n")
  cl <- parallel::makeCluster(Ncores, type="SOCK")
  clusterExport(cl, list("phylo","Descendants", "diff_lineages", "Siblings","n_clade_comb","expand.grid","paste", "clade.size"),
                envir = env.func)
  ALL_comb <- unlist(ParallelLogger::clusterApply(cl, seq_along(n_clade_comb), get.all.comb, progressBar = TRUE))
  stopCluster(cl)
  
  names(ALL_comb) <- NULL
  ALL_comb <- strsplit(ALL_comb, split = "[.]")
  
  # Because of too poor backbones
  comb_to_remove <- NULL
  for(comb in 1:length(ALL_comb)){
    if(c(Ntip(phylo) - length(phylo$tip.label[unlist(Descendants(phylo, as.numeric(ALL_comb[[comb]])))])) < clade.size){
      comb_to_remove <- c(comb_to_remove, comb)
    }
  }
  
  ALL_comb <- ALL_comb[!1:length(ALL_comb) %in% comb_to_remove]
  names(ALL_comb) <- NULL
  
  ALL_clade_names <- rep(list(NULL), length(unique(unlist(all_lineages))))
  
  for(pot_names in 1:length(ALL_clade_names)){
    ALL_clade_names[pot_names] <- list(phylo$tip.label[unlist(Descendants(phylo, as.numeric(all_lineages[pot_names])))])
  }
  names(ALL_clade_names) <- unlist(all_lineages)
  
  # MULTI-BACKBONE ####
  # diff_lineages is crucial because we need to cut a lineage !!
  ALL_bck_comb <- rep(list(NULL),length(ALL_comb))
  names(ALL_bck_comb) <- lapply(ALL_comb, function(x) paste(x, collapse = "."))
  
  message("\n MULTIPLE BACKBONES:\n")
  cl <- parallel::makeCluster(Ncores, type="SOCK")
  clusterExport(cl, list("phylo","Descendants", "Ancestors", "diff_lineages", "Siblings","n_clade_comb","expand.grid","paste",
                         "ALL_comb", "ALL_clade_names", "sampling.fractions","Ntip", "clade.size", "ALL_bck_comb", "drop.tip",
                         "branching.times"),
                envir = env.func)
  
  ALL_bck_comb <- ParallelLogger::clusterApply(cl, seq_along(ALL_comb), function(comb){
    
    message(comb, "/",length(ALL_comb), "\n")
    names(ALL_bck_comb)[comb] <- paste(ALL_comb[[comb]], collapse = ".")
    
    phylo_bck <- drop.tip(phylo, unlist(ALL_clade_names[ALL_comb[[comb]]]))
    
    sub_node <- c(Ancestors(phylo, as.numeric(unlist(ALL_comb[comb])), "parent"),
                  as.numeric(unlist(ALL_comb[comb])),
                  unlist(Descendants(phylo, as.numeric(unlist(ALL_comb[comb])), "all")))
    
    sf_bck <- sampling.fractions[!sampling.fractions$nodes %in% sub_node,c("nodes", "data", "f","to_test")]
    matching <- data.frame(bck.nodes = unique(phylo_bck$edge[,1]), phylo = phylo_bck$node.label)
    sf_bck <- merge(sf_bck, matching, by.x = "nodes", by.y = "phylo", all = TRUE)
    
    sf_bck <- sf_bck[!is.na(sf_bck$to_test),]
    sub_diff_lineages <- diff_lineages[sapply(diff_lineages, function(x) any(ALL_comb[[comb]] %in% x))]
    
    #sf_bck <- sf_bck[!sf_bck$nodes %in% ALL_comb[[comb]]]
    
    # lineage to cut
    same_lineage <- rep(list(NULL),length(ALL_comb[[comb]]))
    
    for(node.c in 1:length(ALL_comb[[comb]])){
      lin <- sub_diff_lineages[sapply(sub_diff_lineages, function(x) ALL_comb[[comb]][node.c] %in% x)]
      
      if(length(lin) > 1){
        for(slin in 1:length(lin)){
          lin[[slin]] <- lin[[slin]][c(which(lin[[slin]] %in% ALL_comb[[comb]][node.c]):length(lin[[slin]]))]
        }
      }
      
      lin <- unique(lin)[[1]]
      
      other_lin <- sub_diff_lineages
      other_lin[sapply(sub_diff_lineages, function(x) ALL_comb[[comb]][node.c] %in% x)] <- list(NULL)
      par_nodes_sub <- c(ALL_comb[[comb]], Ancestors(phylo,as.numeric(ALL_comb[[comb]]),"parent"))
      
      lin <- lin[!lin %in% par_nodes_sub]
      if(length(lin) != 0){
        
        lin.node <- data.frame(node = lin, n.tips = rep(NA, length(lin)))
        for(n.lin in 1:length(lin)){
          desc.n.lin <- length(Descendants(phylo, as.numeric(lin[n.lin]))[[1]])
          # whether this node is present in an other lineage
          if(any(sapply(other_lin, function(x) any(lin[n.lin] %in% x)))){
            sub_to_rem <- unique(c(node.c, which(sapply(other_lin, function(x) any(lin[n.lin] %in% x)))))
            desc.n.lin <- desc.n.lin - length(unlist(Descendants(phylo, as.numeric(ALL_comb[[comb]][sub_to_rem]), "tips")))
          }else{
            desc.n.lin <- desc.n.lin - length(Descendants(phylo, as.numeric(ALL_comb[[comb]][node.c], "tips"))[[1]])
          }
          
          if(desc.n.lin >= clade.size){
            lin.node$n.tips[n.lin] <- desc.n.lin
            same_lineage[[node.c]] <- c(same_lineage[[node.c]], lin[n.lin])
          }
        }
      }
    }
    # combinations of backbones should be done on all nodes independently from the lineage
    
    bck_nodes <- unique(unlist(same_lineage))
    if(is.null(bck_nodes)){
      #next()
      ALL_bck_comb[comb] <- list(NULL)
    } else {
      
      n_c = 0
      n_same_comb <- NULL
      
      for(cut in 1:length(bck_nodes)){
        n_c = n_c + 1
        n_same_comb[n_c] <- list(t(combn(1:length(bck_nodes),cut)))
      }
      
      ALL_bck_comb_n <- NULL
      for(n_cut in 1:length(n_same_comb)){
        for(n_row in 1:nrow(n_same_comb[[n_cut]])){
          ALL_bck_comb_n <- c(ALL_bck_comb_n, list(bck_nodes[n_same_comb[[n_cut]][n_row,]]))
        }
      }
      
      lin.node <- data.frame(node = c(bck_nodes,Ntip(phylo)+1), n.tips = rep(NA, length(bck_nodes)+1))
      lin.node$node <- as.character(lin.node$node)
      for(n.lin in 1:nrow(lin.node)){
        desc.n.lin <- length(Descendants(phylo, as.numeric(lin.node$node[n.lin]))[[1]])
        # whether this node is present in an other lineage
        int.n.lin <- Descendants(phylo, as.numeric(lin.node$node[n.lin]), type = "all")
        int.n.lin <- as.character(int.n.lin[int.n.lin > Ntip(phylo)])
        sub_to_rem <- ALL_comb[[comb]] %in% int.n.lin
        if(any(ALL_comb[[comb]] %in% int.n.lin)){
          desc.n.lin <- desc.n.lin - length(unlist(Descendants(phylo, as.numeric(ALL_comb[[comb]][sub_to_rem]), "tips")))
        }
        lin.node$n.tips[n.lin] <- desc.n.lin
      }
      
      # Selecting comb to remove
      bck_comb_to_remove <- NULL
      
      for(n_bck_comb in 1:length(ALL_bck_comb_n)){
        sub.lin.node <- lin.node[lin.node$node %in% c(ALL_bck_comb_n[[n_bck_comb]], Ntip(phylo)+1),]
        
        sub.lin.node$n.tips_prev <- sub.lin.node$n.tips
        node_order <- names(branching.times(phylo)[order(branching.times(phylo))])
        node_order <- node_order[node_order %in% sub.lin.node$node]
        
        sub.lin.node <- sub.lin.node[match(sub.lin.node$node, node_order),]
        
        for(l.n in 1:nrow(sub.lin.node)){
          int.desc_lin <- unlist(Descendants(phylo, as.numeric(sub.lin.node$node[l.n]), "all"))
          int.desc_lin <- int.desc_lin[int.desc_lin > Ntip(phylo)]
          if(any(sub.lin.node$node %in% int.desc_lin)){
            
            bck_up <- sub.lin.node[which(sub.lin.node$node %in% int.desc_lin),]
            
            ntip_bck_up <- sum(bck_up$n.tips_prev)
            
            sub.lin.node$n.tips_prev[l.n] <- sub.lin.node$n.tips[l.n] - ntip_bck_up
            
          } else {
            sub.lin.node$n.tips_prev[l.n] <- sub.lin.node$n.tips[l.n]
          }
        }
        # from construction and for debugging #
        #if(sum(sub.lin.node$n.tips_prev) != sub.lin.node$n.tips[nrow(sub.lin.node)]){
        #  stop("Error in multi-backbone size calculation...")
        #}
        
        if(any(sub.lin.node$n.tips_prev < clade.size)){
          bck_comb_to_remove <- c(bck_comb_to_remove, n_bck_comb)
        }
      }
      
      if(!is.null(bck_comb_to_remove)){
        ALL_bck_comb_n <- ALL_bck_comb_n[-bck_comb_to_remove]
      }
      
      ALL_bck_comb[[comb]] <- ALL_bck_comb_n
      
      if(length(ALL_bck_comb_n) == 0){
        ALL_bck_comb[comb] <- list(NULL)
      }
      
    }
    return(ALL_bck_comb[[comb]])
  }, progressBar = TRUE)
  stopCluster(cl)
  
  names(ALL_bck_comb) <- sapply(ALL_comb, function(x) paste(x, collapse = "."))
  
  # to transform here
  
  ALL_bck_comb1 <- c()
  for(i in 1:length(ALL_bck_comb)){
    if(is.list(ALL_bck_comb[[i]])){
      ALL_bck_comb1 <- c(ALL_bck_comb1, paste0(names(ALL_bck_comb[i]), "/"))
      for(j in 1:length(ALL_bck_comb[[i]])){
        ALL_bck_comb1 <- c(ALL_bck_comb1, paste(names(ALL_bck_comb)[i], paste(ALL_bck_comb[[i]][[j]], collapse = "."), sep = "/"))
      }
    } else {
      ALL_bck_comb1 <- c(ALL_bck_comb1, paste(names(ALL_bck_comb)[i], ALL_bck_comb[[i]], sep = "/"))
    }
  }
  

  message("\n ONLY SIMPLE BACKBONES:")
  message("\t\n",   sum(sapply(strsplit(ALL_bck_comb1, "/"), length) == 1) ," combination(s) have been detected. \n")
  message("\n WITH MULTIPLE BACKBONES:")
  message("\n", length(ALL_bck_comb1), " combination(s) have been detected. \n")
  
  return(ALL_bck_comb1)
}
