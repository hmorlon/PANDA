simul_posteriors <- function(n = 10000, phylo, sampling.fractions,
                             shift.res, combi = 1, clade.size = 5){
  
  #### argument check ####
  if(!inherits(phylo, "phylo")){
    stop("object \"phylo\" is not of class \"phylo\"")
  } else {
    phylo$node.label <- c(c(Ntip(phylo)+1):c(Ntip(phylo)+Nnode(phylo)))
  }
  if(!inherits(sampling.fractions, "data.frame")){
    stop("object \"sampling.fractions\" is not of class \"data.frame\"")
  }
  if(!is(shift.res)[1] == "list" | any(names(shift.res) != c("whole_tree", "subclades", "backbones", "total"))){
    stop("object \"shift.res\" might be incorrect.")
  }
  if(!is.numeric(combi)){
    stop("argument \"combi\" should be numeric.")
  }
  if(!is.numeric(clade.size)){
    stop("argument \"clade.size\" should be numeric.")
  }
  if(!is.numeric(n)){
    stop("argument \"n\" should be numeric.")
  }
  
  param_equation <- function(model, param, tot_time){
    
    if(grepl("BCST", model)){
      speciation <- param$Lambda
      extinction <- 0
    }
    
    if(grepl("BCST_DCST", model)){
      speciation <- param$Lambda
      extinction <- param$Mu
    }
    
    if(grepl("BVAR", model)){
      t <- c(0:floor(tot_time), tot_time)
      speciation <- function(t) param$Lambda*exp(param$Alpha*tot_time)*exp(-param$Alpha*t)
      extinction <- 0
    }
    
    if(grepl("BVAR_DCST", model)){
      t <- c(0:floor(tot_time), tot_time)
      speciation <- function(t) param$Lambda*exp(param$Alpha*tot_time)*exp(-param$Alpha*t)
      extinction <- param$Mu
    }
    
    if(grepl("BCST_DVAR", model)){
      t <- c(0:floor(tot_time), tot_time)
      speciation <- param$Lambda
      extinction <- function(t) param$Mu*exp(param$Beta*tot_time)*exp(-param$Beta*t)
    }
    return(list(speciation = speciation, extinction = extinction))
  }
  
  # test posterior predictive with Cetacea
  comb <- shift.res$total$Combination[combi]
  if(length(grep("/", comb)) == 1){
    if(length(strsplit(comb, "/")[[1]]) > 1){
      comb.sub <- strsplit(sapply(strsplit(comb, "/"), "[[", 1), "[.]")[[1]]
      comb.bck <- strsplit(sapply(strsplit(comb, "/"), "[[", 2), "[.]")[[1]]
    } else{
      comb.sub <- strsplit(sapply(strsplit(comb, "/"), "[[", 1), "[.]")[[1]]
      comb.bck <- NULL
    }
  } else {
    comb.sub <- strsplit(sapply(strsplit(comb, "/"), "[[", 1), "[.]")[[1]]
    comb.bck <- NULL
  }
  ages <- branching.times(phylo)[comb.sub]
  comb.sub <- comb.sub[order(ages)]
  ages <- ages[order(ages)]
  
  if(length(comb.sub)+length(comb.bck)+1 > 26){
    let <- apply(expand.grid(letters[1], letters), 1, paste, collapse = "")
  }else{
    let <- letters[1:length(comb.sub)]
  }
  
  # 
  crown_age <- max(branching.times(phylo))
  ages_bck <- ifelse(!is.null(comb.bck), branching.times(phylo)[comb.bck], NA)
  comb.bck <- comb.bck[order(ages_bck)]
  ages_bck <- c(sort(ages_bck), crown_age)
  
  names(ages_bck) <- c(comb.bck, Ntip(phylo)+1)
  
  fs1 <- sampling.fractions[sampling.fractions$nodes %in% c(comb.sub, comb.bck), c("nodes", "f")]
  fs <- fs1$f
  names(fs) <- fs1$nodes
  
  # Sampling fractions backbone ####
  
  lin.node <- data.frame(node = c(comb.sub,comb.bck, Ntip(phylo)+1), n.tips = rep(NA, length(comb.sub) + length(comb.bck)+1))
  lin.node$node <- as.character(lin.node$node)
  lin.node <- merge(lin.node, sampling.fractions[sampling.fractions$nodes %in% lin.node$node, c("nodes", "sp_tt"),],
                    by.x = "node", by.y = "nodes")
  
  node_order <- names(branching.times(phylo)[order(branching.times(phylo))])
  node_order <- node_order[node_order %in% lin.node$node]
  
  lin.node <- lin.node[match(node_order, lin.node$node),]
  
  for(n.lin in 1:nrow(lin.node)){
    desc.n.lin <- length(Descendants(phylo, as.numeric(lin.node$node[n.lin]))[[1]])
    # whether this node is present in an other lineage
    int.n.lin <- Descendants(phylo, as.numeric(lin.node$node[n.lin]), type = "all")
    int.n.lin <- as.character(int.n.lin[int.n.lin > Ntip(phylo)])
    # Ntip
    if(any(comb.sub %in% int.n.lin)){
      lin.node$n.tips[n.lin] <- desc.n.lin - sum(lin.node$n.tips[lin.node$node %in% comb.sub[comb.sub %in% int.n.lin]])
      lin.node$sp_tt[n.lin] <- lin.node$sp_tt[n.lin] - sum(lin.node$sp_tt[lin.node$node %in% comb.sub[comb.sub %in% int.n.lin]])
    } else{
      lin.node$n.tips[n.lin] <- desc.n.lin
    }
  }
  
  lin.node$n.tips_prev <- lin.node$n.tips
  lin.node$sp_tt_prev <- lin.node$sp_tt
  
  lin.node_bck <- lin.node[!lin.node$node %in% comb.sub,]
  
  for(l.n in c(1:nrow(lin.node_bck))){
    int.desc_lin <- unlist(Descendants(phylo, as.numeric(lin.node_bck$node[l.n]), "all"))
    int.desc_lin <- int.desc_lin[int.desc_lin > Ntip(phylo)]
    
    if(any(lin.node_bck$node %in% int.desc_lin)){
      
      bck_up <- lin.node_bck[which(lin.node_bck$node %in% int.desc_lin),]
      
      ntip_bck_up <- sum(bck_up$n.tips_prev)
      ntaxo_bck_up <- sum(bck_up$sp_tt_prev)
      
      lin.node_bck$n.tips_prev[l.n] <- lin.node_bck$n.tips[l.n] - ntip_bck_up
      lin.node_bck$sp_tt_prev[l.n] <- lin.node_bck$sp_tt[l.n] - ntaxo_bck_up  
      
    }
  }
  lin.node[lin.node$node %in% lin.node_bck$node,] <- lin.node_bck
  
  lin.node <- lin.node[-(1:length(comb.sub)),]
  
  f <- as.list(lin.node$n.tips_prev/lin.node$sp_tt_prev)
  
  names(f) <- c(comb.bck, Ntip(phylo)+1)
  fs <- c(fs, unlist(f))

  fs <- fs[match(c(names(ages), names(ages_bck)), names(fs))]
  
  anc_sub <- Ancestors(phylo, as.numeric(comb.sub))
  comb.sub_bybck <-c()
  backbones <- rep(list(NULL),length(comb.bck)+1)
  cat("Backbone simulations...\n")
  for(cb in 1:length(backbones)){
    if(length(backbones) > 1){
      comb.sub_bybck[[cb]] <- comb.sub[sapply(anc_sub, function(x) comb.bck[cb] %in% x)]
      model_backbone <- shift.res$backbones[shift.res$total$Combination[combi]][[1]][[cb]]$Models[1]
      param_backbone <- shift.res$backbones[shift.res$total$Combination[combi]][[1]][[cb]][1, c("Lambda", "Alpha", "Mu", "Beta")]
      equation_backbone <- param_equation(model_backbone, param_backbone, ages_bck[cb])
      
      if(cb > 1){
        comb.sub_bybck[[cb]] <- c(comb.sub_bybck[[cb]], comb.bck[comb.bck %in% Descendants(phylo, as.numeric(comb.bck[cb]), "all")])
      }
      
      bck_cb <- tess.sim.age(n = n, lambda = equation_backbone$speciation,
                             mu = equation_backbone$extinction, age = ages_bck[cb])
      bck_cb <- bck_cb[sapply(bck_cb, Ntip) > clade.size+length(comb.sub_bybck[[cb]])*2] # five tips in the backbone (maybe more)
      backbones[[cb]] <- lapply(bck_cb, ladderize)
    }
    
    # deep backbone
    if(cb == length(backbones)){
      comb.sub_bybck[[cb]] <- c(comb.sub[!comb.sub %in% unlist(comb.sub_bybck[[cb]])], comb.bck)
      model_backbone <- shift.res$backbones[shift.res$total$Combination[combi]][[1]][[cb]]$Models[1]
      param_backbone <- shift.res$backbones[shift.res$total$Combination[combi]][[1]][[cb]][1, c("Lambda", "Alpha", "Mu", "Beta")]
      equation_backbone <- param_equation(model_backbone, param_backbone, ages_bck[cb])
      
      bck_cb <- tess.sim.age(n = n, lambda = equation_backbone$speciation,
                             mu = equation_backbone$extinction, age = ages_bck[cb])
      
      bck_cb <- bck_cb[sapply(bck_cb, Ntip) > clade.size+length(comb.sub_bybck[[cb]])*2]
      # five tips in the backbone + the one we will prune (maybe more)
      backbones[[cb]] <- lapply(bck_cb, ladderize)
    }
  }
  if(length(backbones) > 1){
    names(backbones) <- c(comb.bck, Ntip(phylo)+1)
  } else {
    names(backbones) <- as.character(c(Ntip(phylo)+1))
  }
  
  # subclades 
  # subclades
  cat("Subclade simulations...\n")
  all_subclades <- c()
  for(i in 1:length(comb.sub)){
    cat("\t subclade", i, "/", length(comb.sub),"\n")
    # param
    param_sub <- shift.res$subclades[comb.sub[i]][[1]][1, c("Lambda", "Alpha", "Mu", "Beta")]
    model_sub <- shift.res$subclades[comb.sub[i]][[1]]$Models[1]
    
    equation_subclade_i <- param_equation(model_sub, param_sub, ages[i])
    
    all_subclades[[i]] <- tess.sim.age(n = n, lambda = equation_subclade_i$speciation,
                                       age = ages[i], mu = equation_subclade_i$extinction,
                                       samplingProbability = fs[i])
    for(j in 1:length(all_subclades[[i]])){
      all_subclades[[i]][[j]]$tip.label <- gsub("t", let[i], all_subclades[[i]][[j]]$tip.label)
    }
    
    all_subclades[[i]] <- all_subclades[[i]][sapply(all_subclades[[i]], Ntip) > clade.size]
    
  }
  names(all_subclades) <- comb.sub
  
  
  all_parts <- c(all_subclades, backbones)
  all_parts <- lapply(1:min(sapply(all_parts, length)), function(j) lapply(all_parts, function(x) x[[j]]))
  
  all_comb.sub_bybck <- c(rep(list(NULL), length(all_subclades)), comb.sub_bybck)
  all_ages <- c(ages, ages_bck)
  
  # all shifts together
  
  all_shift_nodes <- c()
  for(j in 1:length(all_parts)){
    #cat(j, "/", length(all_parts), "\n")
    shift_nodes <- c()
    shift_nodes_and_anc <- c()
    for(i in 1:length(all_parts[[j]])){
      subtree_bck <- all_parts[[j]][[i]]
      if(i > length(all_subclades)){
        subtree_bck$tip.label <- gsub("t", letters[length(letters) - length(all_parts[[j]]) + i], subtree_bck$tip.label)
        all_parts[[j]][[i]] <- subtree_bck
      }
      if(i != length(all_parts[[j]])){
        subtree_bck_loc <- which(sapply(all_comb.sub_bybck, function(x) names(all_parts[[j]])[i] %in% x))
        tree <- all_parts[[j]][subtree_bck_loc][[1]]
        
        # selecting nodes
        
        edge_ages <- apply(tree$edge, 2, function(x) ifelse(x %in% 1:Ntip(tree), 0, branching.times(tree)[as.character(x)]))
        edge_selection <- tree$edge[apply(edge_ages, 1, function(x) all_ages[i] < x[1] & all_ages[i] > x[2]),]
        edge_selection <- edge_selection[edge_selection[,1] != Ntip(tree)+1,]
        # not always possible (with the less descendants)
        edge_selection <- edge_selection[order(sapply(edge_selection[,2], function(x) length(Descendants(tree, x)[[1]]))),]
        
        # CHANGE THIS
        if(i > 1){
          if(subtree_bck_loc %in% as.numeric(names(shift_nodes))){
            shift_nodes1 <- shift_nodes[names(shift_nodes) == as.character(subtree_bck_loc)] # of this tree
            shift_nodes_and_anc <- c(shift_nodes1, unique(unlist(sapply(shift_nodes1[!is.na(shift_nodes1)], function(x) Ancestors(tree, x)))))
          }
        }
        
        edge_selection <- edge_selection[!edge_selection[,2] %in% shift_nodes_and_anc,]
        
        if(is.vector(edge_selection)){
          if(!edge_selection[2] %in% shift_nodes_and_anc){
            node <- edge_selection[2]
          } else {
            node <- NA
          }
        } else {
          if(nrow(edge_selection) != 0){
            node <- edge_selection[1,2]
          } else {
            node <- NA
          }
        }
        
        shift_nodes <- c(shift_nodes, node)
        names(shift_nodes)[i] <- subtree_bck_loc
      }
      
    }
    all_shift_nodes[[j]] <- shift_nodes
  }
  
  sum(sapply(all_shift_nodes, anyNA) == F)
  
  all_shift_nodes_ok <- all_shift_nodes[sapply(all_shift_nodes, anyNA) == F]
  all_parts_ok <- all_parts[sapply(all_shift_nodes, anyNA) == F]
  
  # grafting on a backbone
  cat("Shifts grafting...\n")
  all_new_tree <- c()
  pb = txtProgressBar(min = 0, max = length(all_parts_ok), initial = 0, style = 3)
  for(j in 1:length(all_parts_ok)){
    setTxtProgressBar(pb,j)
    all_loc_j <- unique(as.numeric(names(all_shift_nodes_ok[[j]])))
    
    for(bck in all_loc_j){
      tree <- all_parts_ok[[j]][[bck]]
      tree$node.label <- c(Ntip(tree)+1):c(Ntip(tree)+Nnode(tree))
      tree_pruned1 <- tree
      for(i in 1:length(all_shift_nodes_ok[[j]])){
        selected_node <- all_shift_nodes_ok[[j]][i]
        if(as.numeric(names(selected_node)) == bck){
          if(Ntip(tree) != Ntip(tree_pruned1)){
            if(length(Descendants(tree, selected_node)[[1]]) > 1){
              selected_node <- getMRCA(tree_pruned1, tree$tip.label[Descendants(tree, selected_node)[[1]]])
            }else{
              selected_node <- which(tree_pruned1$tip.label %in% tree$tip.label[selected_node])
            }
            
          }
          tree_pruned2<-splitTree(tree_pruned1,split=list(node=selected_node,
                                                          bp=branching.times(tree)[as.character(Ancestors(tree, all_shift_nodes_ok[[j]][i], "parent"))] - all_ages[i]))[[1]]
          tree_pruned2$tip.label[tree_pruned2$tip.label=="NA"]<-paste0("shift.", strsplit(all_parts_ok[[j]][[i]]$tip.label, "")[[1]][1])
          tree_pruned1 <- tree_pruned2
        }
      }
      
      for(i in 1:length(all_shift_nodes_ok[[j]])){
        selected_node <- all_shift_nodes_ok[[j]][i]
        if(as.numeric(names(selected_node)) == bck){
          new_tree<-ladderize(bind.tree(tree_pruned1,all_parts_ok[[j]][[i]],where=which(tree_pruned1$tip.label==paste0("shift.", strsplit(all_parts_ok[[j]][[i]]$tip.label, "")[[1]][1]))))
          tree_pruned1 <- new_tree
        }
      }
      all_parts_ok[[j]][[bck]] <- new_tree
    }
    all_new_tree[[j]] <- ladderize(new_tree)
  } 
  close(pb)
  
  to_keep <- sapply(all_new_tree, function(x){
    
    ntip_by_group <- table(sapply(strsplit(x$tip.label, ""), "[[", 1))
    backbone_letters <- letters[(length(letters)-length(backbones)+1):length(letters)]
    monophyly <- sapply(names(ntip_by_group), function(y) is.monophyletic(x, x$tip.label[grepl(y, x$tip.label)]))
    monophyly_backbones <- monophyly[backbone_letters]
    checks <- c(all(backbone_letters %in% names(ntip_by_group)),
                all(monophyly_backbones == F),
                all(ntip_by_group >= clade.size))
    ifelse(any(checks == F), F, T)
  })
  all_new_tree <- all_new_tree[to_keep]
  
  return(all_new_tree)
  
}