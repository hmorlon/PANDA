all_comb_models <- function(to){
  
  # splitting combination into subclades and backbones
  
  comb1 <- strsplit(comb.shift[to], "/")[[1]]
  comb.sub <- strsplit(comb1[[1]], "[.]")[[1]]
  if(length(comb1) == 2){
    comb.bck <- strsplit(comb.shift[to], "/")[[1]][2]
    comb.bck <- strsplit(comb.bck, "[.]")[[1]]
  } else {
    comb.bck <- NULL
  }
  
  message("\n", to, "/", length(comb.shift))
  
  # plot to illustrate
  # plot_phylo_comb(phylo, data, sampling.fractions, comb = comb.shift[to], cex = 0.8, label.offset = 0.2)
  # nodes <- sampling.fractions$nodes[!is.na(sampling.fractions$to_test)]
  # nodelabels(as.character(nodes), nodes)
  
  # create the backbone
  # new way 
 
  int_nodes <- comb.bck
  # order from present to past
  int_nodes <- names(branching.times(phylo)[order(branching.times(phylo))])[names(branching.times(phylo)[order(branching.times(phylo))]) %in% int_nodes]
  branch_times_to_bck <- rep(list(NULL), length(comb.bck)+1)
  phylo_backbone_cut <- rep(list(NULL), length(comb.bck)+1)
  phylo_backbone_core <- drop.tip(phylo, unlist(ALL_clade_names[comb.sub]))
  
  res_bck <- rep(list(NULL), length(comb.bck)+1)
  
  sb.tips <- rep(list(NULL), length(int_nodes))
  sb.desc <- rep(list(NULL), length(int_nodes))
  names(sb.desc) <- int_nodes
  
  for(sb in 1:length(res_bck)){
    
    if(is.null(comb.bck)){ # simple backbone
      
      phylo_backbone_cut <- list(phylo_backbone_core)
      names(phylo_backbone_cut) <- paste0(paste0(comb.sub, collapse = "."),"_bck")
      
      branch_time_sb <- get.branching.nodes(comb.sub, phylo = phylo,
                                            ALL_branch_times_clades = ALL_branch_times_clades,
                                            ALL_clade_names = ALL_clade_names)
      branch_times_to_bck <- list(branch_time_sb)
      names(branch_times_to_bck) <- paste0(comb.sub, collapse = ".")
      
      # check the root? seems ok with parnassiinae
      
    } else { # multibackbone
      
      if(sb < length(phylo_backbone_cut)){ # before deep backbone
        
        sb.desc[[sb]] <- Descendants(phylo, as.numeric(int_nodes[sb]), "all")
        
        if(sb > 1){ # removing descendant in previous int_nodes
          sb.desc[[sb]] <- sb.desc[[sb]][!sb.desc[[sb]] %in% unlist(sb.desc[1:c(sb-1)])]
        }
        
        sb.desc_sb_sp <- phylo$tip.label[sb.desc[[sb]][sb.desc[[sb]] < Ntip(phylo)]]
        sb.desc_sb_sp <- intersect(sb.desc_sb_sp, phylo_backbone_core$tip.label)
        phylo_backbone_cut[[sb]] <- subtree(phylo_backbone_core, sb.desc_sb_sp)
        names(phylo_backbone_cut)[sb] <- paste0(int_nodes[sb],"_sub")
        
        comb.multibackbone <- c(comb.sub[comb.sub %in% sb.desc[[sb]]], int_nodes[int_nodes %in% sb.desc[[sb]]])
        
        branch_time_sb <- get.branching.nodes(comb.multibackbone, phylo = phylo,
                                              ALL_branch_times_clades = ALL_branch_times_clades,
                                              ALL_clade_names = ALL_clade_names)
        
        # check that root of phylo_backbone_cut[[sb]] is int_node
        if(phylo_backbone_cut[[sb]]$node.label[1] != int_nodes[sb] &
           !phylo_backbone_cut[[sb]]$node.label[1] %in% names(branch_time_sb)){
          
          root_sb_to_int_nodes <- c(phylo_backbone_cut[[sb]]$node.label[1], Ancestors(phylo, phylo_backbone_cut[[sb]]$node.label[1]))
          root_sb_to_int_nodes <- root_sb_to_int_nodes[1:c(which(root_sb_to_int_nodes == int_nodes[sb])-1)]
          missed_sb_nodes <- root_sb_to_int_nodes[!root_sb_to_int_nodes %in% as.numeric(names(branch_time_sb))]
          
          for(msb in 1:length(missed_sb_nodes)){
            branch_time_missing_sb <- list(c(missed_sb_nodes[msb], Ancestors(phylo, missed_sb_nodes[msb], "parent")))
            names(branch_time_missing_sb) <- missed_sb_nodes[msb]
            
            branch_time_sb[length(branch_time_sb)+1] <- branch_time_missing_sb
            names(branch_time_sb)[length(branch_time_sb)]<- as.character(missed_sb_nodes[msb])
          }
        }
        
        branch_times_to_bck[sb] <- list(branch_time_sb)
        names(branch_times_to_bck)[sb] <- paste(comb.multibackbone, collapse = ".")
        
      } else {  # deep backbone
        
        tips_up_bck <- unlist(lapply(phylo_backbone_cut, function(x) x$tip.label))
        # remaining comb.sub in the deep backbone
        #tips_last_bck <- unlist(ALL_clade_names[comb.sub[!comb.sub %in% unlist(sb.desc, use.names = FALSE)]])
        
        phylo_backbone_cut[[sb]] <- drop.tip(phylo_backbone_core, tips_up_bck)
        names(phylo_backbone_cut)[sb] <- paste(int_nodes[sb-1],"bck", sep = "_")
        
        int_nodes_deep_backbone <- int_nodes[!int_nodes %in% unlist(sapply(branch_times_to_bck, names), use.names = FALSE)]
 
        comb_deep_backbone <- c(comb.sub[!comb.sub %in% unlist(sb.desc, use.names = FALSE)], int_nodes_deep_backbone)
  
        branch_time_sb <- get.branching.nodes(comb_deep_backbone, phylo = phylo,
                                              ALL_branch_times_clades = ALL_branch_times_clades,
                                              ALL_clade_names = ALL_clade_names)
  
        branch_times_to_bck[sb] <- list(branch_time_sb)
        names(branch_times_to_bck)[sb] <- paste(comb_deep_backbone, collapse = ".")
        
      } # deep backbone
    } # multi backbone
  }
  
  branch_nodes_to_bck <- branch_times_to_bck
  for(bck in 1:length(branch_times_to_bck)){
    for(nodeID in 1:length(branch_nodes_to_bck[[bck]])){
      branch_times_to_bck[[bck]][[nodeID]] <- sapply(branch_nodes_to_bck[[bck]][[nodeID]], get.node.ages, phylo = phylo)
    }  
  }
  
  # checking all branches are taken into account
  desc_comb.sub <-  Descendants(phylo, as.numeric(comb.sub), "all")
  desc_comb.sub <- lapply(desc_comb.sub, function(x) x[x > Ntip(phylo)])
  
  nodes_backbone_th <- setdiff(phylo$node.label, unlist(desc_comb.sub))
  
  nodes_backbone_obs <- unlist(lapply(phylo_backbone_cut, function(x) x$node.label), use.names = FALSE)
  all_branching_nodes_to <- unlist(lapply(branch_nodes_to_bck, function(x) unique(sapply(x, "[[", 2))), use.names = FALSE)
  branch_nodes_to_bck <- unlist(lapply(branch_nodes_to_bck, names), use.names = FALSE)
  
  nodes_backbone_obs <- as.numeric(c(nodes_backbone_obs,
                                     branch_nodes_to_bck,
                                     all_branching_nodes_to))
  
  if(!all(nodes_backbone_th %in% unique(nodes_backbone_obs))){
    stop("\n#### Some branches are missing... ####\n")
  }
  
  # Sampling fractions ####
  
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
  names(f) <- names(phylo_backbone_cut)
  
  for(btb in 1:length(phylo_backbone_cut)){
    
    # by default backbone.option = "crown.shift"
    backbone <- backbone.option
    spec_times <- NULL
    cond <- "crown"
    
    # CHECKED!
    tot_time3 <- max(c(node.age(phylo_backbone_cut[[btb]])$ages, unlist(branch_times_to_bck[[btb]])))
    
    # for converting in stem.shift
    if(backbone.option == "stem.shift"){
      
      spec_times <- sapply(branch_times_to_bck[[btb]], "[[", 2)
      cond <- "stem"
      
      if(!is.null(phylo_backbone_cut[[btb]]$root.edge)){
        tot_time3 <- max(node.age(phylo_backbone_cut[[btb]])$ages) + phylo_backbone_cut[[btb]]$root.edge
      }
      
      # if deep backbone, conditioning backbone at crown 
      if(length(grep("_bck", names(phylo_backbone_cut[btb]))) == 1){
        cond <- "crown"
      }
      branch_times_to_bck[[btb]] <- rep(list(NULL),1)
    }
    
    ##################################### models
    
    results <- div.models(phylo = phylo_backbone_cut[[btb]], tot_time = tot_time3, f = f[[btb]],
                          backbone = backbone, spec_times = spec_times, branch_times = branch_times_to_bck[[btb]],
                          cond = cond, models = models, n.max = n.max, rate.max = rate.max, verbose = TRUE)
    if(btb < length(phylo_backbone_cut)){
      # cond has to be changed to properly estimate likelihood of each part if they are not the last part
      results1 <- div.models(phylo = phylo_backbone_cut[[btb]], tot_time = tot_time3, f = f[[btb]],
                             backbone = backbone, spec_times = spec_times, branch_times = branch_times_to_bck[[btb]],
                             cond = FALSE, models = models, n.max = n.max, rate.max = rate.max, verbose = FALSE)
      
      results2 <- merge(results1[,c(1:4)], results[,c(1,5:8)], by="Models")
      results <- results2[match(results$Models, results2$Models),]
      
      # adding a parameter for the location of the shift (to modify for the printing)
      results$AICc <- 2 * -results$logL + 2 * (results$Parameters+1) + (2 * (results$Parameters+1) * ((results$Parameters+1) + 1))/(Ntip(phylo_backbone_cut[[btb]]) - (results$Parameters+1) - 1)
      results$Parameters <- results$Parameters+1
    }
    
    results[,-1] <- apply(results[,-1], 2, as.numeric)
    res_bck[btb] <- list(results)
  }
  
  names(res_bck) <- names(phylo_backbone_cut)

  return(res_bck)
  # Multi merge
}
