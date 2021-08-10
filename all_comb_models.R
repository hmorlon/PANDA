all_comb_models <- function(to){
  
  if(multi.backbone == F){
    names_ALL_bck_comb <- names(ALL_bck_comb)
    ALL_bck_comb <- rep(list(NULL), length(ALL_bck_comb))
    names(ALL_bck_comb) <- names_ALL_bck_comb
  }
  ALL_comb <- strsplit(names(ALL_bck_comb), "[.]")
  
  ALL_bck_comb_to <- ALL_bck_comb[names(ALL_bck_comb) %in% paste(ALL_comb[[to]],collapse = ".")][[1]]
  
  clade_to_shift <- ALL_comb[[to]]
  cat("\n", to, "/", length(ALL_comb))
  shift <- ALL_branch_times_clades[c(clade_to_shift,unlist(ALL_bck_comb_to))]
  
  # USELESS??? ####
  
  # to check if working with multiple backbone
  branch_times_to <- get.branching.nodes(shift)
  #if(length(branch_times_to) > 1){
  branch_to_bck <- lapply(all_tested_nodes, function(x) get.branching.nodes(ALL_branch_times_clades[x]))
  #} else {
  #  branch_to_bck <- lapply(all_tested_nodes, function(x) get.branching.nodes(ALL_branch_times_clades[x]))
  #}
  
  # ~ Branching TIMES ##### 
  # Branching times should be created from branch_times because there are created at each combination
  branch_times <- branch_times_to
  for(nodeID in 1:length(branch_times_to)){
    branch_times[[nodeID]] <- sapply(branch_times_to[[nodeID]], get.node.ages)
  }
  
  ALL_branch_times_to_bck <- branch_to_bck
  
  for(bID in 1:length(ALL_branch_times_to_bck)){
    branch_times_nodeID <- NULL
    for(nodeID in 1:length(ALL_branch_times_to_bck[[bID]])){
      branch_times_nodeID[[nodeID]] <- sapply(ALL_branch_times_to_bck[[bID]][[nodeID]], get.node.ages)
    }
    ALL_branch_times_to_bck[[bID]] <- branch_times_nodeID
  }
  
  #ALL_branch_times_to_bck[length(ALL_branch_times_to_bck) + 1] <- list(branch_times)
  names(ALL_branch_times_to_bck) <- all_tested_nodes
  
  # Muli backbone
  
  # create the backbone 
  # MULTI-BACKBONE
  if(is.null(ALL_bck_comb_to)){
    ALL_multi_bck_to <- rep(list(NULL), 1)
    phylo_backbone <- rep(list(NULL), 1)
    branch_times_backbone <- rep(list(NULL), 1)
    
  } else {
    ALL_multi_bck_to <- rep(list(NULL), length(ALL_bck_comb_to))
    names(ALL_multi_bck_to) <- sapply(ALL_bck_comb_to, function(x) paste(x, collapse = "."))
    phylo_backbone <- rep(list(NULL), length(ALL_bck_comb_to))
    names(phylo_backbone) <- sapply(ALL_bck_comb_to, function(x) paste(x, collapse = "."))
    branch_times_backbone <- rep(list(NULL), length(ALL_bck_comb_to))
    names(branch_times_backbone) <- sapply(ALL_bck_comb_to, function(x) paste(x, collapse = "."))
    
  }
  
  for(mb in 1:length(ALL_multi_bck_to)){
    
    if(is.null(ALL_bck_comb_to[[mb]])){
      # No split in the backbone
      
      phylo_backbone_cut <- list(drop.tip(phy, unlist(ALL_clade_names[clade_to_shift])))
      names(phylo_backbone_cut) <- paste0(paste(c(clade_to_shift),collapse = "."),"_bck")
      #branch_times_to_bck <- list(unlist(ALL_branch_times_to_bck[clade_to_shift], recursive = F))
      
      branch_times_to_bck <- list(branch_times) # reuse branch_times for other parts
      names(branch_times_to_bck) <- paste(clade_to_shift, collapse = ".")
      int_nodes <- NULL
      
    } else {
      # Split in the backbone
      int_nodes <- ALL_bck_comb_to[[mb]]
      # order from present to past
      int_nodes <- names(branching.times(phy)[order(branching.times(phy))])[names(branching.times(phy)[order(branching.times(phy))]) %in% int_nodes]
      branch_times_to_bck <- rep(list(NULL), length(int_nodes)+1)
      phylo_backbone_cut <- rep(list(NULL), length(int_nodes)+1)
      phylo_backbone_core <- drop.tip(phy, unlist(ALL_clade_names[clade_to_shift]))
      sb.tips <- rep(list(NULL), length(int_nodes))
      
      for(sb in 1:length(int_nodes)){
        
        sb.tips[[sb]] <- phylo_backbone_core$tip.label[phylo_backbone_core$tip.label %in% unlist(ALL_clade_names[int_nodes[sb]])]
        #sb.tips[[sb]] <- sb.tips[[sb]][!sb.tips[[sb]]  %in% sb.tips[[-sb]]]
        phylo_backbone_cut[[sb]] <- subtree(phylo_backbone_core, sb.tips[[sb]])
        names(phylo_backbone_cut)[[sb]] <- paste(int_nodes[sb],"sub", sep = "_")
        
        sb.desc <- Descendants(phy, as.numeric(int_nodes[sb]), type = "all")
        desc_sub_up_bck <- c(clade_to_shift[clade_to_shift %in% sb.desc],int_nodes[int_nodes %in% sb.desc])
        
        if(sb > 1){
          desc_sub_up_bck <- desc_sub_up_bck[!desc_sub_up_bck %in% unlist(strsplit(names(branch_times_to_bck)[!is.na(names(branch_times_to_bck))], "[.]"))]
          phylo_backbone_cut[[sb]] <- drop.tip(phylo_backbone_cut[[sb]], unlist(sb.tips[-sb]))
          child <- Children(phy, as.numeric(int_nodes[sb]))
          
          if(as.numeric(int_nodes[sb-1]) %in% child){ # Direct Descendant
            
            sis <- Siblings(phy, as.numeric(int_nodes[sb-1]))
            
            branch_time_sb <- get.branching.nodes(ALL_branch_times_clades[c(desc_sub_up_bck, as.character(sis))], root_ID = int_nodes[sb])
            
            for(nodeID in 1:length(branch_time_sb)){
              branch_time_sb[[nodeID]] <- sapply(branch_time_sb[[nodeID]], get.node.ages)
            }
            
            branch_times_to_bck[sb] <- list(branch_time_sb)
            
            names(branch_times_to_bck)[sb] <- paste(c(desc_sub_up_bck), collapse = ".")
          } else {
            
            branch_time_sb <- get.branching.nodes(ALL_branch_times_clades[desc_sub_up_bck], root_ID = int_nodes[sb])
            
            for(nodeID in 1:length(branch_time_sb)){
              branch_time_sb[[nodeID]] <- sapply(branch_time_sb[[nodeID]], get.node.ages)
            }
            
            branch_times_to_bck[sb] <- list(branch_time_sb)
            names(branch_times_to_bck)[sb] <- paste(c(desc_sub_up_bck), collapse = ".")
          }
          
        } else{
          
          branch_time_sb <- get.branching.nodes(ALL_branch_times_clades[desc_sub_up_bck], root_ID = int_nodes[sb])
          
          # check that root of phylo_backbone_cut[[sb]] is int_node
          if(phylo_backbone_cut[[sb]]$node.label[1] != int_nodes[sb] &
             !phylo_backbone_cut[[sb]]$node.label[1] %in% names(branch_time_sb)){
            
            root_sb_to_int_nodes <- c(phylo_backbone_cut[[sb]]$node.label[1], Ancestors(phy, phylo_backbone_cut[[sb]]$node.label[1]))
            root_sb_to_int_nodes <- root_sb_to_int_nodes[1:c(which(root_sb_to_int_nodes == int_nodes[sb])-1)]
            missed_sb_nodes <- root_sb_to_int_nodes[!root_sb_to_int_nodes %in% as.numeric(names(branch_time_sb))]
            
            for(msb in 1:length(missed_sb_nodes)){
              branch_time_sb[length(branch_time_sb)+msb] <- unlist(ALL_branch_times_clades[as.character(missed_sb_nodes[msb])], recursive = F)
              names(branch_time_sb)[length(branch_time_sb)]<- as.character(missed_sb_nodes[msb])
            }
          }
          
          for(nodeID in 1:length(branch_time_sb)){
            branch_time_sb[[nodeID]] <- sapply(branch_time_sb[[nodeID]], get.node.ages)
          }
          
          branch_times_to_bck[sb] <- list(branch_time_sb)
          names(branch_times_to_bck)[sb] <- paste(c(clade_to_shift[clade_to_shift %in% Descendants(phy, as.numeric(int_nodes[sb]), "all")]), collapse = ".")
        }
        
        if(sb == length(int_nodes)){
          tips_up_bck <- unlist(lapply(phylo_backbone_cut, function(x) x$tip.label))
          tips_last_bck <- unlist(ALL_clade_names[clade_to_shift[!clade_to_shift %in% sb.desc]])
          
          phylo_backbone_cut[[sb+1]] <- drop.tip(phylo_backbone_core, tips_up_bck)
          names(phylo_backbone_cut)[sb+1] <- paste(int_nodes[sb],"bck", sep = "_")
          
          sb1.desc <- Descendants(phy, Ntip(phy)+1, type = "all")
          sb1.desc <- sb1.desc[!sb1.desc %in% sb.desc]
          sb1.desc <- sb1.desc[sb1.desc > Ntip(phy)]
          
          if(!is.null(tips_last_bck)){  
            # subgroup(s) in the deep backbone
            btt_bck <- c(clade_to_shift[clade_to_shift %in% sb1.desc], ALL_bck_comb_to[[mb]][ALL_bck_comb_to[[mb]] %in% sb1.desc])
            btt_bck <- btt_bck[!btt_bck %in% unlist(strsplit(names(branch_times_to_bck), "[.]"))]
            
            # Using get.branching.nodes to be sure to get all branches
            
            branch_time_deep <- get.branching.nodes(ALL_branch_times_clades[btt_bck], root_ID = Ntip(phy)+1)
            
            for(nodeID in 1:length(branch_time_deep)){
              branch_time_deep[[nodeID]] <- sapply(branch_time_deep[[nodeID]], get.node.ages)
            }
            
            branch_times_to_bck[sb+1] <- list(branch_time_deep)
            
            names(branch_times_to_bck)[sb+1] <- paste(btt_bck, collapse = ".")
            
          } else {
            # no subgroup in the deep backbone
            branch_times_to_bck[sb+1] <- list(ALL_branch_times_to_bck[[ALL_bck_comb_to[[mb]][ALL_bck_comb_to[[mb]] %in% sb1.desc]]])
            names(branch_times_to_bck)[sb+1] <- paste(ALL_bck_comb_to[[mb]][ALL_bck_comb_to[[mb]] %in% sb1.desc], collapse = ".")
            
          }
          
        }
      }
    }
    
    if(multi.backbone == T & !is.null(ALL_bck_comb_to)){
      cat("\n",mb, "/",length(ALL_bck_comb_to), "multibackbone(s)")
    } else {
      cat("\n",mb, "/",length(phylo_backbone_cut), "backbone")
    }
    
    # tot time et spec time
    res_bck <- rep(list(NULL),length(phylo_backbone_cut))
    
    # Sampling fractions ####
    
    lin.node <- data.frame(node = c(ALL_comb[[to]],int_nodes, Ntip(phy)+1), n.tips = rep(NA, length(ALL_comb[[to]]) + length(int_nodes)+1))
    lin.node$node <- as.character(lin.node$node)
    lin.node <- merge(lin.node, sampling.fractions[sampling.fractions$nodes %in% lin.node$node, c("nodes", "sp_tt"),],
                      by.x = "node", by.y = "nodes")
    
    node_order <- names(branching.times(phy)[order(branching.times(phy))])
    node_order <- node_order[node_order %in% lin.node$node]
    
    lin.node <- lin.node[match(node_order, lin.node$node),]
    
    for(n.lin in 1:nrow(lin.node)){
      desc.n.lin <- length(Descendants(phy, as.numeric(lin.node$node[n.lin]))[[1]])
      # whether this node is present in an other lineage
      int.n.lin <- Descendants(phy, as.numeric(lin.node$node[n.lin]), type = "all")
      int.n.lin <- as.character(int.n.lin[int.n.lin > Ntip(phy)])
      # Ntip
      if(any(ALL_comb[[to]] %in% int.n.lin)){
        lin.node$n.tips[n.lin] <- desc.n.lin - sum(lin.node$n.tips[lin.node$node %in% ALL_comb[[to]][ALL_comb[[to]] %in% int.n.lin]])
        lin.node$sp_tt[n.lin] <- lin.node$sp_tt[n.lin] - sum(lin.node$sp_tt[lin.node$node %in% ALL_comb[[to]][ALL_comb[[to]] %in% int.n.lin]])
      } else{
        lin.node$n.tips[n.lin] <- desc.n.lin
      }
    }
    
    lin.node$n.tips_prev <- lin.node$n.tips
    lin.node$sp_tt_prev <- lin.node$sp_tt
    
    lin.node_bck <- lin.node[!lin.node$node %in% ALL_comb[[to]],]
    
    for(l.n in c(1:nrow(lin.node_bck))){
      int.desc_lin <- unlist(Descendants(phy, as.numeric(lin.node_bck$node[l.n]), "all"))
      int.desc_lin <- int.desc_lin[int.desc_lin > Ntip(phy)]
      
      if(any(lin.node_bck$node %in% int.desc_lin)){
        
        bck_up <- lin.node_bck[which(lin.node_bck$node %in% int.desc_lin),]
        
        ntip_bck_up <- sum(bck_up$n.tips_prev)
        ntaxo_bck_up <- sum(bck_up$sp_tt_prev)
        
        lin.node_bck$n.tips_prev[l.n] <- lin.node_bck$n.tips[l.n] - ntip_bck_up
        lin.node_bck$sp_tt_prev[l.n] <- lin.node_bck$sp_tt[l.n] - ntaxo_bck_up  
        
      }
    }
    lin.node[lin.node$node %in% lin.node_bck$node,] <- lin.node_bck
    
    lin.node <- lin.node[-(1:length(ALL_comb[[to]])),]
    
    f <- as.list(lin.node$n.tips_prev/lin.node$sp_tt_prev)
    names(f) <- names(phylo_backbone_cut)
    
    for(btb in 1:length(phylo_backbone_cut)){
      if(multi.backbone == T & !is.null(ALL_bck_comb_to)){
        cat("\n -",btb, "/",length(phylo_backbone_cut), "sub backbones\n")
      }
      
      # by default backbone.option = backbone2
      backbone <- backbone.option
      spec_times <- NULL
      cond <- "crown"
      
      # CHECKED!
      tot_time3 <- max(c(node.age(phylo_backbone_cut[[btb]])$ages, unlist(branch_times_to_bck[[btb]])))
      
      # for converting in backbone1: not fully working yet
      if(backbone.option == "backbone1"){
        
        spec_times <- sapply(branch_times_to_bck[[btb]], "[[", 2)
        cond <- "stem"
        
        # conditioning backbone at root
        if(length(grep("_bck", names(phylo_backbone_cut[btb]))) == 1){
          cond <- "crown"
        }
        
        branch_times_to_bck[[btb]] <- rep(list(NULL),1)
        
        if(!is.null(phylo_backbone_cut[[btb]]$root.edge)){
          tot_time3 <- tot_time3[[1]] + phylo_backbone_cut[[btb]]$root.edge
        }
      }
      
      #####################################
      
      results <- div.models(phylo = phylo_backbone_cut[[btb]], tot_time = tot_time3, f = f[[btb]],
                            backbone = backbone, spec_times = spec_times, branch_times = branch_times_to_bck[[btb]],
                            cond = cond, models = models, n.max = n.max, rate.max = rate.max, verbose = T)
      if(btb < length(phylo_backbone_cut)){
        # cond has to be changed to properly estimate likelihood of each part if they are not the last part
        results1 <- div.models(phylo = phylo_backbone_cut[[btb]], tot_time = tot_time3, f = f[[btb]],
                               backbone = backbone, spec_times = spec_times, branch_times = branch_times_to_bck[[btb]],
                               cond = F, models = models, n.max = n.max, rate.max = rate.max, verbose = F)
        
        results2 <- merge(results1[,c(1:4)], results[,c(1,5:8)], by="Models")
        results <- results2[match(results$Models, results2$Models),]
        
      }
      
      results[,-1] <- apply(results[,-1], 2, as.numeric)
      res_bck[btb] <- list(results)
    }
    
    names(res_bck) <- names(phylo_backbone_cut)
    ALL_multi_bck_to[mb] <- list(res_bck)
    
  } # END BCK_COMB
  #names(ALL_multi_bck_to) <- paste(clade_to_shift, collapse = ".")
  return(ALL_multi_bck_to)
  # Multi merge
}