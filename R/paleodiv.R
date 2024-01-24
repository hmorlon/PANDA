paleodiv <- function(phylo, data, sampling.fractions, shift.res,
                     backbone.option = "crown.shift", combi = 1,
                     split.div = F){

  # Checking arguments ####
  # phylo
  if(!inherits(phylo, "phylo")){
    stop("object \"phylo\" is not of class \"phylo\"")
  } else {
    phylo$node.label <- c(c(Ntip(phylo)+1):c(Ntip(phylo)+Nnode(phylo)))
  }
  
  # data
  if(!inherits(data, "data.frame")){
    stop("object \"data\" is not of class \"data.frame\"")
  }
  # sampling.fractions
  if(phylo$Nnode + Ntip(phylo) != nrow(sampling.fractions) | is(sampling.fractions)[1]!="data.frame"){
    stop("object \"sampling.fractions is not of class \"data.frame\" or is do not correspond to the provided phylogeny")
  }
  # shift.res
  if(!is(shift.res)[1] == "list" | any(names(shift.res) != c("whole_tree", "subclades", "backbones", "total"))){
    stop("object \"shift.res\" might be incorrect.")
  }
  if(!is.numeric(combi)){
    stop("object \"combi\" should be numeric.")
  }
  if(is(split.div)[1] != "logical"){
    stop("object \"split.div\" should be logical.")
  }
  
  if(!backbone.option %in% c("stem.shift", "crown.shift")){
    cat("\nArgument \"backbone.option\" is incorrect.")
    stop()
  }
  
  best_subclades_df <- do.call(rbind.data.frame, lapply(shift.res$subclades, function(x) x[1,]))
  best_subclades_df$Clades <- row.names(best_subclades_df)
  row.names(best_subclades_df) <- NULL
  best_subclades_df <- best_subclades_df[,c(10,1:8)]
  
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
  
  tot_time <- max(branching.times(phylo))
  totalsp <- list(nrow(data))
  
  time.seq <- c(tot_time, seq(floor(tot_time),0,by=-1))
  globaldiv <- matrix(NA,length(comb.sub)+length(comb.bck)+ifelse(comb.sub != "whole_tree",1,0), length(time.seq)) #matrix(NA, Number of clades, Crown age of the whole tree + 1)
  
  if(any(comb.sub == "whole_tree")){
    
    best_whole_tree_combi <- shift.res$whole_tree
    best_whole_tree_combi <- best_whole_tree_combi[best_whole_tree_combi$AICc == min(best_whole_tree_combi$AICc),]
    
    model <- as.character(best_whole_tree_combi$Models[1])
    values <- as.numeric(best_whole_tree_combi[1,-c(1,2,9)])
    names(values) <- names(best_whole_tree_combi[1,-c(1,2,9)])
    
    lamb_pari <- as.numeric(c(values["Lambda"],values["Alpha"]))
    mu_pari <- as.numeric(c(values["Mu"],values["Beta"]))
    agei <- tot_time
    sizei <- totalsp[[1]]
    time_seq <- c(agei, seq(floor(agei),0,by=-1))
    
    if (grepl("BCST", model)){
      div <-sizei*exp(-abs(lamb_pari[1])*time_seq)
    }
    
    if (grepl("BCST_DCST", model)){
      div<-sizei*exp(-abs(lamb_pari[1])*time_seq+abs(mu_pari[1])*time_seq)
    }
    
    if (grepl("BVAR", model)){
      div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq)))
    }
    
    if (grepl("BVAR_DCST", model)){
      div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq))+abs(mu_pari[1])*time_seq)
    }
    
    if (grepl("BCST_DVAR", model)){
      div<-sizei*exp(-abs(lamb_pari[1])*time_seq-abs(mu_pari[1])/mu_pari[2]*(1-exp(mu_pari[2]*time_seq)))
    }
    
    if (grepl("BVAR_DVAR", model)){
      div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq))-abs(mu_pari[1])/mu_pari[2]*(1-exp(mu_pari[2]*time_seq)))
    }
    
    div2 <- div[length(div):1]
    globaldiv[1,1:length(div)]<-div2
    
    
  } else {
    
    best_subclades_df_combi <- best_subclades_df[best_subclades_df$Clades %in% as.numeric(comb.sub),]
    best_subclades_df_combi <- best_subclades_df_combi[match(comb.sub, best_subclades_df_combi$Clades), ]
    
    if(backbone.option == "stem.shift"){
      parental_nodes <- Ancestors(phylo, as.numeric(best_subclades_df_combi$Clades), type = "parent")
      tot_time2 <- as.list(branching.times(phylo)[as.character(parental_nodes)])
    } else {
      tot_time2 <- as.list(branching.times(phylo)[best_subclades_df_combi$Clades])
    }
    
    totalsp2 <- as.list(sampling.fractions$sp_tt[sampling.fractions$nodes %in% as.numeric(comb.sub)])
    names(totalsp2) <- comb.sub
    names(tot_time2) <- comb.sub
  }
  
  # Subclades diversity (RPANDA FONCTIONS !!!!)
  
  
  if(all(comb.sub != "whole_tree")){
    # Backbone diversity
    
    for(i in 1:nrow(best_subclades_df_combi)){
      clade <- as.character(best_subclades_df_combi$Clades[i])
      model <- as.character(best_subclades_df_combi$Models[i])
      values <- as.numeric(best_subclades_df_combi[i,-c(1,2,10)])
      names(values) <- names(best_subclades_df_combi[i,-c(1,2,10)])
      
      lamb_pari <- as.numeric(c(values["Lambda"],values["Alpha"]))
      mu_pari <- as.numeric(c(values["Mu"],values["Beta"]))
      agei <- tot_time2[[clade]]
      sizei <- totalsp2[[clade]]
      time_seq <- c(agei, seq(floor(agei),0,by=-1))
      
      if (grepl("BCST", model)){
        div <-sizei*exp(-abs(lamb_pari[1])*time_seq)
      }
      
      if (grepl("BCST_DCST", model)){
        div<-sizei*exp(-abs(lamb_pari[1])*time_seq+abs(mu_pari[1])*time_seq)
      }
      
      if (grepl("BVAR", model)){
        div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq)))
      }
      
      if (grepl("BVAR_DCST", model)){
        div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq))+abs(mu_pari[1])*time_seq)
      }
      
      if (grepl("BCST_DVAR", model)){
        div<-sizei*exp(-abs(lamb_pari[1])*time_seq-abs(mu_pari[1])/mu_pari[2]*(1-exp(mu_pari[2]*time_seq)))
      }
      
      if (grepl("BVAR_DVAR", model)){
        div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq))-abs(mu_pari[1])/mu_pari[2]*(1-exp(mu_pari[2]*time_seq)))
      }
      
      div2 <- div[length(div):1]
      globaldiv[i,1:length(div)]<-div2
    }
    
    best_backbones <- shift.res$backbones[paste(paste(comb.sub, collapse = "."), paste(comb.bck, collapse = "."), sep = "/")][[1]]
    best_backbones_df <- do.call(rbind.data.frame, lapply(best_backbones, function(x) x[1,]))
    best_backbones_df$parts <- row.names(best_backbones_df)
    row.names(best_backbones_df) <- NULL
    best_backbones_df <- best_backbones_df[,c(10,1:8)]
    all_tested_nodes <- c(comb.sub, comb.bck)
    
    ALL_clade_names <- rep(list(NULL), length(all_tested_nodes))
    
    for(pot_names in 1:length(ALL_clade_names)){
      ALL_clade_names[pot_names] <- list(phylo$tip.label[unlist(Descendants(phylo, as.numeric(all_tested_nodes[pot_names])))])
    }
    names(ALL_clade_names) <- all_tested_nodes
    ALL_nodes_ages <- as.data.frame(apply(data.frame(nodesID=names(branching.times(phylo)),ages=branching.times(phylo)), 2, as.numeric))
    
    ALL_branch_times_clades <- rep(list(NULL),length(all_tested_nodes))
    names(ALL_branch_times_clades) <- all_tested_nodes
    
    for(clade in 1:length(all_tested_nodes)){
      
      parental_node <- Ancestors(phylo, as.numeric(all_tested_nodes[clade]), type = "parent")
      branch_times_clade <- unlist(list(rep(list(NULL),1)),recursive = F)
      
      bt_cl <- as.numeric(c(all_tested_nodes[clade], parental_node))
      branch_times_clade[1] <- list(bt_cl)
      
      ALL_branch_times_clades[[clade]] <- branch_times_clade
    }
    
    int_nodes <- comb.bck
    # order from present to past
    int_nodes <- names(branching.times(phylo)[order(branching.times(phylo))])[names(branching.times(phylo)[order(branching.times(phylo))]) %in% int_nodes]
    branch_times_to_bck <- rep(list(NULL), length(comb.bck)+1)
    phylo_backbone_cut <- rep(list(NULL), length(comb.bck)+1)
    phylo_backbone_core <- drop.tip(phylo, unlist(ALL_clade_names[comb.sub]))
    
    sb.tips <- rep(list(NULL), length(int_nodes))
    sb.desc <- rep(list(NULL), length(int_nodes))
    names(sb.desc) <- int_nodes
    
    for(sb in 1:length(phylo_backbone_cut)){
      
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
          tips_last_bck <- unlist(ALL_clade_names[comb.sub[!comb.sub %in% unlist(sb.desc, use.names = F)]])
          
          phylo_backbone_cut[[sb]] <- drop.tip(phylo_backbone_core, tips_up_bck)
          names(phylo_backbone_cut)[sb] <- paste(int_nodes[sb-1],"bck", sep = "_")
          
          int_nodes_deep_backbone <- int_nodes[!int_nodes %in% unlist(sapply(branch_times_to_bck, names), use.names = F)]
          
          comb_deep_backbone <- c(comb.sub[!comb.sub %in% unlist(sb.desc, use.names = F)], int_nodes_deep_backbone)
          
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
    
    lin.node <- data.frame(node = c(comb.sub, comb.bck,Ntip(phylo)+1), n.tips = rep(NA, length(comb.sub) + length(comb.bck)+1))
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
        ndata_bck_up <- sum(bck_up$sp_tt_prev)
        
        lin.node_bck$n.tips_prev[l.n] <- lin.node_bck$n.tips[l.n] - ntip_bck_up
        lin.node_bck$sp_tt_prev[l.n] <- lin.node_bck$sp_tt[l.n] - ndata_bck_up  
        
      } else {
        lin.node_bck$n.tips_prev[l.n] <- lin.node_bck$n.tips[l.n]
        lin.node_bck$sp_tt_prev[l.n] <- lin.node_bck$sp_tt[l.n]
      }
      
    }
    lin.node[lin.node$node %in% lin.node_bck$node,] <- lin.node_bck
    
    lin.node <- lin.node[!lin.node$node %in% comb.sub,]
    
    for(j in 1:nrow(best_backbones_df)){
      
      model <- as.character(best_backbones_df$Models[j])
      lamb_pari <- as.numeric(best_backbones_df[j, c("Lambda","Alpha")])
      mu_pari <- as.numeric(best_backbones_df[j, c("Mu","Beta")])
      
      if(backbone.option == "stem.shift"){
        parental_nodes <- Ancestors(phylo, as.numeric(lin.node$node[j]), type = "parent")
        if(parental_nodes == 0){
          parental_nodes <- Ntip(phylo)+1
        }
        agei <- as.numeric(branching.times(phylo)[as.character(parental_nodes)])
      } else { # if branching times contain older branches
        agei <- max(branching.times(phylo)[lin.node$node[j]], max(unlist(branch_times_to_bck[[j]], use.names = F)))
      }
      sizei <- lin.node$sp_tt_prev[j]
      
      time_seq <- c(agei, seq(floor(agei),0,by=-1))
      
      if (grepl("BCST", model)){
        div <-sizei*exp(-abs(lamb_pari[1])*time_seq)
      }
      
      if (grepl("BCST_DCST", model)){
        div<-sizei*exp(-abs(lamb_pari[1])*time_seq+abs(mu_pari[1])*time_seq)
      }
      
      if (grepl("BVAR", model)){
        div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq)))
      }
      
      if (grepl("BVAR_DCST", model)){
        div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq))+abs(mu_pari[1])*time_seq)
      }
      
      if (grepl("BCST_DVAR", model)){
        div<-sizei*exp(-abs(lamb_pari[1])*time_seq-abs(mu_pari[1])/mu_pari[2]*(1-exp(mu_pari[2]*time_seq)))
      }
      
      if (grepl("BVAR_DVAR", model)){
        div<-sizei*exp(abs(lamb_pari[1])/lamb_pari[2]*(1-exp(lamb_pari[2]*time_seq))-abs(mu_pari[1])/mu_pari[2]*(1-exp(mu_pari[2]*time_seq)))
      }
      
      div2 <- div[length(div):1]
      globaldiv[i+j,1:length(div)]<-div2
    }
    if(is.null(comb.bck)){
      row.names(globaldiv) <- c(comb.sub, comb.bck, "backbone")
    } else {
      row.names(globaldiv) <- c(comb.sub, paste("Backbone of", comb.bck),"Deep backbone")
    }
  } else {
    row.names(globaldiv) <- "whole_tree"
  }
  
  globaldiv <- globaldiv[,ncol(globaldiv):1]
  
  if(comb != "whole_tree"){
    past.div.curve<-apply(globaldiv,2,function(x)sum(x,na.rm=T))  
  } else {
    past.div.curve <- globaldiv
  }
  
  if(split.div == T){
    return(globaldiv)
  }else{
    return(past.div.curve)
  }
}
