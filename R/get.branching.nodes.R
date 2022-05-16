get.branching.nodes <- function(shift, root_ID = Ntip(phylo)+1){
  
  root_ID <- as.numeric(root_ID)
  
  root_clade <- 0
  root_node <-  NULL
  
  # GET PARENTAL NODES WITH ALL ANCESTORS
  if(length(names(shift)) > 1){
    desc_shift <- Descendants(phylo, as.numeric(names(shift)), type = "all")
    shiftsub <- names(shift)[!sapply(desc_shift, function(x) any(x %in% as.numeric(names(shift))))]
  } else {
    shiftsub <- shift
  }
  
  # account for poor backbone resulting in a subclade
  phylo_backbone_cut <- drop.tip(phylo, unlist(ALL_clade_names[shiftsub]))
  sibling_shift_nodes <- unlist(Siblings(phylo, as.numeric(shiftsub)))
  
  which_sibling_shift_nodes <- sibling_shift_nodes %in% phylo_backbone_cut$node.label[1]
  
  if(any(which_sibling_shift_nodes)){
    sibling_shift_node <- sibling_shift_nodes[which_sibling_shift_nodes]
    sibling_shift_node_bt <- list(list(c(sibling_shift_node, Ancestors(phylo, sibling_shift_node, type = "parent"))))
    names(sibling_shift_node_bt) <- sibling_shift_node
    shift <- c(shift, sibling_shift_node_bt)
  }
  
  desc_shift <- Descendants(phylo, as.numeric(names(shift)), type = "all")
  shiftsub <- names(shift)[!sapply(desc_shift, function(x) any(x %in% as.numeric(names(shift))))]
  
  # account for poor backbone resulting in a subclade
  phylo_backbone_cut <- drop.tip(phylo, unlist(ALL_clade_names[shiftsub]))
  sibling_shift_nodes <- unlist(Siblings(phylo, as.numeric(clade_to_shift)))
  
  which_sibling_shift_nodes <- sibling_shift_nodes %in% phylo_backbone_cut$node.label[1]
  
  if(any(which_sibling_shift_nodes)){
    sibling_shift_node <- sibling_shift_nodes[which_sibling_shift_nodes]
    sibling_shift_node_bt <- list(list(c(sibling_shift_node, Ancestors(phylo, sibling_shift_node, type = "parent"))))
    names(sibling_shift_node_bt) <- sibling_shift_node
    shift <- c(shift, sibling_shift_node_bt)
  }
  
  df_ALL <- as.data.frame(sapply(unlist(shift,recursive = F), function(m) m[2]))
  colnames(df_ALL) <- "node"
  row.names(df_ALL) <- 1:nrow(df_ALL)
  
  # detect the root in the clades TO REMOVE BECAUSE ONLY ON PARENTAL NODES
  #if(any(df_ALL$node == Ntip(phylo) + 1)){
  # root_clade <- 0
  #root_node <-  NULL # because already in df_all
  #} 
  
  df_ALL <- data.frame(node = df_ALL[which(!df_ALL$node %in% c(root_ID)),])
  
  if(nrow(df_ALL) > 1){
    
    all_ancestors <- unlist(list(rep(list(NULL), nrow(df_ALL))),recursive = F)
    
    for(df_l in 1:nrow(df_ALL)){
      
      all_ancestors[df_l] <- list(c(df_ALL$node[df_l],Ancestors(phylo, df_ALL$node[df_l], type = "all")))
      
    }
    
    # removing root node
    all_ancestors <- lapply(all_ancestors, function(x) x[1:c(which(x == root_ID)-1)])
    
    # counting parental nodes
    ALL_par_nodes <- NULL
    
    coal <- as.data.frame(table(sapply(all_ancestors, function(m) m[1])))
    
    while(any(coal$Freq == 2) & is.null(all_ancestors) == F){
      
      if(any(coal$Freq == 2)){
        ALL_par_nodes <- c(ALL_par_nodes,as.numeric(as.character(coal$Var1[coal$Freq == 2])))
        all_ancestors <- unique(all_ancestors)
        
        all_ancestors <- lapply(all_ancestors, function(x) x[x %in% ALL_par_nodes == F])
        
        coal <- as.data.frame(table(sapply(all_ancestors, function(m) m[1])))
         
      } else {
        ALL_par_nodes <- NULL
        all_ancestors <- NULL
      }
    }
    
  } else {ALL_par_nodes <- NULL}
  
  
  if(length(ALL_par_nodes) != 0){
    
    parental_nodes <- unlist(list(rep(list(NULL),length(ALL_par_nodes))),recursive = F)
    
    # ALL OTHER NODES
    if(length(parental_nodes) != 0){
      
      for(df_l in 1:length(parental_nodes)){
        
        if(ALL_par_nodes[df_l] != root_ID){
          
          parental_nodes[[df_l]] <- c(ALL_par_nodes[df_l], Ancestors(phylo, ALL_par_nodes[df_l], type = "parent"))
        }
      }
      # WHETHER PARENTAL NODES ARE THE ROOT
      for(p in 1:length(parental_nodes)){
        
        if(parental_nodes[[p]][2] == root_ID){
          root_node <- parental_nodes[[p]][1]
          root_clade <- 1
        }
        
      }
    }
    
  } else {
    parental_nodes <- NULL
  }
  
  df_ALL <- t(as.data.frame(unlist(shift,recursive = F)))
  
  branches_df_all <- apply(df_ALL, 1, paste, collapse = ".")
  if(!is.null(parental_nodes)){
    branches_parental <- apply(do.call(rbind, parental_nodes), 1, paste, collapse = ".")
    parental_nodes <- parental_nodes[!branches_parental %in% branches_df_all]
  }

  #
  branch_times_to <- unlist(list(rep(list(NULL),nrow(df_ALL) + length(parental_nodes) + root_clade)),recursive = F)
  
  bt_1 <-  unlist(shift,recursive = F)
  for(bt in 1:length((bt_1))){
    branch_times_to[bt] <- bt_1[bt]
  }
  
  p = 0
  if(length(parental_nodes) != 0){
    for(p in 1:length(parental_nodes)){
      branch_times_to[bt + p] <- parental_nodes[p]
    }
  }

  branch_root <- c(Siblings(phylo, root_node),Ancestors(phylo, root_node, type = "parent"))
  
  if(root_clade == 1 & paste(branch_root, collapse = ".") %in% branches_df_all == F){
      branch_times_to[bt + p + root_clade] <- list(branch_root)
  }
  
  #names(branch_times_to) <- c(names(bt_1),rep("parental_node",length(parental_nodes)),rep("root",length(root_node)))
  names(branch_times_to) <- c(names(bt_1), sapply(parental_nodes, function(x) ifelse(!is.null(x), x[1], NULL)), Siblings(phylo, root_node))
  branch_times_to <- branch_times_to[!sapply(branch_times_to, is.null)]
  
  return(branch_times_to)
}
