paleodiv <- function(phylo, data, sampling.fractions, shift.res,
                     backbone.option = "backbone2", combi = 1, split.div = F){

  # Checking arguments ####
  # phylo
  if(!inherits(phylo, "phylo")){
    stop("object \"phylo\" is not of class \"phylo\"")
  }
  
  phylo$node.label <- NULL
  
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
  
  if(!backbone.option %in% c("backbone1", "backbone2")){
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
  
  #time.seq <- c(tot_time, seq(floor(tot_time),0,by=-1))
  time.seq <- unlist(ifelse(round(tot_time) == floor(tot_time), list(seq(floor(tot_time),0,by=-1)), list(c(tot_time, seq(floor(tot_time),0,by=-1)))))
  globaldiv <- matrix(NA,length(comb.sub)+length(comb.bck)+ifelse(comb.sub != "whole_tree",1,0), length(time.seq)) #matrix(NA, Number of clades, Crown age of the whole tree + 1)
  
  if(any(comb.sub == "whole_tree")){
    
    best_subclades_df_combi <- shift.res$whole_tree
    best_subclades_df_combi <- best_subclades_df_combi[best_subclades_df_combi$AICc == min(best_subclades_df_combi$AICc),]
    best_subclades_df_combi$Clades <- comb.sub
    
    tot_time2 <- tot_time
    names(tot_time2) <- comb.sub
    totalsp2 <- totalsp
    names(totalsp2) <- comb.sub
    
  } else {
    
    best_subclades_df_combi <- best_subclades_df[best_subclades_df$Clades %in% as.numeric(comb.sub),]
    best_subclades_df_combi <- best_subclades_df_combi[match(comb.sub, best_subclades_df_combi$Clades), ]
    
    if(backbone.option == "backbone1"){
      parental_nodes <- Ancestors(phylo, as.numeric(best_subclades_df_combi$Clades), type = "parent")
      tot_time2 <- as.list(branching.times(phylo)[as.character(parental_nodes)])
    } else {
      tot_time2 <- as.list(branching.times(phylo)[best_subclades_df_combi$Clades])
    }
    
    totalsp2 <- as.list(sampling.fractions$sp_tt[sampling.fractions$nodes %in% as.numeric(comb.sub)])
    names(totalsp2) <- comb.sub
  }
  
  # Subclades diversity (RPANDA FONCTIONS !!!!)
  
  for(i in 1:nrow(best_subclades_df_combi)){
    clade <- as.character(best_subclades_df_combi$Clades[i])
    model <- as.character(best_subclades_df_combi$Models[i])
    values <- as.numeric(best_subclades_df_combi[i,-c(1,2,9)])
    names(values) <- names(best_subclades_df_combi[i,-c(1,2,9)])
    
    lamb_pari <- as.numeric(c(values["Lambda"],values["Alpha"]))
    mu_pari <- as.numeric(c(values["Mu"],values["Beta"]))
    agei <- tot_time2[[clade]]
    sizei <- totalsp2[[clade]]
    #time_seq <- c(agei, seq(floor(agei),0,by=-1))
    time_seq <- unlist(ifelse(round(agei) == floor(agei), list(seq(floor(agei),0,by=-1)), list(c(agei, seq(floor(agei),0,by=-1)))))
    
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
  
  if(all(comb.sub != "whole_tree")){
    # Backbone diversity
    
    if(is.null(comb.bck)){
      best_backbones <- shift.res$backbones[paste(comb.sub, collapse = ".")][[1]][[1]][[1]]
      best_backbones_df <- best_backbones[1,]
      best_backbones_df$parts <- paste0(paste(comb.sub, collapse = "."), "_bck")
      row.names(best_backbones_df) <- NULL
      best_backbones_df <- best_backbones_df[,c(10,1:8)]
      
    } else {
      best_backbones <- shift.res$backbones[paste(comb.sub, collapse = ".")][[1]][paste(comb.bck, collapse = ".")][[1]]
      #best_backbones <- best_backbones[match(comb.bck, names(best_backbones))]
      #best_backbones <- best_backbones[[1]]
      best_backbones_df <- do.call(rbind.data.frame, lapply(best_backbones, function(x) x[1,]))
      best_backbones_df$parts <- row.names(best_backbones_df)
      row.names(best_backbones_df) <- NULL
      best_backbones_df <- best_backbones_df[,c(10,1:8)]
      
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
      
      if(backbone.option == "backbone1"){
        parental_nodes <- Ancestors(phylo, as.numeric(lin.node$node[j]), type = "parent")
        if(parental_nodes == 0){
          parental_nodes <- Ntip(phylo)+1
        }
        agei <- as.numeric(branching.times(phylo)[as.character(parental_nodes)])
      } else {
        agei <- as.numeric(branching.times(phylo)[lin.node$node[j]])
      }
      sizei <- lin.node$sp_tt_prev[j]
      
      #time_seq <- c(agei, seq(floor(agei),0,by=-1))
      
      time_seq <- unlist(ifelse(round(agei) == floor(agei), list(seq(floor(agei),0,by=-1)), list(c(agei, seq(floor(agei),0,by=-1)))))
      
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
      row.names(globaldiv) <- c(comb.sub, paste("Backbone from", comb.bck),"Deep backbone")
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

