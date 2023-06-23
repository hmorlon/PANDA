apply_prob_dtt <- function(phylo, data, sampling.fractions, shift.res,
                           combi = 1, backbone.option = "crown.shift",
                           time.interval = 1, m = NULL){
  
  # some checks ####
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
  
  if(any(!phylo$tip.label %in% data$Species)){
    cat("The following tips are not in the database.\n \n")
    cat(phylo$tip.label[!phylo$tip.label %in% data$Species], "\n")
    stop()
  }
  
  if(!backbone.option %in% c("stem.shift","crown.shift")){
    stop("\"backbone.option\" argument is incorrect." )
  }
  
  if(!is(shift.res)[1] == "list" | any(names(shift.res) != c("whole_tree", "subclades", "backbones", "total"))){
    stop("object \"shift.res\" might be incorrect.")
  }
  if(!is.numeric(combi)){
    stop("object \"combi\" should be numeric.")
  }
  
  # because backbone.option = "crown.shift" by default
  type <- "crown" 
  
  # function ####
  mimic_fit.bd <- function(res){
    
    mod <- res$Models[1]
    mod <- strsplit(mod, split = "_")
    
    if(grepl("BCST", mod[[1]][1])){
      f.lamb <- function(x,y){y}
      cst.lamb = T
      expo.lamb = F
      
    } else {
      f.lamb <- function(x,y){y[1]*exp(y[2]*x)}
      cst.lamb = F
      expo.lamb = T
    }
    
    lamb_par <- res[1,c("Lambda", "Alpha")]
    lamb_par <- lamb_par[!is.na(lamb_par)]
    
    f.lamb.par <- function(t) {
      abs(f.lamb(t, lamb_par))
    }
    
    if(length(mod[[1]]) > 1){
      if(grepl("DCST", mod[[1]][2])){
        f.mu <- function(x,y){y}
        cst.mu = T
        expo.mu = F
      } else {
        f.mu <- function(x,y){y[1]*exp(y[2]*x)}
        cst.mu = F
        expo.mu = T
      }
      
      mu_par <- res[1,c("Mu", "Beta")]
      mu_par <- mu_par[!is.na(mu_par)]
      
      f.mu.par <- function(t) {
        abs(f.mu(t, mu_par))
      }
    } else {
      f.mu <- function(x,y){0}
      cst.mu = T
    }
    
    if(length(mod[[1]]) == 1){
      fit.bd <- list(model = "birth.death",  LH = res$logL[1], aicc = res$AICc[1],
                     lamb_par = lamb_par, f.lamb = Vectorize(f.lamb.par))
    } else {
      fit.bd <- list(model = "birth.death",  LH = res$logL[1], aicc = res$AICc[1],
                     lamb_par = lamb_par, mu_par = mu_par,
                     f.lamb = Vectorize(f.lamb.par), f.mu = Vectorize(f.mu.par))
    }
    class(fit.bd) <- "fit.bd"
    return(fit.bd)
    
  }
  no_decline <- function(x){
    dif <- c()
    for(i in 2:length(x)-1){
      dif <- c(dif, x[i+1] - x[i] > 0)
    }
    return(all(dif))
  }
  
  # core script ####
  comb <- shift.res$total$Combination[combi]

  if(comb == "whole_tree"){
    
    whole_df <- shift.res$whole_tree[shift.res$whole_tree$AICc == min(shift.res$whole_tree$AICc),]
    whole_fit.bd <- mimic_fit.bd(whole_df)
    tot_time <- max(branching.times(phylo))

    N0 <- sampling.fractions$sp_tt[sampling.fractions$nodes == Ntip(phylo)+1]
    l <- sampling.fractions$sp_in[sampling.fractions$nodes == Ntip(phylo)+1]    
    
    whole_diversity <- paleodiv(phylo = phylo, data = data, split.div = F,
                                sampling.fractions = sampling.fractions,
                                shift.res = shift.res, combi = combi)
    min_sumprob <- c()
    check_prob <- F
    if(no_decline(whole_diversity)){
      m_range <- 1
    } else{
      m_range <- c(2, 3, 5, 7, 10) 
    }
    
    while(check_prob == F){
      
      if(no_decline(whole_diversity)){
        prob_whole <- list(prob_dtt(whole_fit.bd, tot_time, seq(1,tot_time, time.interval),
                                    N0 = N0, l = l, type = type, prec = 10000,
                                    m = 1:max(whole_diversity)*m_range))
      } else {
        
        prob_whole <- list(prob_dtt(whole_fit.bd, tot_time, seq(1,tot_time, time.interval),
                                    N0 = N0, l = l, type = type, prec = 10000,
                                    m = 1:max(whole_diversity)*m_range[1]))
      }
      # adding diversities at present
      prob_whole[[1]] <- cbind(prob_whole[[1]],rep(0, nrow(prob_whole[[1]])))
      prob_whole[[1]][row.names(prob_whole[[1]]) == as.character(whole_diversity[length(whole_diversity)]),ncol(prob_whole[[1]])] <- 1
      colnames(prob_whole[[1]])[ncol(prob_whole[[1]])] <- "0"
      
      m_range <- m_range[-1]
      
      min_sumprob <- c(min_sumprob, min(colSums(prob_whole[[1]])))
      cat(" -> minimum value of the sum of probabilities/Myr=", min_sumprob[length(min_sumprob)], "\n")
      
      if(min_sumprob[length(min_sumprob)] >= 0.95){
        check_prob <- T
      
      } else {
        if(length(m_range) == 0){
          check_prob <- T
          cat("\nWarnings: the sum of probabilities for each time point did not reach 95%.
                You should use another range of m for the backbone.")
        }
      }
      if(length(min_sumprob) > 1){
        if(min_sumprob[length(min_sumprob)] == min_sumprob[length(min_sumprob)-1]){
          check_prob <- T
          cat("\nWarnings: the sum of probabilities did not reach 95% for each time Myr.\n")
        }
      }
    }
    
    return(prob_whole)
    
  } else { # combination is not the whole tree
    
    comb <- unlist(strsplit(comb, split = "/"))
    comb.sub <- unlist(strsplit(comb[1], split = "[.]"))
    
    if(length(comb) > 1){
      comb.bck <- unlist(strsplit(comb[2], split = "[.]"))
    } else {
      comb.bck <- NULL
    }
    
    # Deterministic diversity to set the limit
    diversities <- paleodiv(phylo = phylo, data = data, split.div = T,
                            sampling.fractions = sampling.fractions, shift.res = shift.res, combi = combi)
    
    if(is.null(m)){
      row.names(diversities)[!row.names(diversities) %in% comb.sub] <- unlist(ifelse(!is.null(comb.bck), list(c(comb.bck, as.character(Ntip(phylo)+1))), Ntip(phylo)+1))
      max_diversities <- ceiling(round(apply(diversities, 1, max, na.rm = T))/10)*10
    } else {
      
      row.names(diversities)[!row.names(diversities) %in% comb.sub] <- unlist(ifelse(!is.null(comb.bck), list(c(comb.bck, as.character(Ntip(phylo)+1))), Ntip(phylo)+1))
      max_diversities <- m
      names(max_diversities) <- row.names(diversities)
      cat("\n#### Maximum values of m are manually specified: ####\n m =", paste0(m[-length(m)], ","), m[length(m)], "\n")
    }
    
    # subclade(s) ####
    tips_sub <- Descendants(phylo, as.numeric(comb.sub))
    
    subclades_trees <- list()
    for(i in 1:length(tips_sub)){
      subclades_trees[[i]] <- subtree(phylo, phylo$tip.label[tips_sub[[i]]])
    }
    
    # fit_bd values for subclade(s) ####
    subclades_df <- list()
    for(i in 1:length(comb.sub)){
      subclades_df[[i]] <- shift.res$subclades[names(shift.res$subclades) == comb.sub[i]][[1]]
    }
    
    subclades_fit.bd <- lapply(subclades_df, mimic_fit.bd)
    names(subclades_fit.bd) <- names(subclades_trees) <- comb.sub
    
    # can be done without subclades trees
    subclades_tot_times <- sapply(subclades_trees, function(x) max(branching.times(x)))
    if(backbone.option == "stem.shift"){
      subclades_tot_times <- sapply(subclades_trees, function(x) max(branching.times(x))+max(x$root.edge))
      type <- "stem"
    }
    
    subclades_N0 <- sapply(comb.sub, function(x) sampling.fractions$sp_tt[sampling.fractions$nodes == x])
    subclades_l <- sapply(comb.sub, function(x) sampling.fractions$sp_in[sampling.fractions$nodes == x])
    
    cat("\nDTT calculation for subclade(s):\n")
    prob_subclades <- list()
    for(i in 1:length(subclades_trees)){
      cat("\t", i, "/", length(subclades_trees), "\n")
      l <- subclades_l[i]
      N0 <- subclades_N0[i]
      method <- ifelse(l/N0 == 1, "simple", "hard")
      max_div <- max_diversities[names(max_diversities) == names(subclades_fit.bd)[i]]
      if(length(seq(1,subclades_tot_times[i], time.interval)) != 1){
        time <- seq(1,subclades_tot_times[i], time.interval)
      } else {
        time <- c(seq(1,subclades_tot_times[i], time.interval),subclades_tot_times[i]) 
      }
      prob_subclades[[i]] <- prob_dtt(fit.bd = subclades_fit.bd[[i]], tot_time = subclades_tot_times[i],
                                      time = time,
                                      N0 = N0, l = l, prec = 1000,
                                      m = 1:round(max_div))
      # adding diversities at present
      prob_subclades[[i]] <- cbind(prob_subclades[[i]],rep(0, nrow(prob_subclades[[i]])))
      prob_subclades[[i]][row.names(prob_subclades[[i]]) == as.character(diversities[i,ncol(diversities)]),
                          ncol(prob_subclades[[i]])] <- 1
      colnames(prob_subclades[[i]])[ncol(prob_subclades[[i]])] <- "0"
    }
    names(prob_subclades) <- comb.sub
    
    # N0 for backbones ####
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
    lin.node <- lin.node[-(1:length(comb.sub)),]
    
    # backbone trees
    
    backbone_tre <- drop.tip(phylo, phylo$tip.label[unlist(tips_sub)])
    backbone_trees <- list()
    if(!is.null(comb.bck)){
      for(i in 1:length(comb.bck)){
        
        tips_bcki <- phylo$tip.label[Descendants(phylo, as.numeric(comb.bck))[[1]]]
        tips_bcki <- tips_bcki[tips_bcki %in% backbone_tre$tip.label]
        backbone_trees[[i]] <- subtree(backbone_tre, tips_bcki)
        
        if(i == length(comb.bck)){
          backbone_trees[[i+1]] <- drop.tip(backbone_tre, unlist(lapply(backbone_trees, function(x) x$tip.label)))
        }
      }
      names(backbone_trees) <- c(comb.bck, Ntip(phylo)+1)
    } else {
      backbone_trees[[1]] <- backbone_tre
      names(backbone_trees) <- c(Ntip(phylo)+1)
    }
    
    
    # fit_bd values for backbone(s) ####
    backbones_df <- shift.res$backbones[paste(paste(comb.sub, collapse = "."), paste(comb.bck, collapse = "."), sep = "/")][[1]]
    backbone_fit.bd <- lapply(backbones_df, mimic_fit.bd)
    
    names(backbone_fit.bd) <- names(max_diversities[!names(max_diversities) %in% comb.sub])
    type <- rep(type, nrow(lin.node_bck))
    
    backbone_tot_times <- c()
    for(i in 1:nrow(lin.node_bck)){
      if(backbone.option == "stem.shift"){
        parental_nodes <- Ancestors(phylo, as.numeric(lin.node_bck$node[i]), type = "parent")
        # last backbone should be "crown"
        if(i == nrow(lin.node_bck)){
          type[i] <- "crown"
        }
        if(parental_nodes == 0){
          parental_nodes <- Ntip(phylo)+1
        }
        backbone_tot_times[i] <- as.numeric(branching.times(phylo)[as.character(parental_nodes)])
      } else {
        backbone_tot_times[i] <- as.numeric(branching.times(phylo)[lin.node_bck$node[i]])
      }
    }
    
    #branch_times <- list(c(branching.times(phylo)["89"], branching.times(phylo)[as.character(Ancestors(phylo, 89))]))
    #spec_times <- NULL
    #backbone <- "crown.shift"
    #cond = "crown"
    
    #test_backbone <- fit_bd_backbone(backbone_trees[[2]], tot_time = backbone_tot_times[i],
    #                                 f.lamb = f.lamb, f.mu = f.mu, f = l/N0,
    #                                 lamb_par = c(0.1,0.01), mu_par = c(0.01, 0.001),
    #                                 backbone = backbone, spec_times = spec_times, branch_times = branch_times,
    #                                 cst.lamb = F, cst.mu = F, expo.lamb = T, expo.mu = T,
    #                                 cond = cond, model = "BVAR_DVAR", fix.mu = F)
    
    cat("\nDTT calculation for backbone(s):\n")
    
    prob_backbone <- list()
    for(i in 1:length(backbone_fit.bd)){
      # first attempt to get a minimum of 95% for the sum of the probabilities per Myr
      min_sumprob <- c()
      check_prob <- F
      cat("\t", i, "/", length(backbone_fit.bd), "\n")
      if(is.null(m)){
        m_range <- c(2, 3, 5, 7, 10) 
      } else{
        m_range <- 1
      }
      
      l <- lin.node_bck$n.tips_prev[lin.node_bck$node == names(backbone_fit.bd)[i]]
      N0 <- lin.node_bck$sp_tt_prev[lin.node_bck$node == names(backbone_fit.bd)[i]]
      method <- ifelse(l/N0 == 1, "simple", "hard")
      max_div <- max_diversities[names(max_diversities) == names(backbone_fit.bd)[i]]
      
      while(check_prob == F){
        cat("with a maximum value of m =", max_div*m_range[1], paste0("(max. deterministic value x",m_range[1],")"), "\n")
        prob_backbone[[i]] <- prob_dtt(fit.bd = backbone_fit.bd[[i]], tot_time = backbone_tot_times[i],
                                       time = seq(1,backbone_tot_times[i], time.interval),
                                       type = type[i], prec = 1000,
                                       method = method, l = l, N0 = N0,
                                       m = 1:round(max_div*m_range[1]))
        m_range <- m_range[-1]
        
        # adding diversities at present
        prob_backbone[[i]] <- cbind(prob_backbone[[i]],rep(0, nrow(prob_backbone[[i]])))
        prob_backbone[[i]][row.names(prob_backbone[[i]]) == as.character(diversities[length(prob_subclades)+i,ncol(diversities)]),
                            ncol(prob_backbone[[i]])] <- 1
        colnames(prob_backbone[[i]])[ncol(prob_backbone[[i]])] <- "0"
        min_sumprob <- c(min_sumprob, min(colSums(prob_backbone[[i]])))
        
        if(min_sumprob[length(min_sumprob)] >= 0.95){
          check_prob <- T
          cat(" -> minimum value of the sum of probabilities/Myr=", min_sumprob[length(min_sumprob)], "\n")
        } else {
          if(length(m_range) == 0){
            check_prob <- T
            cat("\nWarnings: the sum of probabilities for each time point did not reach 95%.\nYou should use another range of m for the backbone.")
          }
        }
        
        if(length(min_sumprob) > 1){
          if(min_sumprob[length(min_sumprob)] == min_sumprob[length(min_sumprob)-1]){
            check_prob <- T
            cat("\nWarnings: the sum of probabilities for each time point did not reach 95%.\nThe minimum value reached is m =", paste0(as.character(round(min_sumprob[length(min_sumprob)]*100, 2)), "%."))
          }
        }
      }
    }
    names(prob_backbone) <- names(backbone_fit.bd)
    
    return(list(backbones = prob_backbone, subclades = prob_subclades))
  }
}
