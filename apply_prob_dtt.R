# This function applies prob_dtt to a combination of shifts from the output of shift.estimates.
# 
# N.B.: last time points is not present because of prob_dtt.
# Two options: (1) modify prob_dtt (2) or this function

apply_prob_dtt <- function(phy, data, sampling.fractions, shifts,
                           combi = 1, backbone.option = "backbone2", scale = 1){
  
  # some checks ####
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
  
  if(any(!phy$tip.label %in% data$Species)){
    cat("The following tips are not in the database.\n \n")
    cat(phy$tip.label[!phy$tip.label %in% data$Species], "\n")
    stop()
  }
  
  if(!backbone.option %in% c("backbone1","backbone2")){
    stop("\"backbone.option\" argument is incorrect." )
  }
  
  if(!is(shifts)[1] == "list" | any(names(shifts) != c("whole_tree", "subclades", "backbones", "total"))){
    stop("object \"shift.estimates.res\" might be incorrect.")
  }
  if(!is.numeric(combi)){
    stop("object \"combi\" should be numeric.")
  }
  
  # because backbone.option = "backbone2" by default
  type <- "crown" 
  
  # functions ####
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
  
  # core script ####
  
  comb <- shifts$total$Combination[combi]
  
  # Deterministic diversity to set the limit
  
  if(comb == "whole_tree"){
    
    whole_df <- shifts$whole_tree[shifts$whole_tree$AICc == min(shifts$whole_tree$AICc),]
    whole_fit.bd <- mimic_fit.bd(whole_df)

    N0 <- sampling.fractions$sp_tt[sampling.fractions$nodes == Ntip(phy)+1]
    l <- sampling.fractions$sp_in[sampling.fractions$nodes == Ntip(phy)+1]    
    
    whole_diversity <- paleodiv(phy = phy, data = data, split.div = F,
                            sampling.fractions = sampling.fractions, shift.estimates.res = shifts, combi = combi)
    prob_whole <- list(prob_dtt(whole_fit.bd, tot_time, 1:tot_time,
                                N0 = N0, l = l, type = type, prec = 10000,
                                m = 1:max(whole_diversity)))
    
    colSums(prob_whole[[1]])
    
    return(prob_whole)
    
  }  else { # combination is not the whole tree
    
    comb <- strsplit(comb, split = "/")[[1]]
    comb <- sapply(comb, function(x) strsplit(x, split = "[.]"))
    comb.sub <- comb[[1]]
    if(length(comb) > 1){
      comb.bck <- comb[[2]]
    } else {
      comb.bck <- NULL
    }
    
    # Deterministic diversity to set the limit
    diversities <- paleodiv(phy = phy, data = data, split.div = T,
                            sampling.fractions = sampling.fractions, shift.estimates.res = shifts, combi = combi)
    row.names(diversities)[!row.names(diversities) %in% comb.sub] <- unlist(ifelse(!is.null(comb.bck), list(c(comb.bck, as.character(Ntip(phy)+1))), Ntip(phy)+1))
    max_diversities <- ceiling(round(apply(diversities, 1, max, na.rm = T))/10)*10
   
    # subclade(s) ####
    tips_sub <- Descendants(phy, as.numeric(comb.sub))
    
    subclades_trees <- list()
    for(i in 1:length(tips_sub)){
      subclades_trees[[i]] <- subtree(phy, phy$tip.label[tips_sub[[i]]])
    }
    
    # fit_bd values for subclade(s) ####
    subclades_df <- list()
    for(i in 1:length(comb.sub)){
      subclades_df[[i]] <- shifts$subclades[names(shifts$subclades) == comb.sub[i]][[1]]
    }
    
    subclades_fit.bd <- lapply(subclades_df, mimic_fit.bd)
    names(subclades_fit.bd) <- names(subclades_trees) <- comb.sub
    
    # can be done without subclades trees
    subclades_tot_times <- sapply(subclades_trees, function(x) max(node.age(x)$ages))
    if(backbone.option == "backbone1"){
      subclades_tot_times <- sapply(subclades_trees, function(x) max(node.age(x)$ages)+max(x$root.edge))
      type <- "stem"
    }
    
    subclades_N0 <- sapply(comb.sub, function(x) sampling.fractions$sp_tt[sampling.fractions$nodes == x])
    subclades_l <- sapply(comb.sub, function(x) sampling.fractions$sp_in[sampling.fractions$nodes == x])
  
  
    prob_subclades <- list()
    for(i in 1:length(subclades_trees)){
      l <- subclades_l[i]
      N0 <- subclades_N0[i]
      method <- ifelse(l/N0 == 1, "simple", "hard")
      max_div <- max_diversities[names(max_diversities) == names(subclades_fit.bd)[i]]
      if(max_div > 10000){
        max_div <- 1000
      }
      prob_subclades[[i]] <- prob_dtt(fit.bd = subclades_fit.bd[[i]], tot_time = subclades_tot_times[i],
                                      time = 1:subclades_tot_times[i],
                                      N0 = N0, l = l, prec = 1000,
                                      m = 1:round(max_div))
    }
    names(prob_subclades) <- comb.sub
    
    lapply(prob_subclades, colSums)
    
    # N0 for backbones ####
    lin.node <- data.frame(node = c(comb.sub, comb.bck,Ntip(phy)+1), n.tips = rep(NA, length(comb.sub) + length(comb.bck)+1))
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
      int.desc_lin <- unlist(Descendants(phy, as.numeric(lin.node_bck$node[l.n]), "all"))
      int.desc_lin <- int.desc_lin[int.desc_lin > Ntip(phy)]
      
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
    
    backbone_tre <- drop.tip(phy, phy$tip.label[unlist(tips_sub)])
    backbone_trees <- list()
    if(!is.null(comb.bck)){
      for(i in 1:length(comb.bck)){
        
        tips_bcki <- phy$tip.label[Descendants(phy, as.numeric(comb.bck))[[1]]]
        tips_bcki <- tips_bcki[tips_bcki %in% backbone_tre$tip.label]
        backbone_trees[[i]] <- subtree(backbone_tre, tips_bcki)
        
        if(i == length(comb.bck)){
          backbone_trees[[i+1]] <- drop.tip(backbone_tre, unlist(lapply(backbone_trees, function(x) x$tip.label)))
        }
      }
      names(backbone_trees) <- c(comb.bck, Ntip(phy)+1)
    } else {
      backbone_trees[[1]] <- backbone_tre
      names(backbone_trees) <- c(Ntip(phy)+1)
    }
    
    
    # fit_bd values for backbone(s) ####
    if(!is.null(comb.bck)){
      backbones_df <- shifts$backbones[paste(comb.sub, collapse = ".")][[1]][paste(comb.bck, collapse = ".")][[1]]
    } else {
      backbones_df <- shifts$backbones[paste(comb.sub, collapse = ".")][[1]][[1]][[1]]
    }
    
    backbone_fit.bd <- list()
    if(!is.null(comb.bck)){
      backbone_fit.bd <- lapply(backbones_df, mimic_fit.bd)
    } else {
      backbone_fit.bd <- list(mimic_fit.bd(backbones_df))
    }
  
    if(!is.null(comb.bck)){
      names(backbone_fit.bd) <- c(comb.bck, Ntip(phy)+1)
    } else {
      names(backbone_fit.bd) <- Ntip(phy)+1
    }
    
    type <- rep(type, nrow(lin.node))
    
    backbone_tot_times <- c()
    for(i in 1:nrow(lin.node)){
      if(backbone.option == "backbone1"){
        parental_nodes <- Ancestors(phy, as.numeric(lin.node$node[i]), type = "parent")
        # last backbone should be "crown"
        if(i == nrow(lin.node)){
          type[i] <- "crown"
        }
        if(parental_nodes == 0){
          parental_nodes <- Ntip(phy)+1
        }
        backbone_tot_times[i] <- as.numeric(branching.times(phy)[as.character(parental_nodes)])
      } else {
        backbone_tot_times[i] <- as.numeric(branching.times(phy)[lin.node$node[i]])
      }
    }
    
    #branch_times <- list(c(branching.times(phy)["89"], branching.times(phy)[as.character(Ancestors(phy, 89))]))
    #spec_times <- NULL
    #backbone <- "backbone2"
    #cond = "crown"
    
    #test_backbone <- fit_bd_backbone(backbone_trees[[2]], tot_time = backbone_tot_times[i],
    #                                 f.lamb = f.lamb, f.mu = f.mu, f = l/N0,
    #                                 lamb_par = c(0.1,0.01), mu_par = c(0.01, 0.001),
    #                                 backbone = backbone, spec_times = spec_times, branch_times = branch_times,
    #                                 cst.lamb = F, cst.mu = F, expo.lamb = T, expo.mu = T,
    #                                 cond = cond, model = "BVAR_DVAR", fix.mu = F)
    
    prob_backbone <- list()
    for(i in 1:length(backbone_fit.bd)){
      l <- lin.node$n.tips_prev[lin.node$node == names(backbone_fit.bd)[i]]
      N0 <- lin.node$sp_tt_prev[lin.node$node == names(backbone_fit.bd)[i]]
      method <- ifelse(l/N0 == 1, "simple", "hard")
      
      max_div <- max_diversities[names(max_diversities) == names(backbone_fit.bd)[i]]
      if(max_div > 10000){
        max_div <- 1000
      }
      prob_backbone[[i]] <- prob_dtt(fit.bd = backbone_fit.bd[[i]], tot_time = backbone_tot_times[i],
                                     time = 1:backbone_tot_times[i], type = type[i], prec = 1000,
                                     method = method, l = l, N0 = N0,
                                     m = 1:round(max_div*scale))
    }
    
    summary(colSums(prob_backbone[[1]]))
    
    as.vector(colSums(prob_backbone[[1]]))
    which(colSums(prob_backbone[[1]]) == min(colSums(prob_backbone[[1]])))
    
    #plot_dtt(test_backbone, backbone_tot_times[1], N0)
    #plot_dtt(backbone_fit.bd[[1]], backbone_tot_times, N0)

    names(prob_backbone) <- names(backbone_fit.bd)
    
    return(list(backbones = prob_backbone, subclades = prob_subclades))
    
  }
}

#mat <- prob_backbone[[1]]
#mean.val <- t(as.numeric(rownames(mat))) %*% mat
#max_mean_div <- max(mean.val)

#plot_prob_dtt(prob_backbone[[1]], lwd = 2, grain = 0.1, ylim = c(0, max_mean_div + max_mean_div*0.50))
#par(new = T)
#plot(1:length(diversities["88",]), diversities["88",], type = "l", ylim = c(0, max_mean_div + max_mean_div*0.50),
#     xlab = "", ylab = "", axes = F, lty = 3, col = "deepskyblue4", lwd = 3)

#plot_prob_dtt(prob_backbone[[1]], lwd = 2, grain = 0.05)
#plot_prob_dtt(prob_subclades[[1]], lwd = 2, grain = 0.05, add = T, col.mean = "blue")




