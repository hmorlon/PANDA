# This function plots a phylogenetic tree and represents shifts of diversification rates detected
# by the function shift.estimates.res
#
# Version 1.0 from November 12, 2020
# geotime.scale not emplemented yet
# 

# Version 1.1 from May 13, 2021
# gts is not included in this version but developped as an external function add.geochrono
# 


plot.phy.comb <- function(phy, data, sampling.fractions, shift.estimates.res, combi = 1,
                          backbone.option = "backbone2",
                          col.sub = NULL, col.bck = "black", lad = T, tested_nodes = F, lty.bck = 1,
                          text.cex = 1, pch.cex = 1,
                          leg = T, ...){
  
  # Loading packages ####
  pkgs <- c("ape", "phytools", "phangorn","picante","stats", "strap", "geiger", "RColorBrewer")
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new.pkgs)){install.packages(new.pkgs)}
  lapply(pkgs, require, character.only = T)
  
  # Checking arguments ####
  # phy
  if(!inherits(phy, "phylo")){
    stop("object \"phy\" is not of class \"phylo\"")
  }
  # data
  if(!inherits(data, "data.frame")){
    stop("object \"data\" is not of class \"data.frame\"")
  }
  
  # sampling.fractions
  if(phy$Nnode + Ntip(phy) != nrow(sampling.fractions) | is(sampling.fractions)[1]!="data.frame"){
    stop("object \"sampling.fractions is not of class \"data.frame\" or is do not correspond to the provided phylogeny")
  }
  # shift.estimates.res
  if(!is(shift.estimates.res)[1] == "list" | any(names(shift.estimates.res) != c("whole_tree", "subclades", "backbones", "total"))){
    stop("object \"shift.estimates.res\" might be incorrect.")
  }
  
  if(!is.numeric(combi)){
    stop("argument \"combi\" should be numeric.")
  }
  
  if(is(tested_nodes)[1] != "logical" ){
    stop("argument \"tested_nodes\" should be logical.")
  }
  
  if(is(leg)[1] != "logical"){
    stop("argument \"leg\" should be logical.")
  }

  if(is(lad)[1] != "logical"){
    stop("argument \"lad\" should be logical.")
  }
  
  if(is(lty.bck)[1] != "numeric"){
    stop("argument \"lty.bck\" should be numeric.")
  }
  
  if(!backbone.option %in% c("backbone1", "backbone2")){
    cat("\nArgument \"backbone.option\" is incorrect.")
    stop()
  }
  
  # Script ####
  
  if(any(names(sampling.fractions) == "taxo")){
    names(sampling.fractions)[names(sampling.fractions) == "taxo"] <- "data"
  }
  
  phy1 <- phy
  
  if(lad == T){
    pos_leg <- "bottomleft"
    phy1 <- ladderize(phy1, right = T)
  } else {
    pos_leg <- "topleft"
    phy1 <- ladderize(phy1, F)
  }
  
  phy1$node.label <- c(Ntip(phy1)+1):c(Ntip(phy1)+Nnode(phy1))
  node_legends <- sampling.fractions$data[sampling.fractions$nodes %in% phy1$node.label]
  node_legends <- ifelse(node_legends %in% sampling.fractions$data[sampling.fractions$nodes %in% sampling.fractions$to_test], node_legends, NA)
  
  comb <- shift.estimates.res$total$Combination[combi]
  
  if(comb == "whole_tree"){
    colors_clades <- rep("black", Nedge(phy1))
    lty_clades <- rep(lty.bck, Nedge(phy1))
    
  } else {
    
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
    
    if(!is.null(col.sub)){
      if(length(col.sub) != length(comb.sub)){
        stop("length of argument \"col.sub\" should match the number of subclades.")
      }
    }
    
    names_leg <- sampling.fractions$data[sampling.fractions$nodes %in% comb.sub]
    for(j in 1:length(names_leg)){
      if(is.na(names_leg[j])){
        names_leg_NA <- sampling.fractions$data[unlist(Descendants(phy1, as.numeric(comb.sub[j]), "children"))]
        if(length(names_leg_NA) > 3 | any(is.na(names_leg_NA))){
          names_leg_NA <- paste0("node ", comb.sub[j])
        }
        names_leg[j] <- paste(names_leg_NA, collapse = " + ")
      }
    }
    
    if(is.null(comb.bck)){
      names_leg <- c(names_leg,"Backbone")
    }else{
      names_leg1 <- c(paste("Backbone from", comb.bck),"Deep backbone")
      names_leg <- c(names_leg,names_leg1)
      # to test with an example
    }
    # Multi combi if possible.
    
    if(is.null(col.sub)){
      col.sub <- c(c(brewer.pal(8, "Dark2"),brewer.pal(8, "Set1"),"darkmagenta","dodgerblue2" , "orange", "forestgreen"))[c(1:length(comb.sub))]
    } 
    
    if(!is.null(comb.bck) & length(col.bck) == 1){
      col.bck <- c(c("blue4", "orange4", "red4", "grey40", "coral4", "deeppink4", "khaki4", "darkolivegreen", "darkslategray")[1:c(length(comb.bck))],"black")
      colors_clades <- rep("black", Nedge(phy1))
    } else {
      colors_clades <- rep(col.bck[length(col.bck)], Nedge(phy1))
    }
    
    if(!is.null(comb.bck)){
      if(length(col.bck) != length(comb.bck)+1){
        stop("length of argument \"col.bck\" should match the number of backbones.")
      }
    }
    
    lty_clades <- rep(lty.bck, Nedge(phy1))
    
    for(i in 1:length(comb.sub)){
      clade_edges <- Descendants(phy1, as.numeric(comb.sub[i]), type = "all")
      if(backbone.option == "backbone1"){
        clade_edges <- c(as.numeric(comb.sub[i]),clade_edges)
      }
      colors_clades[which(phy1$edge[,2] %in% clade_edges)] <- col.sub[i]
      lty_clades[which(phy1$edge[,2] %in% clade_edges)] <- 1
    }
    
    if(!is.null(comb.bck)){
      for(j in 1:length(comb.bck)){
        clade_edges <- Descendants(phy1, as.numeric(comb.bck[j]), type = "all")
        colors_clades[which(phy1$edge[,2] %in% clade_edges)] <- ifelse(colors_clades[which(phy1$edge[,2] %in% clade_edges)] == col.bck[length(col.bck)], col.bck[j], colors_clades[which(phy1$edge[,2] %in% clade_edges)])
      }
    }
    
    model_leg <- sapply(shift.estimates.res$subclades[comb.sub], function(x) x$Models[1])
    
    if(!is.null(comb.bck)){
      model_leg_bck <- sapply(shift.estimates.res$backbones[paste(comb.sub, collapse = ".")][[1]][paste(comb.bck, collapse = ".")][[1]], function(x) x$Models[1])
    } else {
      model_leg_bck <- sapply(shift.estimates.res$backbones[paste(comb.sub, collapse = ".")][[1]][[1]], function(x) x$Models[1])
    }
    model_leg <- c(model_leg, model_leg_bck)
    model_leg <- gsub("_", " ", model_leg)
    model_leg <- paste0(names_leg, " (", model_leg, ")")
    
  }
  
  if(shift.estimates.res$total$delta_AICc[combi] > 0){
    AICc_leg <- paste0("Combination ", combi, " (delta AICc = " , round(shift.estimates.res$total$delta_AICc[combi],3),")")
  } else {
    AICc_leg <- "Best combination"
  }

  
  plot(phy1, label.offset = 0.4, edge.color = colors_clades, edge.lty = lty_clades, ...)
  mtext(text = AICc_leg, line = 1, side = 3, cex = text.cex)

  if(comb == "whole_tree"){
    if(leg == T){
      model_leg <- shift.estimates.res$whole_tree$Models[shift.estimates.res$whole_tree$AICc == min(shift.estimates.res$whole_tree$AICc)]
      legend(pos_leg, legend = paste0("whole_tree (", model_leg, ")"), text.col = "black",cex = text.cex, bty = "n")
    }
    
  } else {
    
    if(tested_nodes == T){
      pos_leg_n <- c(par("xaxp")[1]-2, c(par("yaxp")[2]+3))
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      node <- (lastPP$Ntip + 1):length(lastPP$xx)
      XX <- lastPP$xx[node]
      XX_lab <- XX[sort(sampling.fractions$nodes[sampling.fractions$nodes > Ntip(phy1) & !is.na(sampling.fractions$f) & !is.na(sampling.fractions$to_test) & !is.na(sampling.fractions$sp_in)]) - Ntip(phy1)]
      YY <- lastPP$yy[node]
      YY_lab <- YY[sort(sampling.fractions$nodes[sampling.fractions$nodes > Ntip(phy1) & !is.na(sampling.fractions$f) & !is.na(sampling.fractions$to_test) & !is.na(sampling.fractions$sp_in)]) - Ntip(phy1)]
      BOTHlabels(text="", node, XX_lab, YY_lab, adj = c(0.5, 0.5), 
                 frame = "none", pch = 21, thermo = NULL, pie = NULL, 
                 piecol = NULL, col = "black", bg = "red", 
                 horiz = FALSE, width = NULL, height = NULL, cex=pch.cex)
    }
    # bottomleft
    
    if(leg == T){
      legend(pos_leg, legend = paste(model_leg,sep=" "), text.col = c(col.sub, col.bck), title = "Subgroups (Best models)",
             title.col = "black", cex = text.cex, bty = "n")
    }
  }
  
}

