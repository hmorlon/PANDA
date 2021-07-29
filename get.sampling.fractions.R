# This function calculates and represent sampling fractions from group given by a data frame (could be the taxonomy, trait data, etc.).
# 
# Created on January 28, 2021
#

get.sampling.fractions <- function(phy, data, clade.size = 5, plot = F, lad = T, text.cex = 1, pch.cex = 0.8, ...){
  
  # Packages ####
  pkgs <- c("ape", "phytools", "phangorn")
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new.pkgs)){install.packages(new.pkgs)}
  lapply(pkgs, require, character.only = T, quietly = T)
  
  # Checks ####
  if(is.null(phy) | is(phy)[1] != "phylo"){
    stop("Object \"phy\" is NULL or not of class \"phylo\".")
  }
  
  phy$node.label <- c(Ntip(phy) + 1):c(Ntip(phy) + Nnode(phy))
  
  if(is.null(data) | is(data)[1] != "data.frame"){
    stop("Object \"data\" is NULL or not of class \"data.frame\".")
  }
  
  if(any("Species" %in% colnames(data)) == F){
    stop("No column named \"Species\" in the database
         \nPlease rename the corresponding column with the name \"Species\".")
  }
  
  if(is.null(data$Species)){
    stop("You need to name the column with species names as \"Species\".")
  }
  
  if(any(table(data$Species) > 1)){
    stop("The object data contains some species more than once.
         \nPlease check your database.")
  }
  
  # Species should be the first column
  data <- data[,c("Species", names(data)[names(data) != "Species"])]

  if(any(!phy$tip.label %in% data$Species)){
    stop("Some species of the phylogeny are not in the database \"data\".")
  }
  
  data[] <- lapply(data, factor)
  
  # Core function ####
  
  phy.df <- data.frame(nodes = 1:c(Ntip(phy) + Nnode(phy)), data = NA, f = NA, sp_in = NA, sp_tt = NA)
  
  data_phylo <- data[data$Species %in% phy$tip.label, ]
  
  data_loop <- data.frame(data[,!colnames(data) %in% "Species"])
  names(data_loop) <- colnames(data)[!colnames(data) %in% "Species"]
  data_loop[] <- lapply(data_loop, factor)
  
  data_loop_phylo <- data.frame(data_phylo[,!colnames(data_phylo) %in% "Species"])
  names(data_loop) <- colnames(data)[!colnames(data) %in% "Species"]
  data_loop_phylo[] <- lapply(data_loop_phylo, factor)
  
  for(j in 1:ncol(data_loop)){
    for(i in 1:nlevels(data_loop[,j])){
      
      group <- levels(data_loop[,j])[i]
      
      sp_in <- data_phylo$Species[data_phylo[,c(j+1)] == group]
      sp_tt <- data$Species[data[,c(j+1)] == group]
      
      if(length(sp_in) == 1 | length(sp_tt) == 1){
        next()
      }
      
      MRCA <- getMRCA(phy, as.character(sp_in))
      
      if(is.na(phy.df$data[phy.df$nodes == MRCA]) | phy.df$data[phy.df$nodes == MRCA] == "not_mono"){
        if(length(Descendants(phy, MRCA)[[1]]) != length(sp_in)){
          phy.df[phy.df$nodes == MRCA,] <- cbind(MRCA, "not_mono", NA, NA, NA)
          
        } else {
          phy.df[phy.df$nodes == MRCA,] <- cbind(MRCA, group, length(sp_in)/length(sp_tt), length(sp_in), length(sp_tt))
          
        }
      }
      
    }
  }
  phy.df[!colnames(phy.df) %in% "data"][] <- lapply(phy.df[!colnames(phy.df) %in% "data"], as.numeric)
  
  phy.df1 <- phy.df[phy.df$sp_in >= clade.size,]
  phy.df1 <- na.omit(phy.df1)
  
  phy.df$to_test <- ifelse(phy.df$nodes %in% phy.df1$nodes, phy.df$nodes, NA)
  phy.df$f[!phy.df$nodes %in% phy.df$to_test] <- NA
  phy.df[phy.df$nodes == Ntip(phy)+1, c("f", "sp_in", "sp_tt")] <- c(Ntip(phy)/nrow(data), Ntip(phy), nrow(data))
  
  # PLOT optionnal ####
  if(plot == T){
    
    if(lad == T){
      pos_leg <- "bottomleft"
      phy <- ladderize(phy, right = T)
    } else {
      pos_leg <- "topleft"
      phy <- ladderize(phy, F)
    }
    
    phy$node.label <- c(Ntip(phy)+1):c(Ntip(phy)+Nnode(phy))
    node_legends <- phy.df$data[phy.df$nodes %in% phy$node.label]
    node_legends <- ifelse(node_legends %in% phy.df$data[phy.df$nodes %in% phy.df$to_test], node_legends, NA)
    #node_legends <- ifelse(grepl("clade", node_legends), NA, node_legends)
    
    
    plot(phy, show.node.label = F, label.offset = 0.4, ...)
    nodelabels(node_legends, adj = c(1.1,-0.5), bg = "white", cex = text.cex, frame = "n")  
    
    axisPhylo(backward = T, cex.axis = text.cex)
    mtext(text = "Time (Myrs)", side = 1, line = 2, at = max(branching.times(phy))/2, cex = text.cex)
    
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    node <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[node]
    XX_lab <- XX[sort(phy.df$nodes[phy.df$nodes > Ntip(phy) & !is.na(phy.df$f) & !is.na(phy.df$to_test) & !is.na(phy.df$sp_in)]) - Ntip(phy)]
    YY <- lastPP$yy[node]
    YY_lab <- YY[sort(phy.df$nodes[phy.df$nodes > Ntip(phy) & !is.na(phy.df$f) & !is.na(phy.df$to_test) & !is.na(phy.df$sp_in)]) - Ntip(phy)]
    BOTHlabels(text="", node, XX_lab, YY_lab, adj = c(0.5, 0.5), 
               frame = "none", pch = 21, thermo = NULL, pie = NULL, 
               piecol = NULL, col = "black", bg = "red", 
               horiz = FALSE, width = NULL, height = NULL, cex= pch.cex)
    
  }
  return(phy.df)
  
}

