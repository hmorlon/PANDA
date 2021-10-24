# This function calculates and represent sampling fractions from group given by a data frame (could be the taxonomy, trait data, etc.).
# 
# Created on January 28, 2021
#

get.sampling.fractions <- function(phylo, data, clade.size = 5, plot = F, lad = T, text.cex = 1, pch.cex = 0.8, ...){
  
  # Packages ####
  pkgs <- c("ape", "phytools", "phangorn")
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new.pkgs)){install.packages(new.pkgs)}
  lapply(pkgs, require, character.only = T, quietly = T)
  
  # Checks ####
  if(is.null(phylo) | is(phylo)[1] != "phylo"){
    stop("Object \"phylo\" is NULL or not of class \"phylo\".")
  }
  
  phylo$node.label <- c(Ntip(phylo) + 1):c(Ntip(phylo) + Nnode(phylo))
  
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

  if(any(!phylo$tip.label %in% data$Species)){
    stop("Some species of the phylogeny are not in the database \"data\".")
  }
  
  data[] <- lapply(data, factor)
  
  # Core function ####
  
  phylo.df <- data.frame(nodes = 1:c(Ntip(phylo) + Nnode(phylo)), data = NA, f = NA, sp_in = NA, sp_tt = NA)
  
  data_phylo <- data[data$Species %in% phylo$tip.label, ]
  
  data_loop <- data.frame(data[,!colnames(data) %in% "Species"])
  names(data_loop) <- colnames(data)[!colnames(data) %in% "Species"]
  data_loop[] <- lapply(data_loop, factor)
  
  data_loop_phylo <- data.frame(data_phylo[,!colnames(data_phylo) %in% "Species"])
  names(data_loop_phylo) <- colnames(data)[!colnames(data) %in% "Species"]
  data_loop_phylo[] <- lapply(data_loop_phylo, factor)
  
  for(j in 1:ncol(data_loop_phylo)){
    for(i in 1:nlevels(data_loop_phylo[,j])){
      
      group <- levels(data_loop_phylo[,j])[i]
      
      sp_in <- data_phylo$Species[data_phylo[,c(j+1)] == group]
      sp_tt <- data$Species[data[,c(j+1)] == group]
      
      if(length(sp_in) == 1 | length(sp_tt) == 1){
        next()
      }
      
      MRCA <- getMRCA(phylo, as.character(sp_in))
      
      if(is.na(phylo.df$data[phylo.df$nodes == MRCA]) | phylo.df$data[phylo.df$nodes == MRCA] == "not_mono"){
        if(length(Descendants(phylo, MRCA)[[1]]) != length(sp_in)){
          phylo.df[phylo.df$nodes == MRCA,] <- cbind(MRCA, "not_mono", NA, NA, NA)
          
        } else {
          phylo.df[phylo.df$nodes == MRCA,] <- cbind(MRCA, group, length(sp_in)/length(sp_tt), length(sp_in), length(sp_tt))
          
        }
      }
      
    }
  }
  phylo.df[!colnames(phylo.df) %in% "data"][] <- lapply(phylo.df[!colnames(phylo.df) %in% "data"], as.numeric)
  
  phylo.df1 <- phylo.df[phylo.df$sp_in >= clade.size,]
  phylo.df1 <- na.omit(phylo.df1)
  
  phylo.df$to_test <- ifelse(phylo.df$nodes %in% phylo.df1$nodes, phylo.df$nodes, NA)
  phylo.df$f[!phylo.df$nodes %in% phylo.df$to_test] <- NA
  phylo.df[phylo.df$nodes == Ntip(phylo)+1, c("f", "sp_in", "sp_tt")] <- c(Ntip(phylo)/nrow(data), Ntip(phylo), nrow(data))
  
  # PLOT optionnal ####
  if(plot == T){
    
    if(lad == T){
      pos_leg <- "bottomleft"
      phylo <- ladderize(phylo, right = T)
    } else {
      pos_leg <- "topleft"
      phylo <- ladderize(phylo, F)
    }
    
    phylo$node.label <- c(Ntip(phylo)+1):c(Ntip(phylo)+Nnode(phylo))
    node_legends <- phylo.df$data[phylo.df$nodes %in% phylo$node.label]
    node_legends <- ifelse(node_legends %in% phylo.df$data[phylo.df$nodes %in% phylo.df$to_test], node_legends, NA)
    #node_legends <- ifelse(grepl("clade", node_legends), NA, node_legends)
    
    
    plot(phylo, show.node.label = F, label.offset = 0.4, ...)
    nodelabels(node_legends, adj = c(1.1,-0.5), bg = "white", cex = text.cex, frame = "n")  
    
    axisPhylo(backward = T, cex.axis = text.cex)
    mtext(text = "Time (Myrs)", side = 1, line = 2, at = max(branching.times(phylo))/2, cex = text.cex)
    
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    node <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[node]
    XX_lab <- XX[sort(phylo.df$nodes[phylo.df$nodes > Ntip(phylo) & !is.na(phylo.df$f) & !is.na(phylo.df$to_test) & !is.na(phylo.df$sp_in)]) - Ntip(phylo)]
    YY <- lastPP$yy[node]
    YY_lab <- YY[sort(phylo.df$nodes[phylo.df$nodes > Ntip(phylo) & !is.na(phylo.df$f) & !is.na(phylo.df$to_test) & !is.na(phylo.df$sp_in)]) - Ntip(phylo)]
    BOTHlabels(text="", node, XX_lab, YY_lab, adj = c(0.5, 0.5), 
               frame = "none", pch = 21, thermo = NULL, pie = NULL, 
               piecol = NULL, col = "black", bg = "red", 
               horiz = FALSE, width = NULL, height = NULL, cex= pch.cex)
    
  }
  return(phylo.df)
  
}

