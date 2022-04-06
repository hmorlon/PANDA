plot.phylo.comb <- function(phylo, data, sampling.fractions, shift.res = NULL,
                            combi, backbone.option = "crown.shift",
                            main = NULL, col.sub = NULL, col.bck = "black",
                            lty.bck = 1, tested_nodes = F, lad = T,
                            leg = T, text.cex = 1, pch.cex = 1, ...){
  
  # Checking arguments ####
  # phylo
  if(!inherits(phylo, "phylo")){
    stop("object \"phylo\" is not of class \"phylo\"")
  }
  # data
  if(!inherits(data, "data.frame")){
    stop("object \"data\" is not of class \"data.frame\"")
  }
  # sampling.fractions
  if(phylo$Nnode + Ntip(phylo) != nrow(sampling.fractions) | is(sampling.fractions)[1]!="data.frame"){
    stop("object \"sampling.fractions is not of class \"data.frame\" or is do not correspond to the provided phylogeny")
  }

  if(!is.null(shift.res)){
    if(!is(shift.res)[1] == "list" | any(names(shift.res) != c("whole_tree", "subclades", "backbones", "total"))){
      stop("object \"shift.res\" might be incorrect.")
    }
    if(!is.numeric(combi)){
      stop("if shift.res is specified, argument \"combi\" should be numeric.")
    }
  } else {
    if(!is.character(combi)){
      stop("if shift.res is not specified, argument \"combi\" should be a character.")
    }
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
  
  if(!backbone.option %in% c("stem.shift", "crown.shift")){
    cat("\nArgument \"backbone.option\" is incorrect.")
    stop()
  }
  
  # Script ####
  
  if(any(names(sampling.fractions) == "taxo")){
    names(sampling.fractions)[names(sampling.fractions) == "taxo"] <- "data"
  }
  
  phylo1 <- phylo
  
  if(lad == T){
    pos_leg <- "bottomleft"
    phylo1 <- ladderize(phylo1, right = T)
  } else {
    pos_leg <- "topleft"
    phylo1 <- ladderize(phylo1, F)
  }
  
  phylo1$node.label <- c(Ntip(phylo1)+1):c(Ntip(phylo1)+Nnode(phylo1))
  node_legends <- sampling.fractions$data[sampling.fractions$nodes %in% phylo1$node.label]
  node_legends <- ifelse(node_legends %in% sampling.fractions$data[sampling.fractions$nodes %in% sampling.fractions$to_test], node_legends, NA)
  
  if(!is.null(shift.res)){
    comb <- shift.res$total$Combination[combi]
  } else {
    comb <- combi
  }
  
  if(comb == "whole_tree"){
    colors_clades <- rep("black", Nedge(phylo1))
    lty_clades <- rep(lty.bck, Nedge(phylo1))
    
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
        names_leg_NA <- sampling.fractions$data[unlist(Descendants(phylo1, as.numeric(comb.sub[j]), "children"))]
        if(length(names_leg_NA) > 3 | any(is.na(names_leg_NA))){
          names_leg_NA <- paste0("node ", comb.sub[j])
        }
        names_leg[j] <- paste(names_leg_NA, collapse = " + ")
      }
    }
    
    if(is.null(comb.bck)){
      names_leg <- c(names_leg,"Backbone")
    }else{
      names_leg1 <- c(paste("Backbone of", sampling.fractions$data[sampling.fractions$nodes %in% comb.bck]),"Deep backbone")
      names_leg <- c(names_leg,names_leg1)
    }
    
    if(is.null(col.sub)){
      col.sub <- c(c(brewer.pal(8, "Dark2"),brewer.pal(8, "Set1"),"darkmagenta","dodgerblue2" , "orange", "forestgreen"))[c(1:length(comb.sub))]
    } 
    
    if(!is.null(comb.bck) & length(col.bck) == 1){
      col.bck <- c(c("blue4", "orange4", "red4", "grey40", "coral4", "deeppink4", "khaki4", "darkolivegreen", "darkslategray")[1:c(length(comb.bck))],"black")
      colors_clades <- rep("black", Nedge(phylo1))
    } else {
      colors_clades <- rep(col.bck[length(col.bck)], Nedge(phylo1))
    }
    
    if(!is.null(comb.bck)){
      if(length(col.bck) != length(comb.bck)+1){
        stop("length of argument \"col.bck\" should match the number of backbones.")
      }
    }
    
    lty_clades <- rep(lty.bck, Nedge(phylo1))
    
    for(i in 1:length(comb.sub)){
      clade_edges <- Descendants(phylo1, as.numeric(comb.sub[i]), type = "all")
      if(backbone.option == "stem.shift"){
        clade_edges <- c(as.numeric(comb.sub[i]),clade_edges)
      }
      colors_clades[which(phylo1$edge[,2] %in% clade_edges)] <- col.sub[i]
      lty_clades[which(phylo1$edge[,2] %in% clade_edges)] <- 1
    }
    
    if(!is.null(comb.bck)){
      for(j in 1:length(comb.bck)){
        clade_edges <- Descendants(phylo1, as.numeric(comb.bck[j]), type = "all")
        colors_clades[which(phylo1$edge[,2] %in% clade_edges)] <- ifelse(colors_clades[which(phylo1$edge[,2] %in% clade_edges)] %in% col.bck, col.bck[j], colors_clades[which(phylo1$edge[,2] %in% clade_edges)])
      }
    }
    
    if(!is.null(shift.res)){
      model_leg <- sapply(shift.res$subclades[comb.sub], function(x) x$Models[1])
      
    if(any(grepl("/", shift.res$total$Combination))){
      model_leg_bck <- unlist(sapply(shift.res$backbones[paste(paste(comb.sub, collapse = "."),paste0(comb.bck, collapse = "."), sep = "/")][[1]], function(x) x$Models[1]))
    }
    
      model_leg <- c(model_leg, model_leg_bck)
      model_leg <- gsub("_", " ", model_leg)
      model_leg <- paste0(names_leg, " (", model_leg, ")")
    } else {
      model_leg <- names_leg
    }
  }
  
  if(is.null(main)){
    if(!is.null(shift.res)){
      if(shift.res$total$delta_AICc[combi] > 0){
        main <- paste0("Combination ", combi, " (delta AICc = " , round(shift.res$total$delta_AICc[combi],3),")")
      } else {
        main <- "Best combination"
      }
    } else {
      main <- ""
    }
  }
  
  plot(phylo1, edge.color = colors_clades,
       edge.lty = lty_clades, main = main, ...)
  
  if(tested_nodes == T){
    pos_leg_n <- c(par("xaxp")[1]-2, c(par("yaxp")[2]+3))
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    node <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[node]
    XX_lab <- XX[sort(sampling.fractions$nodes[sampling.fractions$nodes > Ntip(phylo1) & !is.na(sampling.fractions$f) & !is.na(sampling.fractions$to_test) & !is.na(sampling.fractions$sp_in)]) - Ntip(phylo1)]
    YY <- lastPP$yy[node]
    YY_lab <- YY[sort(sampling.fractions$nodes[sampling.fractions$nodes > Ntip(phylo1) & !is.na(sampling.fractions$f) & !is.na(sampling.fractions$to_test) & !is.na(sampling.fractions$sp_in)]) - Ntip(phylo1)]
    BOTHlabels(text="", node, XX_lab, YY_lab, adj = c(0.5, 0.5), 
               frame = "none", pch = 21, thermo = NULL, pie = NULL, 
               piecol = NULL, col = "black", bg = "red", 
               horiz = FALSE, width = NULL, height = NULL, cex=pch.cex)
  }
  leg_title <- ""
  if(comb == "whole_tree"){
    if(leg == T){
      if(!is.null(shift.res)){
        model_leg <- shift.res$whole_tree$Models[shift.res$whole_tree$AICc == min(shift.res$whole_tree$AICc)]
        legend(pos_leg, legend = paste0("whole_tree (", model_leg, ")"), text.col = "black",cex = text.cex, bty = "n")  
      } else {
        legend(pos_leg, legend = paste0("whole_tree"), text.col = "black",cex = text.cex, bty = "n")  
      }
    } else {
      title <- ""
    }
  } else {
    if(leg == T){
      if(!is.null(shift.res)){
        leg_title <- "Shifts (Best model)"
      } else {
        leg_title <- "Shifts"
      }
      
      legend(pos_leg, legend = paste(model_leg,sep=" "), text.col = c(col.sub, col.bck), title = leg_title,
             title.col = "black", cex = text.cex, bty = "n")
    }
  }
  
}

