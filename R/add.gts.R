add.gts <- function(thickness, quaternary = T, is.phylo = F,
                    xpd.x = T, time.interval = 1,
                    names = NULL, fill = T, cex = 1, direction = "rightwards", padj = -0.5){
  
  # BETA VERSION: SHOULD BE TESTED MORE DEEPLY
  par(xpd = T)
  plot_dim = par("usr")
  
  if(is.phylo){
    plot.obj.phylo<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    present <- plot.obj.phylo$xx[1]
    root.age <- 0
  } else {
    root.age <- plot_dim[1] + plot_dim[2]
    present <- 0
  }

  # GTS
  ages <- data.frame(start = NA, end = NA)
  ages[1,] <- -c(0, 0.0117)
  ages[2,] <- -c(0.0117, 2.58)
  ages[3,] <- -c(2.58, 5.33)
  ages[4,] <- -c(5.33, 23.03)
  ages[5,] <- -c(23.03, 33.9)
  ages[6,] <- -c(33.9, 56)
  ages[7,] <- -c(56, 66)
  ages[8,] <- -c(66, 100.5)
  ages[9,] <- -c(100.5, 145)
  ages[10,] <- -c(145, 163.5)
  ages[11,] <- -c(163.5, 174.1)
  ages[12,] <- -c(174.1, 201.3)
  
  ages[13,] <- -c(201.3,251.9)
  
  ages[14,] <- -c(251.9, 298.9)
  ages[15,] <- -c(298.9, 358.9)
  
  if(direction == "leftwards"){
    stop("NOT IMPLEMENTED YET.")
    ages[,c(1,2)] <- apply(ages[,c(1,2)], 2, function(x) -x)
  }
  
  if(quaternary == T){
    ages[2,] <- c(ages[1,1], ages[2,2])
    ages <- ages[-1,]
    col <- c("#FFF1C4", "#FFF7B2", "#FFED00", "#FBCC98", "#FAC18A",
             "#F8B77D", "#BAD25F", "#A0C96D", "#abe1fa", "#71cfeb", "#00b4eb",
             "plum2", "orangered1", "mediumaquamarine") # to check for exact color
    if(is.null(names)){
      ages$names <- c("Q.", "Pli.", "Miocene", "Oligocene", "Eocene", "Paleocene",
                      "Upper Cretaceous", "Lower Cretaceous", "Upper Jurassic", "Middle Jurassic", "Lower Jurassic",
                      "Triassic", "Permian", "Carboniferous")
      #ages$names <- c("", "", "Mio.", "Oli.", "Eo.", "P.",
      #                "Up.Cret.", "Low.Cret.", "Up.Jur.", "MJ.", "Low.Jur.",
      #                "Trias.", "Perm.", "Carb.")
      ages$col <- col
    } 
    
  } else {
    col <- c("#FEF6F2", "#FFF1C4", "#FFF7B2", "#FFED00", "#FBCC98", "#FAC18A",
             "#F8B77D", "#BAD25F", "#A0C96D", "#abe1fa", "#71cfeb", "#00b4eb",
             "plum2", "orangered1", "mediumaquamarine") # to check for exact color
    if(is.null(names)){
      ages$names <- c("Holocene", "Pleistocene", "Pliocene", "Miocene", "Oligocene", "Eocene", "Paleocene",
                      "Upper Cretaceous", "Lower Cretaceous", "Upper Jurassic", "Middle Jurassic", "Lower Jurassic",
                      "Triassic", "Permian", "Carboniferous")
      ages$col <- col
      
    }
  }
  
  thickness <- plot_dim[3] - abs(thickness)
  #thickness <- thickness
  
  # dealing with x
  y2 <- plot_dim[3]
  
  if(is.phylo){
    time.seq <- c(root.age, seq(present - floor(present), present))
  } else {
    if(root.age == 0){
      time.seq <- c(seq(root.age, present))
      present <- time.seq[length(time.seq)]
    } else {
      if(root.age > 0){
        time.seq <- unique(c(root.age, floor(root.age):present))
      } else {
        time.seq <- c(root.age, ceiling(root.age):present)
      }
    }
  }
  
  if(all(time.seq >= 0)){
    if(time.seq[1] > time.seq[2]){
      ages$start <- -ages$start
      ages$end <- -ages$end
    }
  }
    
  ages$start  <- ages$start + present
  ages$end  <- ages$end + present  

  if(direction == "leftwards"){
    stop("NOT IMPLEMENTED YET.")
    xlimit <- ceiling(c(max(plot.obj.phylo$xx)+present)/5)*5
    ages <- ages[which(ages$start < xlimit),]
  } else {
    if(all(time.seq >= 0)){
      if(time.seq[1] > time.seq[2]){
        ages <- ages[c(1:which(ages$end > root.age)[1]),]
      } else {
        ages <- ages[c(1:which(ages$end < root.age)[1]),]
      }
    } else {
      if(root.age > 0){
        ages <- ages[c(1:which(ages$end<0)[1]),]  
      } else {
        ages <- ages[which(ages$start > root.age),]  
      }
    }
  }
  
  # working for rightwards only
  if(xpd.x == F){
    ages$end[nrow(ages)] <- plot_dim[1]
  }
  
  if(!is.null(names)){
    if(length(names) == nrow(ages)){
      ages$names <- names
      ages$col <- col[1:length(names)]
    } else {
      stop("Argument \'names\' should be as long as the number of periods.")
    }
  }
  
  if(direction == "leftwards"){
    stop("NOT IMPLEMENTED YET.")
    labels = as.character(-c(seq_time-present))
  } else{
    if(all(time.seq >= 0)){
      if(time.seq[1] > time.seq[2]){
        labels <- as.character(-time.seq)
      } else {
        labels <- rev(as.character(-time.seq))
      }
    } else {
      labels <- as.character(time.seq)
    }
    
    #min_time <- ifelse(xpd.x, min(ages[,c(1,2)]),0)
    #labels <- seq(-ceiling(present-min_time),0,1)
    #root_time <- c(present - max(abs(labels)))
    #seq_time <- present
    #labels <- as.character(labels)
    #for(i in 1:(length(labels)-1)){
    #  seq_time <- c(seq_time, seq_time[length(seq_time)]-1)
    #}
    #seq_time <- rev(seq_time)
    
  }
  if(is.phylo == F){
    if(any(round(time.seq) != time.seq)){
      labels <- labels[round(time.seq) == time.seq]
      time.seq <- time.seq[round(time.seq) == time.seq]
    }
  } else {
    time.seq2 <- time.seq - time.seq[2]
    labels2 <- time.seq2 - floor(present)
    
    time.seq <- time.seq[round(time.seq2) == time.seq2]
    labels <- as.character(labels2[round(time.seq2) == time.seq2])
  }
  
  if(time.interval == 1){
    axis(side = 1, at = time.seq, labels = rep("", length(time.seq)), cex = cex, col.ticks = "grey", line = -0.5, pos = thickness, cex.axis = cex, padj = padj)
    
    time.seq5 <- time.seq[seq(length(labels),1,-5)]
    labels5 <- labels[seq(length(labels),1,-5)]
    
    axis(side = 1, at = time.seq5, labels = labels5, cex = cex, line = -0.5, pos = thickness, cex.axis = cex, padj = padj)
  } else {
    axis(side = 1, at = time.seq[seq(length(labels),1,-time.interval)], labels = labels[seq(length(labels),1,-time.interval)],
         cex = cex, pos = thickness, cex.axis = cex, padj = padj)
  }
  
  for(i in 1:nrow(ages)){
    
    if(fill){
      col.period <- rep(c("grey92","white"), nrow(ages)/2+1)
      rect(xleft = ages[i, 2], xright = ages[i, 1], ybottom = y2, ytop = plot_dim[4], col = col.period[i], border = col.period[i])
    } else {
      arrows(x0 = ages[i, 2], y0 = y2, x1 = ages[i, 2], plot_dim[4], lty = 2, col = "grey50", length = 0)
    }
    polygon(unlist(rep(ages[i, c(1,2)], each = 2)), c(thickness,y2,y2,thickness), col = ages$col[i], lwd = 0.5)
    
    if(i != nrow(ages)){
      mean_i <- mean(as.numeric(ages[i, c(1,2)]))
      text(mean_i, thickness, labels = ages$names[i], pos = 3,
           cex = cex, adj = c(0.5, 0.5))
    } else {
      if(direction == "leftwards"){
        text(mean(as.numeric(c(ages[i, c(1)],plot_dim[2]))), thickness, labels = ages$names[i], pos = 3,
             cex = cex, adj = c(0.5, 0.5))
      } else {
        # to improve
        if(!is.null(names)){
          x.text <- ifelse(xpd.x, 0, mean(c(plot_dim[1],ages$start[i])))
        }else {
          x.text <- ifelse(xpd.x, mean(c(0,ages$start[i])), mean(c(plot_dim[1],ages$start[i])))
        }
        text(x.text, thickness, labels = ages$names[i], pos = 3, cex = cex, adj = c(0.5, 0.5))
      }
    } 
  }
}
