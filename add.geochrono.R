add.geochrono <- function(Y1, quaternary = F, plot_dim = par("usr"), plot.obj.phylo = NULL, present = NULL,
                          xpd.x = T, names = NULL,
                          cex = 1, root.age = NULL, direction = "rightwards"){
   # to add as a argument
  if(is.null(plot.obj.phylo)){
    plot.obj.phylo<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  } 
  
  if(is.null(present)){ # to counter plot.obj.phylo
    if(is.null(root.age)){
      if("xx" %in% names(plot.obj.phylo)){
        present <- plot.obj.phylo$present[1]
      } else {
        present <- 40 # to adjust for other plots
      }
    } else {
      present <- root.age
    }
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
  
  y1 <- plot_dim[3] - abs(Y1)
  #y1 <- Y1
  
  y2 <- plot_dim[3]

  ages$start  <- ages$start + present
  ages$end  <- ages$end + present  

  if(direction == "leftwards"){
    xlimit <- ceiling(c(max(plot.obj.phylo$xx)+present)/5)*5
    ages <- ages[which(ages$start < xlimit),]
  } else {
    ages <- ages[c(1:which(ages$end<0)[1]),]
  }
  
  col.period <- rep(c("grey92","white"), nrow(ages)/2+1)
  
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
    seq_time
    labels = as.character(-c(seq_time-present))
  } else{
    min_time <- ifelse(xpd.x, min(ages[,c(1,2)]),0)
    labels <- seq(-ceiling(present-min_time),0,1)
    #root_time <- c(present - max(abs(labels)))
    seq_time <- present
    labels <- as.character(labels)
    for(i in 1:(length(labels)-1)){
      seq_time <- c(seq_time, seq_time[length(seq_time)]-1)
    }
    seq_time <- rev(seq_time)
    
  }
  
  axis(side = 1, at = seq_time, labels = rep("", length(seq_time)), cex = cex, col.ticks = "grey", line = -0.5, pos = y1, cex.axis = cex)
  axis(side = 1, at = seq_time[seq(length(labels),1,-5)], labels = labels[seq(length(labels),1,-5)], cex = cex, line = -0.5, pos = y1, cex.axis = cex)
  for(i in 1:nrow(ages)){
    
    rect(xleft = ages[i, 2], xright = ages[i, 1], ybottom = y2, ytop = plot_dim[4], col = col.period[i], border = col.period[i])
    polygon(unlist(rep(ages[i, c(1,2)], each = 2)), c(y1,y2,y2,y1), col = ages$col[i], lwd = 0.5)
    
    if(i != nrow(ages)){
      mean_i <- mean(as.numeric(ages[i, c(1,2)]))
      text(mean_i, y1, labels = ages$names[i], pos = 3,
           cex = cex)
    } else {
      if(direction == "leftwards"){
        text(mean(as.numeric(c(ages[i, c(1)],plot_dim[2]))), y1, labels = ages$names[i], pos = 3,
             cex = cex)
      } else {
        # to improve
        if(!is.null(names)){
          x.text <- ifelse(xpd.x, 0, mean(c(plot_dim[1],ages$start[i])))
        }else {
          x.text <- ifelse(xpd.x, mean(c(0,ages$start[i])), mean(c(plot_dim[1],ages$start[i])))
        }
        text(x.text, y1, labels = ages$names[i], pos = 3, cex = cex)
      }
    } 
  }

  #mtext(text = "Time (Myrs)", side = 1, line = 3.5, at = present/2, cex = cex*2/3)
}