resume_switches <-
function(name,index,est_ksi,nb_tree,host_tree){
  if(est_ksi>0){
    
    transparent_theme_no <- theme(panel.grid = element_blank(),
                                  axis.line = element_line("black"),
                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                  plot.background = element_rect(fill = "transparent",colour = NA),
                                  legend.key=element_blank(),legend.background=element_blank())
    
    #### RESUME ESTIMATED SWITCHES
    
    trees_ll <- read.table(paste("results/optim_ll_",name,"_",index,"_",est_ksi,".txt",sep=""))
    trees_ll <- cbind(1:length(trees_ll$V1),trees_ll$V1)
    colnames(trees_ll) <- c("index","minloglik")
    load(file=paste("simulated_trees/simulated_switches_",name,"_",est_ksi,".RData",sep=""))
    
    max_likelihood <- max(-trees_ll$minloglik)
    proba_trees <- exp(-trees_ll$minloglik-max_likelihood)*exp(max_likelihood)
    if (all(proba_trees==0)){proba_trees <- exp(-trees_ll$minloglik-max_likelihood)}
    
    list_tree <- sample(x = trees_ll$index, prob=proba_trees, size = nb_tree, replace = T)
    
    list_switches_selected <- list_switches
    for (i in 1:length(list_tree)){
      index_tree <- list_tree[i]
      list_switches_selected[[i]] <- rbind(rbind(list_switches[[index_tree]],rep(trees_ll$minloglik[index_tree],ncol(list_switches[[index_tree]]))),rep(index_tree,ncol(list_switches[[index_tree]])))}
    selected_switches <- matrix(unlist(list_switches_selected),nrow=5)    
    
    table_inferred_switches <- data.frame(t(selected_switches[c(1,2,3,4,5),]))
    colnames(table_inferred_switches) <- c("branch_depature","branch_arrival","position","likelihood","tree")
    
    load(file=paste("simulated_trees/simulated_switches_",name,"_",est_ksi,".RData",sep=""))
    table_simul_switches <- data.frame(t(matrix(unlist(list_switches),nrow=3)))
    
    colnames(table_simul_switches) <- c("branch_depature","branch_arrival","position")
    
    #######  PLOT MAIN SWITCHES
    tree <- host_tree
    plot(host_tree)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    
    table_inferred_switches$name <- paste(table_inferred_switches[,1],table_inferred_switches[,2],sep="_")
    
    combinations <-  strsplit(unique(table_inferred_switches$name),split="_")
    deviation_uniform <- c()
    for (i in 1:length(combinations)){
      deviation_uniform[i] <- length(which(table_inferred_switches$name==paste(combinations[[i]][1],combinations[[i]][2],sep="_")))
    }
    combinations_order <- combinations[order(deviation_uniform)]
    
    t_col <- function(color, percent = 50, name = NULL) {rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],alpha = (100-percent)*255/100,names = name)
    return(t.col)}
    
    pdf(paste("infer_PHS/host_tree_inferred_switches_uniform_",name,"_",index,".pdf",sep=""),width = (10+Ntip(host_tree))/2 ,height = (10+Ntip(host_tree))/2 )
    plot(tree,edge.width=4)
    nodes <- 1:(2*Ntip(tree)-2)
    for (i in 1:length(nodes)){nodes[i] <- paste(paste(rep("0",nchar(nodes[length(nodes)])-nchar(nodes[i])),collapse = ""),nodes[i],sep="")}
    edgelabels(nodes,frame="circle",col="#1b2631",bg="#dc7633")
    add.scale.bar()
    
    Colors = colorRampPalette(c(t_col("#f9e79f",50),t_col("salmon1",90),t_col("darkorange",30), t_col("red",10),"red4"))(100)
    col = sort(round( (deviation_uniform - min(deviation_uniform)) / diff(range(deviation_uniform))*99   )+1)
    Width_Arrow <- seq(0.5,5,0.1)
    width_arrow = sort(round( (deviation_uniform - min(deviation_uniform)) / diff(range(deviation_uniform))*(length(Width_Arrow)-1)   )+1)
    for (i in 1:length(combinations_order)){
      b0=as.numeric(combinations_order[[i]][1])
      b1=as.numeric(combinations_order[[i]][2])
      y0=lastPP$yy[lastPP$edge[b0,2]]
      y1=lastPP$yy[lastPP$edge[b1,2]]
      min_x=max(lastPP$xx[lastPP$edge[b1,1]],lastPP$xx[lastPP$edge[b0,1]])
      max_x=min(lastPP$xx[lastPP$edge[b1,2]],lastPP$xx[lastPP$edge[b0,2]])
      x=min_x+(max_x-min_x)/2
      x=rnorm(1,x,abs((max_x-min_x)/10))
      fancyarrows(x0=x, y0=y0, x1=x, y1=y1, length = 0.25,lwd = Width_Arrow[width_arrow[i]],type = "harpoon",col=Colors[col[i]])
      fancyarrows(x0=x, y0=y1, x1=x, y1=y0, length = 0.25,lwd = Width_Arrow[width_arrow[i]],type = "harpoon",col=Colors[col[i]])
    }
    dev.off()
    
    #######  PLOT POSITION
    
    pdf(paste("infer_PHS/positions_inferred_switches_uniform_",name,"_",index,".pdf",sep=""),width=7,height=5)
    print(ggplot() +
            geom_density(data=table_inferred_switches,aes(x=position), fill = "#d35400",color = "#d35400",alpha=0.5) +
            geom_density(data=table_simul_switches,aes(x=position), fill = "#1b4f72",color = "#1b4f72",alpha=0.5) +
            annotate("text", x = max(table_simul_switches$position)/100, y = Inf, label = "Estimated",color="#d35400",fontface="bold",hjust = 0,size=4, vjust=2)+
            annotate("text", x = max(table_simul_switches$position)/100, y = Inf, label = "Uniform",color="#1b4f72",fontface="bold",hjust = 0,size=4, vjust=4)+
            transparent_theme_no+labs(x = "Switches positions", y="Density"))
    dev.off()
  }
}
