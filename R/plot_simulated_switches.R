plot_simulated_switches <-
function(n,host_tree,name,index,switches,...){ 
  if(!exists("path")) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  setwd(path)
  pdf(paste(path,"/figures/host_tree_switches_",name,"_",index,".pdf",sep=""))
  tree <- ladderize(host_tree)
  for(i in 1:n){tree$tip.label[i] <- paste("    ",tree$tip.label[i],sep="")}
  print(plot(tree,edge.width=4))
  nodes <- 1:(2*n-2)
  for (i in 1:length(nodes)){nodes[i] <- paste(rep("0",nchar(nodes[length(nodes)])-nchar(nodes[i])),nodes[i],sep="")}
  edgelabels(nodes,frame="circle",col="#78281f",bg="#f39c12")
  add.scale.bar()
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  for (i in 1:ncol(switches)){
    b0=as.numeric(switches[1,i])
    b1=as.numeric(switches[2,i])
    y0=lastPP$yy[lastPP$edge[b0,2]]
    y1=lastPP$yy[lastPP$edge[b1,2]]
    x=as.numeric(switches[3,i])
    print(fancyarrows(x0=x, y0=y0, x1=x, y1=y1, length = 0.25,lwd = 4,type = "harpoon",col="#1e8449"))
  }
  dev.off()
}
