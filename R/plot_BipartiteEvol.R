

######## Plot the result of the model ############################

# gen an oject returned by make.gen
# spec the corresponding object returned by define.species
#


plot_div.BipartiteEvol=function(gen,spec,trait.id,lwdgen=1,lwdsp=lwdgen,scale=NULL){
  
  oldpar <- par(no.readonly = TRUE) # current par
  on.exit(par(oldpar)) # restore the previous par at the end of the function
  
  # cut the plot window
  layout(mat = matrix(c(1,1,8,2,2,3,3,8,4,4,5,7,7,7,6),ncol=5,byrow = TRUE),
         widths = c(3,2,1.5,2,3), heights = c(3,3,2))
  
  # define the scale range
  if(is.null(scale)) {
    scale=range(c(gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }else{
    scale=range(c(scale,gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }
  if(scale[1]==scale[2]){scale[2]=scale[1]+1}
  
  # plot the genealogies
  #for P
  plot_with_trait(gen$P,gen$P$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("P",gen$P$ini))
  # for H
  plot_with_trait(gen$H,gen$H$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("H",gen$H$ini))
  
  # plot the species trees
  # for P
  if(!is.null(spec$Pphylo$tree)){
    plot_with_trait(spec$Pphylo$tree,spec$Pphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Pphylo$tree$tip.label))
  }else{plot(1,spec$Pphylo$mean.trait[trait.id])}
  # for H
  if(!is.null(spec$Hphylo$tree)){
    plot_with_trait(spec$Hphylo$tree,spec$Hphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Hphylo$tree$tip.label))
  }else{plot(1,spec$Hphylo$mean.trait[trait.id])}
  
  # plot the trait densities
  densplot(mcmc(gen$P$x.tip[trait.id,]))
  densplot(mcmc(gen$H$x.tip[trait.id,]))
  
  # plot the traits with their correlation 
  plot(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]}),xlab="P trait",ylab="H trait")
  lines(c(scale[1],scale[2]),c(scale[1],scale[2]),type = 'l',col="red")
  cor=try(cor(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]})))
  if (inherits(cor,"try_error")){
    title(main="_") 
  }else{
    title(main=cor)  
  }
  
  # add legend
  ma = par("mar")
  par(mar = 0.5+c(0,0,0,0))
  plot.new()
  rasterImage(as.raster(colorRampPalette(rev(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4")))( 100 )),0,0,0.5,1)
  text(x=0.7, y = seq(0,1,l=5), labels = round(seq(scale[1],scale[2],l=5),digits = 2))
  axis(side=4,at=seq(0,1,l=5),pos=0.5,labels = FALSE)
  
  par(mar = ma)
  
  
}



######## Plot a phylogeny with colored branches ##############################

plot_with_trait=function(phylo,rate,scale=NULL,lwd=1,direction="rightwards"){
  
  # define scale
  scale=range(c(scale,rate))
  if(scale[1]==scale[2]){scale[2]=scale[1]+1}
  
  # define color palette
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
  
  # define branch colors
  col = round( (rate - min(scale)) / diff(range(scale))*99   )+1
  
  plot(phylo, edge.color = Colors[col], show.tip.label = FALSE,edge.width =lwd,direction=direction)
}

######## Plot the colored spatial matrix ################################

spatial.plot=function(out,trait.id,scale=NULL,nx=NULL, sort_trait = FALSE){
  
  oldpar <- par(no.readonly = TRUE) # current par
  on.exit(par(oldpar)) # restore the previous par at the end of the function

  if(is.null(nx))   nx=sqrt(length(out$P[1,]))
  scale=range(c(scale,out$P[trait.id,],out$H[trait.id,]))
  if(scale[1]==scale[2]){
    scale[1]=scale[1]-0.5
    scale[2]=scale[2]+0.5
  }
  MP=max(out$P[trait.id,])
  mP=min(out$P[trait.id,])
  MH=max(out$H[trait.id,])
  mH=min(out$H[trait.id,])
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
  par(mfrow=c(1,2))
  if (sort_trait){
    order_P = order(out$P[trait.id,])
    image(matrix(out$P[trait.id,][order_P],ncol=nx,byrow = TRUE),col=Colors[(round( (mP - min(scale)) / diff(range(scale))*99   )+1):(round( (MP - min(scale)) / diff(range(scale))*99   )+1)],
          main="",axes=FALSE)
    image(matrix(out$H[trait.id,][order_P],ncol=nx,byrow = TRUE),col=Colors[(round( (mH - min(scale)) / diff(range(scale))*99   )+1):(round( (MH - min(scale)) / diff(range(scale))*99   )+1)]
          ,main="",axes=FALSE)
  }else{
    image(matrix(out$P[trait.id,],ncol=nx,byrow = TRUE),col=Colors[(round( (mP - min(scale)) / diff(range(scale))*99   )+1):(round( (MP - min(scale)) / diff(range(scale))*99   )+1)],
          main="",axes=FALSE)
    image(matrix(out$H[trait.id,],ncol=nx,byrow = TRUE),col=Colors[(round( (mH - min(scale)) / diff(range(scale))*99   )+1):(round( (MH - min(scale)) / diff(range(scale))*99   )+1)]
          ,main="",axes=FALSE)
  }
}

######## Plot the network ##################################

plot_network=function(link,spec,trait.id=1,method="bipartite",order=TRUE,scale=c()){
  
  if(method=="bipartite"){
    if(is.vector(spec$Pphylo$mean.trait)){
      spec$Pphylo$mean.trait = matrix(spec$Pphylo$mean.trait, nrow = 1)
      spec$Hphylo$mean.trait = matrix(spec$Hphylo$mean.trait, nrow = 1)
      
    }
    
    Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 )
    scale=range(c(scale,spec$Pphylo$mean.trait[trait.id,],spec$Hphylo$mean.trait[trait.id,]))
    col.P=rep(0,length(spec$Pphylo$tree$tip.label))
    col.H=rep(0,length(spec$Hphylo$tree$tip.label))
    col.P[spec$Pphylo$tree$tip.label]=Colors[(round( (spec$Pphylo$mean.trait[trait.id,] - min(scale)) / diff(range(scale))*99   )+1)]
    col.H[spec$Hphylo$tree$tip.label]=Colors[(round( (spec$Hphylo$mean.trait[trait.id,] - min(scale)) / diff(range(scale))*99   )+1)]
    try(plotweb(as.matrix(link),lower_color=col.P,higher_color = col.H,link_color = "lightgray",empty=TRUE))
  }else{
    Mat=as.matrix(link)
    if(order) Mat=Mat[order(rowSums(Mat),decreasing = TRUE),order(colSums(Mat),decreasing = FALSE)]
    image(log(Mat+1), axes = FALSE, col = grey(seq(1, 0.2, length = 256)))
    title(xlab = "P", ylab="H",outer = FALSE,line=1)
  }
}



######## Plot the results with the network ###############

plot_net.BipartiteEvol=function(gen,spec,trait.id, link,out,lwdgen=1,lwdsp=lwdgen,scale=NULL,nx=NULL,cor=FALSE,network.method="bipartite",spatial=FALSE){
  
  oldpar <- par(no.readonly = TRUE) # current par
  on.exit(par(oldpar)) # restore the previous par at the end of the function
  
  # cut the plot window
  if(spatial){
    layout(mat = matrix(c(1,1,8,2,2,2,2,3,3,8,4,4,4,4,5,7,7,7,6,6,6,9,9,9,10,10,11,11),ncol=7,byrow = TRUE),
           widths = c(3,2,1.5,2,0.5,0.5,2), heights = c(3,3,2.5,4))
  }else{
    layout(mat = matrix(c(1,1,8,2,2,2,2,3,3,8,4,4,4,4,5,7,7,7,6,6,6,9,9,9,9,9,9,9),ncol=7,byrow = TRUE),
           widths = c(3,2,1.5,2,0.5,0.5,2), heights = c(3,3,2.5,4))
  }
  
  # define the scale range
  if(is.null(scale)) {
    scale=range(c(gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }else{
    scale=range(c(scale,gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }
  if(scale[1]==scale[2]){scale[2]=scale[1]+1}
  
  # plot the genealogies
  #for P
  plot_with_trait(gen$P,gen$P$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("P",gen$P$ini))
  # for H
  plot_with_trait(gen$H,gen$H$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("H",gen$H$ini))
  
  # plot the species trees
  # for P
  if(!is.null(spec$Pphylo$tree)){
    plot_with_trait(spec$Pphylo$tree,spec$Pphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Pphylo$tree$tip.label))
  }else{plot(1,spec$Pphylo$mean.trait[trait.id])}
  # for H
  if(!is.null(spec$Hphylo$tree)){
    plot_with_trait(spec$Hphylo$tree,spec$Hphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Hphylo$tree$tip.label))
  }else{plot(1,spec$Hphylo$mean.trait[trait.id])}
  
  # plot the trait densities
  densplot(mcmc(gen$P$x.tip[trait.id,]))
  densplot(mcmc(gen$H$x.tip[trait.id,]))
  
  # plot the traits with their correlation 
  if(cor){
    plot(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]}),xlab="P trait",ylab="H trait")
    lines(c(scale[1],scale[2]),c(scale[1],scale[2]),type = 'l',col="red")
  }else{
    if(is.null(spec)){
      plot(1,1)
    }else{
      plot_species(spec,cex=2,xlab=c(),ylab=c(), net = link)}
  }
  cor=try(cor(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]})))
  if (inherits(cor,"try_error")){
    title(main="_") 
  }else{
    title(main=cor)  
  }
  
  # add legend
  ma = par("mar")
  par(mar = 0.5+c(0,0,0,0))
  plot.new()
  rasterImage(as.raster(colorRampPalette(rev(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4")))( 100 )),0,0,0.5,1)
  text(x=0.7, y = seq(0,1,l=5), labels = round(seq(scale[1],scale[2],l=5),digits = 2))
  axis(side=4,at=seq(0,1,l=5),pos=0.5,labels = FALSE)
  
  par(mar = ma)
  
  # plot.new()
  t=try(plot_network(link,spec,trait.id,method=network.method))
  if(inherits(t,"try-error")) plot.new()
  
  if(spatial){
    if(is.null(nx))   nx=sqrt(length(out$P[1,]))
    scale=range(c(scale,out$P[trait.id,],out$H[trait.id,]))
    MP=max(out$P[trait.id,])
    mP=min(out$P[trait.id,])
    MH=max(out$H[trait.id,])
    mH=min(out$H[trait.id,])
    Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
      image(matrix(out$P[trait.id,],ncol=nx,byrow = TRUE),col=Colors[(round( (mP - min(scale)) / diff(range(scale))*99   )+1):(round( (MP - min(scale)) / diff(range(scale))*99   )+1)],
            main="",axes=FALSE)
      image(matrix(out$H[trait.id,],ncol=nx,byrow = TRUE),col=Colors[(round( (mH - min(scale)) / diff(range(scale))*99   )+1):(round( (MH - min(scale)) / diff(range(scale))*99   )+1)]
            ,main="",axes=FALSE)
  }
}



######## Plot species in the phenotype space ##############



plot_species=function(spec,trait.id=1:3,net = NULL,...){
  
  D = max(1,nrow(spec$Pphylo$mean.trait))
  trait.id[trait.id>D] = D
  gr = "gray50"
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 )
  scale=range(c(spec$Pphylo$mean.trait[trait.id[3],],spec$Hphylo$mean.trait[trait.id[3],]))
  col.P=Colors[(round( (spec$Pphylo$mean.trait[trait.id[3],] - min(scale)) / diff(range(scale))*99   )+1)]
  col.H=Colors[(round( (spec$Hphylo$mean.trait[trait.id[3],] - min(scale)) / diff(range(scale))*99   )+1)]
  
  scale_x=range(c(spec$Pphylo$mean.trait[trait.id[1],],spec$Hphylo$mean.trait[trait.id[1],]))
  scale_y=range(c(spec$Pphylo$mean.trait[trait.id[2],],spec$Hphylo$mean.trait[trait.id[2],]))
  plot(c(),c(),xlim=scale_x,ylim=scale_y)
  if(!is.null(net)){
    for (i in 1:nrow(net)){
      for (j in 1:ncol(net)){
        if (net[i,j]>0){
          lines(c(spec$Pphylo$mean.trait[trait.id[1],i],spec$Hphylo$mean.trait[trait.id[1],j]),
                c(spec$Pphylo$mean.trait[trait.id[2],i],spec$Hphylo$mean.trait[trait.id[2],j]),
                col = gr)
        }
      }
    }
  }
  points(spec$Pphylo$mean.trait[trait.id[1],],spec$Pphylo$mean.trait[trait.id[2],],pch=16,col=col.P,...)
  points(spec$Hphylo$mean.trait[trait.id[1],],spec$Hphylo$mean.trait[trait.id[2],],pch=18,col=col.H,...)
}
