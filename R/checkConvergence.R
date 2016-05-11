norm2=function(x){sqrt(sum(x^2))}

check.convergence=function(DEopt,ylim1=NULL,ylim2=NULL,i0=1,phylo,rel=T){
  ancestors=get_ancestors(phylo)
  best=DEopt$optim$bestmem
  if(rel){dif=sapply(2:nrow(DEopt$member$bestmemit),function(i){norm2(get_rates(lambda = DEopt$member$bestmemit[i,],Ancestors = ancestors,phylo=phylo)[-1]-
                                                                get_rates(lambda = DEopt$member$bestmemit[i-1,],Ancestors = ancestors,phylo=phylo)[-1])})
    par(mfrow=c(1,2))
    plot((i0+1):length(dif),dif[(i0+1):length(dif)],ylim=ylim1,xlab="iteration",ylab="distance between iterations")
  }else{
    dif=sapply(1:nrow(DEopt$member$bestmemit),function(i){norm2(DEopt$member$bestmemit[i,]-best)})
    par(mfrow=c(1,2))
    plot(i0:length(dif),dif[i0:length(dif)],ylim=ylim1,xlab="iteration",ylab="distance to best")
    }
  # 
 
  plot(i0:length(dif),DEopt$member$bestvalit[i0:length(dif)],ylim=ylim2,xlab="iteration",ylab="-logLik")
  par(mfrow=c(1,1))
}

distance.to.true=function(DEopt,true.rate,i0=1,ylim=NULL,rel=F,phylo=NULL,mcmc=F,not.lambda=c(1)){
  if(mcmc){
    if(rel){
      not.l=c(not.lambda,ncol(DEopt)-1,ncol(DEopt))
      nnode=phylo$Nnode
      Ancestors=get_ancestors(phylo)
      dif=sapply(1:nrow(DEopt),function(i){norm2(get_rates(phylo = phylo,lambda = DEopt[i,-not.l],Ancestors=Ancestors)[-1]-true.rate)})
    }else{
      not.l=c(not.lambda,ncol(DEopt)-1,ncol(DEopt))
      dif=sapply(1:nrow(DEopt),function(i){norm2(DEopt[i,-not.l]-true.rate)})}
  }else{
    if(rel){
      Ancestors=get_ancestors(phylo)
      dif=sapply(1:nrow(DEopt$member$bestmemit),function(i){norm2(get_rates(phylo = phylo,lambda = DEopt$member$bestmemit[i,],Ancestors=Ancestors)[-1]-true.rate)})
    }else{
      dif=sapply(1:nrow(DEopt$member$bestmemit),function(i){norm2(DEopt$member$bestmemit[i,]-true.rate)})}}
  plot(i0:length(dif),dif[i0:length(dif)],ylim=ylim,xlab="iteration",ylab="distance to true.rate")
}

plot.with.rate=function(phylo,rate1,rate2=NULL,same.scale=T,main=NULL,lwd=1){
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
  if(is.null(rate2)){
    if(isTRUE(all.equal(rep(rate1[1],length(rate1)),rate1))){
      col=rep(1,length(rate1))
      plot(phylo, edge.color = Colors[col], show.tip.label = F,main=main,edge.width =lwd)
      image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
    }else{
      col = round( (rate1 - min(rate1)) / diff(range(rate1))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = F,main=main,edge.width =lwd)
      image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)}
  }else{
    if(same.scale){
      min=min(min(rate1),min(rate2))
      max=max(max(rate1),max(rate2))
      par(mfrow=c(1,2))
      col = round(( (rate1 - min) / (max-min))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
      col = round(( (rate2 - min) / (max-min))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
      par(mfrow=c(1,1))
      image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T)
    }else{
      par(mfrow=c(1,2))
      if(isTRUE(all.equal(rep(rate1[1],length(rate1)),rate1))){
        col=rep(1,length(rate1))
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
      }else{
        col = round(( (rate1 - min(rate1)) / (max(rate1)-min(rate1)))*99   )+1
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        image.plot(z = c(min(rate1),max(rate1)),col = Colors, horizontal=T,legend.only = T)
      }
      if(isTRUE(all.equal(rep(rate2[1],length(rate2)),rate2))){
        col=rep(1,length(rate2))
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        image.plot(z = c(rate2[1],2*rate2[1]),col = Colors, horizontal=T,legend.only = T)
      }else{
        col = round(( (rate2 - min(rate2)) / (max(rate2)-min(rate2)))*99   )+1
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        image.plot(z = c(min(rate2),max(rate2)),col = Colors, horizontal=T,legend.only = T)
      }
    }
    par(mfrow=c(1,1))
    }
}

add.iterations=function(DEopt,target,control=list(maxit=100)){
  initialpop=DEopt$member$pop
  Control=control
  Control$initialpop=initialpop
  res=DEoptim(target,upper=DEopt$member$upper,lower=DEopt$member$lower,control=Control)
  res$member$bestmemit=rbind(DEopt$member$bestmemit,res$member$bestmemit)
  res$member$bestvalit=c(DEopt$member$bestvalit,res$member$bestvalit)
  return(res)
}
