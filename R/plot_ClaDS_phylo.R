
plot_ClaDS_phylo=function(phylo, rates, rates2=NULL, same.scale=TRUE, main=NULL, lwd=2, log=TRUE, show.tip.label=FALSE, ...){
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
  if(is.null(rates2)){
    if(log) rates=log(rates)
    if(isTRUE(all.equal(rep(as.numeric(rates[1]),length(rates)),as.numeric(rates)))){
      col=rep(1,length(rates))
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
      if(log){
        image.plot(z = c(exp(rates[1]),2*exp(rates[1])),col = Colors, horizontal=T,legend.only = T)
      }else{
        image.plot(z = c(rates[1],2*rates[1]),col = Colors, horizontal=T,legend.only = T)
      }
    }else{
      col = round( (rates - min(rates)) / diff(range(rates))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
      if(log){
        min=min(rates)
        max=max(rates)
        m10=floor(min/log(10))
        M10=ceiling(max/log(10))
        if((M10-m10)<4){
          ticks=c(1,2,5)
        }else{
          ticks=1
        }
        ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
        lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
        if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
        image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
      }else{
        image.plot(z = as.matrix(rates),col = Colors, horizontal=T,legend.only = T)
      }
    }
  }else{
    if(log){
      rates=log(rates)
      rates2=log(rates2)
    }
    if(same.scale){
      min=min(min(rates),min(rates2))
      max=max(max(rates),max(rates2))
      par(mfrow=c(1,2))
      col = round(( (rates - min) / (max-min))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
      col = round(( (rates2 - min) / (max-min))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
      par(mfrow=c(1,1))
      if(log){
        m10=floor(min/log(10))
        M10=ceiling(max/log(10))
        if((M10-m10)<4){
          ticks=c(1,2,5)
        }else{
          ticks=1
        }
        ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
        lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
        if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
        # ticks=seq(min,max,length.out = 5)
        image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
      }else{
        image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T)
      }
    }else{
      par(mfrow=c(1,2))
      if(isTRUE(all.equal(rep(rates[1],length(rates)),rates))){
        col=rep(1,length(rates))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
        if(log){
          
          image.plot(z = c(exp(rates[1]),2*exp(rates[1])),col = Colors, horizontal=T,legend.only = T)
        }else{
          image.plot(z = c(rates[1],2*rates[1]),col = Colors, horizontal=T,legend.only = T)
        }
      }else{
        col = round(( (rates - min(rates)) / (max(rates)-min(rates)))*99   )+1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
        if(log){
          min=min(rates)
          max=max(rates)
          m10=floor(min/log(10))
          M10=ceiling(max/log(10))
          if((M10-m10)<4){
            ticks=c(1,2,5)
          }else{
            ticks=1
          }
          ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
          lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
          if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
          image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
        }else{
          image.plot(z = as.matrix(rates),col = Colors, horizontal=T,legend.only = T)
        }
      }
      if(isTRUE(all.equal(rep(rates2[1],length(rates2)),rates2))){
        col=rep(1,length(rates2))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
        if(log){
          image.plot(z = c(exp(rates2[1]),2*exp(rates2[1])),col = Colors, horizontal=T,legend.only = T)
        }else{
          image.plot(z = c(rates2[1],2*rates2[1]),col = Colors, horizontal=T,legend.only = T)
        }
      }else{
        col = round(( (rates2 - min(rates2)) / (max(rates2)-min(rates2)))*99   )+1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, main=main, edge.width = lwd, ...)
        if(log){
          min=min(rates2)
          max=max(rates2)
          m10=floor(min/log(10))
          M10=ceiling(max/log(10))
          if((M10-m10)<4){
            ticks=c(1,2,5)
          }else{
            ticks=1
          }
          ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
          lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
          if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
          image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
        }else{
          image.plot(z = as.matrix(rates2),col = Colors, horizontal=T,legend.only = T)
        }
      }
    }
    par(mfrow=c(1,1))
  }
}
