.ReconcileGeoObjectSMatrix<-function(geo.object,S.matrix,phylo,rnd=5){

	geot<-round(geo.object$times,rnd)
	clat<-round(S.matrix$times,rnd)
	nodeDist<-sort(unique(c(geot,clat)))
	nodeDiff<-diff(nodeDist)
	if(any(nodeDiff<= (2*(10^-rnd)))){
		nodes<-round(max(branching.times(phylo))-branching.times(phylo),rnd)
		gch<-which(!nodes%in%geot)
		chg<-which(round(geot,rnd-2)%in%round(nodes[gch],rnd-2))
		if(length(chg)!=length(gch)){
			stop("potential rounding error, two time bins very similar, try changing rnd digits")
			}else{
			geot[chg]<-nodes[gch]
			}
	
		nodeDist<-sort(unique(c(geot,clat)))
		nodeDiff<-diff(nodeDist)
		if(any(nodeDiff<= (2*(10^-rnd)))){
		
		sch<-which(!nodes%in%clat)
		chs<-which(round(clat,rnd-2)%in%round(nodes[sch],rnd-2))
		if(length(chs)!=length(sch)){
			stop("potential rounding error, two time bins very similar, try changing rnd digits")
			}else{
			clat[chs]<-nodes[sch]
			}
		
		nodeDist<-sort(unique(c(geot,clat)))
		nodeDiff<-diff(nodeDist)
		if(any(nodeDiff<= (2*(10^-rnd)))){
			stop("potential rounding error, two time bins very similar, try changing rnd digits")
		}
		}
	
	}
	
geo.new<-list()
sma.new<-list()
#initialize counter for class object and geo object
u<-0 #class object
y<-0 #geo object

for(i in 1:length(nodeDiff)){

	if((nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #if timing is the same for both
		u = u+1
		y = y+1
		geo.new[[i]]<-geo.object$geography.object[[y]]
		sma.new[[i]]<-S.matrix$S.matrix[[u]]
	}
	if((nodeDist[i]%in%geot) && (!nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
		y = y+1
		geo.new[[i]]<-geo.object$geography.object[[y]]
		sma.new[[i]]<-S.matrix$S.matrix[[u]]
	}
	if((!nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
		u = u+1
		geo.new[[i]]<-geo.object$geography.object[[y]]
		sma.new[[i]]<-S.matrix$S.matrix[[u]]
	}
	if(!all(colnames(sma.new[[i]])==colnames(geo.new[[i]]))){stop("names don't correspond between matrices")}
	}

return(list(geo.object=list(geography.object=geo.new,times=nodeDist,spans=nodeDiff),S.matrix=list(S.matrix=sma.new,times=nodeDist,spans=nodeDiff)))


}