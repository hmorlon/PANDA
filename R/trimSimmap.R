#This is a script that trims a stochastic map to preserve only those branches containing a given category (according to the algorithm in Drury et al. MS)
#The inputs are:
#map (simmap object that is a stochastic map of discrete traits used to identify group of competitors) 
#trim.class (discrete category identifying species in 'data' (e.g., herbivores))

.trimSimmap<-function(map,trim.class){
		trc=trim.class
		smap<-map
		tot.len<-length(smap$tip.label)
		node = tot.len + 1
		
		while(node <= (tot.len+smap$Nnode)){
			descL<-smap$edge[which(smap$edge[,1]==node),2][1]
			Lclade<-unique(c(which(smap$edge[,1]==node)[1],which(smap$edge[,2]%in%getDescendants(smap,descL))))
			
			if(!any(sapply(smap$maps[Lclade],function(x)any(names(x)%in%trc)))){ #left descendants DO NOT contain target trait
		
				if(length(Lclade)==1){ 	#if only descendant is tip, drop entire branch
					smap<-drop.tip.simmap(smap,smap$tip.label[descL])
					tot.len<-length(smap$tip.label)
					node = tot.len + 1
					
				} else{		#if descendant is clade, drop all but one branch (randomly)
					Ltips<-smap$edge[Lclade[which(smap$edge[Lclade,2]<=tot.len)],2]	
					smap<-drop.tip.simmap(smap,smap$tip.label[sample(Ltips,length(Ltips)-1)])
					tot.len<-length(smap$tip.label)
					node = tot.len + 1
				}
		
			} else {
			
			node = node + 1
			
			}
			}	
		
		node = tot.len + 1
		
		
		while(node <= (tot.len+smap$Nnode)){
			descR<-smap$edge[which(smap$edge[,1]==node),2][2]
			Rclade<-unique(c(which(smap$edge[,1]==node)[2],which(smap$edge[,2]%in%getDescendants(smap,descR))))
			
			if(!any(sapply(smap$maps[Rclade],function(x)any(names(x)%in%trc)))){ #right descendants DO NOT contain target trait
		
				if(length(Rclade)==1){ 	#if only descendant is tip, drop entire branch
					smap<-drop.tip.simmap(smap,smap$tip.label[descR])
					tot.len<-length(smap$tip.label)
					node = tot.len + 1
					
				} else{		#if descendant is clade, drop all but one branch (randomly)
					Rtips<-smap$edge[Rclade[which(smap$edge[Rclade,2]<=tot.len)],2]	
					smap<-drop.tip.simmap(smap,smap$tip.label[sample(Rtips,length(Rtips)-1)])
					tot.len<-length(smap$tip.label)
					node = tot.len + 1
				}
		
			} else {
			
			node = node + 1
			
			}
			}	
	return(smap)
	}	

