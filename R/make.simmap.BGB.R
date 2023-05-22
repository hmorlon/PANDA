.get.lineage.vec<-function(tree,node){
	if(is.root(node,tree)){
	return(node)
	}else{
	if(length(node)>1){stop("node should be a single integer value corresponding to node number")}
	E<-tree$edge
	anc.vec = E[match(node,E[,2]),1]
	hold = anc.vec
	while(hold%in%E[,2]){
		hold=E[match(hold,E[,2]),1]
		anc.vec = c(anc.vec,hold)
	}
	return(c(node,unique(anc.vec)))
	}
}


.make.simmap.BGB<-function(anc.phylo,subclade.phylo,ana.events,clado.events,return.mat=FALSE,rnd=rnd){
	
	if(!all(subclade.phylo$tip.label%in%anc.phylo$tip.label)){stop("ERROR: subclade tree contains tips missing from full phylogeny")}
	phylo<-subclade.phylo
	subclade.tips<-phylo$tip.label
	paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
	paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
	totlen<-length(phylo$tip.label)
	root <-totlen  + 1
	heights<-nodeHeights(phylo)
	totheight=max(nodeHeights(phylo))
	totheight.anc=max(nodeHeights(anc.phylo))
	subclade.root=getMRCA(anc.phylo,subclade.tips)
	tips = clado.events$node[which(clado.events$label%in%subclade.tips)]
	if(length(tips)>450){
		nn<-subclade.root
		for(i in 1:(length(getDescendants(anc.phylo,subclade.root))-length(anc.phylo$tip.label))){
			node=subclade.root+i
			node_right=anc.phylo$edge[anc.phylo$edge[,1]==node,2][1]
			node_left=anc.phylo$edge[anc.phylo$edge[,1]==node,2][2]
			dscndnts_left=sum(getDescendants(anc.phylo,node_left)%in%tips)
			dscndnts_right=sum(getDescendants(anc.phylo,node_right)%in%tips)
			if(dscndnts_left >0 && dscndnts_right>0){
				nn<-c(nn,node)
			}
		}
	} else{
		nn<-unique(apply(combn(tips,2),2,function(x)getMRCA(anc.phylo,x)))
	}
	allnodes<-c(tips,nn)
	desc<-vector()
	for(i in 1:length(allnodes)){
		desc<-c(desc,.get.lineage.vec(anc.phylo,allnodes[i]))
	}
	desc=unique(desc)
	desc=desc[c(which(desc>=subclade.root),which(desc<=length(anc.phylo$tip.label)))]
	subnodes=c(desc[which(desc>length(anc.phylo$tip.label))])
	ghost.nodes=desc[which(!desc%in%c(nn,tips))]
	ana.events<-ana.events[which(ana.events$node%in%desc[which(desc!=subclade.root)]),]
	clado.events<-clado.events[which(clado.events$node%in%subnodes),]
	nodes=unique(round(clado.events$node_ht,8))
	events=unique(round(totheight.anc-ana.events$abs_event_time,8))
	nodeDist=sort(c(nodes,events,totheight.anc))
	nodeDiff=diff(nodeDist)
	old.phylo<-phylo
	old.edge<-anc.phylo$edge[which(anc.phylo$edge[,1]%in%nn),]
	other.edge<-phylo$edge
	if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
	old.labels<-as.numeric(names(sort(branching.times(phylo),decreasing=TRUE)))
	if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
		checkmat<-cbind(old.labels,seq(root,length(phylo$tip.label)+phylo$Nnode))
		for(j in 1:phylo$Nnode){phylo$edge[which(other.edge==checkmat[j,1])]<-checkmat[j,2]}
		}	
	mat<-matrix(nrow=0, ncol=3)
	counter_three_letters <- 0
	for(i in 1:phylo$Nnode){
		other<-phylo$edge[phylo$edge[,1]==i+totlen, 2]
		for(b in other){
			int<-matrix(ncol=3)
			int[1]<-i+totlen
			if(b>totlen){
				counter_three_letters <- counter_three_letters + 1
				int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
				int[3]<-b
				} else {
				int[2]<-phylo$tip.label[[b]]
				int[3]<-0 
				}
			mat<-rbind(mat,int)
			}
		}
	desc.mat<-matrix(nrow=length(desc), ncol=2)
	desc.mat[,1]<-desc
	for(i in 1:length(desc)){
		if(desc.mat[i,1]%in%ghost.nodes){
			flag=0
			j=1
			while(flag==0){
			term=anc.phylo$edge[which(anc.phylo$edge[,1]==desc.mat[i,1]),2][j]
			if(term%in%desc){
			flag=1
			}else{
			j=2
			}}
			if(term%in%tips){
				desc.mat[i,2]<-which(anc.phylo$tip.label[term]==phylo$tip.label)
			} 
			if(term%in%nn){
				desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==term),1][1]
			}
			term2=term 
			while(term2%in%ghost.nodes){
				flag=0
				j=1
				while(flag==0){
				term=anc.phylo$edge[which(anc.phylo$edge[,1]==term2),2][j]
				if(term%in%desc){
				flag=1
				}else{
				j=2
				}}
				if(term%in%tips){
					desc.mat[i,2]<-which(anc.phylo$tip.label[term]==phylo$tip.label)
				} 
				if(term%in%nn){
					desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==term),1][1]
				} 
				term2=term
				}
			} 
		if(desc.mat[i,1]%in%tips){
			desc.mat[i,2]<-which(anc.phylo$tip.label[desc.mat[i,1]]==phylo$tip.label)
		}
		if(desc.mat[i,1]%in%nn){
				desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==desc.mat[i,1]),1][1]
		}
	}
	
	
	#while(any(desc.mat[,2]%in%ghost.nodes)){
	#	nxt<-desc.mat[which(desc.mat[,2]%in%ghost.nodes),2]
	#	desc.mat[which(desc.mat[,2]%in%ghost.nodes),2]<-desc.mat[match(nxt ,desc.mat[,1]),2]
	#}
	
	nat<-list()
	nodecount=1
	for(i in 1:length(nodeDiff)){
		if(nodeDist[i]%in%nodes){ #a new species appears, this is a node
			if(clado.events$node[which(round(clado.events$node_ht,8)==round(nodeDist[i],8))]%in%nn){## node has two descendants because of subtree trimming
				IN<-vector()
				node<-totlen+nodecount
				P<-mat[as.numeric(mat[,1])<=(node),c(2,3)]
				IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(node),1])
				#need to find a way to look up old node
				#oldnode=old.edge[which(phylo$edge[,1]==node)][1]
				delbr<-clado.events$node[which(round(clado.events$node_ht,8)==round(nodeDist[i],8))]
				left=clado.events$left_desc_nodes[which(clado.events$node==delbr)] 
				left=desc.mat[which(desc.mat[,1]==left),2]
				right=clado.events$right_desc_nodes[which(clado.events$node==delbr)] 
				right=desc.mat[which(desc.mat[,1]==right),2]
				m<-regexpr("->",clado.events$clado_event_txt[which(clado.events$node==delbr)]) 
				regs<-regmatches(clado.events$clado_event_txt[which(clado.events$node==delbr)],m,invert=TRUE)[[1]][2] 
				p<-regexpr(",",regs)	
				sides<-regmatches(regs,p,invert=TRUE)
				lside<-sides[[1]][1]
				rside<-sides[[1]][2]
				if(i==1){
					geo.vector<-vector(length=length(IN))
					if(left<=totlen){
						geo.vector[which(IN==phylo$tip.label[left])]<-lside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==left])]<-lside
					}
					if(right<=totlen){
						geo.vector[which(IN==phylo$tip.label[right])]<-rside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==right])]<-rside
					}
					nat[[i]]<-cbind(IN,geo.vector)
				} else {
					ogv<-geo.vector			##need to update old geo.vector
					geo.vector<-vector(length=length(IN))
					terms<-which(nat[[i-1]][,1]%in%IN) #this should give the item numbers of old things
					geo.vector[match(nat[[i-1]][,1][terms],IN)]<-ogv[terms]
					if(left<=totlen){
						geo.vector[which(IN==phylo$tip.label[left])]<-lside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==left])]<-lside
					}
					if(right<=totlen){
						geo.vector[which(IN==phylo$tip.label[right])]<-rside
					}else{
						geo.vector[which(IN==mat[,2][mat[,3]==right])]<-rside
					}
					nat[[i]]<-cbind(IN,geo.vector)			
				}} else{
				IN<-nat[[i-1]][,1]
				delbr<-clado.events[which(round(clado.events$node_ht,8)==round(nodeDist[i],8)),1]
				left=clado.events$left_desc_nodes[which(clado.events$node==delbr)] 
				right=clado.events$right_desc_nodes[which(clado.events$node==delbr)] 
				m<-regexpr("->",clado.events$clado_event_txt[which(clado.events$node==delbr)]) 
				regs<-regmatches(clado.events$clado_event_txt[which(clado.events$node==delbr)],m,invert=TRUE)[[1]][2] 
				p<-regexpr(",",regs)	
				sides<-regmatches(regs,p,invert=TRUE)
				lside<-sides[[1]][1]
				rside<-sides[[1]][2]
				if(left%in%desc){ #left branch should be recorded
					hold<-desc.mat[which(desc.mat[,1]==left),2]
					if(hold<=totlen){
						elno<-which(IN==phylo$tip.label[hold])
					}else{
						elno<-which(IN==mat[,2][mat[,3]==hold])
					}
					geo.vector[elno]<-lside
				} else{#right branch should be recorded
					hold<-desc.mat[which(desc.mat[,1]==right),2]
					if(hold<=totlen){
						elno<-which(IN==phylo$tip.label[hold])
					}else{
						elno<-which(IN==mat[,2][mat[,3]==hold])
					}
					geo.vector[elno]<-rside
				}
				nat[[i]]<-cbind(IN,geo.vector)			
		}} else{ # this is a branch change
				IN<-nat[[i-1]][,1]
				#identify which branch changes
				##note should eventually move this up so subtraction doesn't happen every time(for speed up)
				delbr<-which(nodeDist[i]==round(totheight.anc-ana.events$abs_event_time,8))
				if(length(delbr)>1){stop("more than one events at same time")}
				hold<-desc.mat[which(desc.mat[,1]==ana.events[delbr,]$nodenum_at_top_of_branch),2]
				if(hold<=totlen){
					elno<-which(IN==phylo$tip.label[hold])
				}else{
					elno<-which(IN==mat[,2][mat[,3]==hold])
				}
				#identify new state
				geo.vector[elno]<-ana.events[delbr,]$new_rangetxt
				nat[[i]]<-cbind(IN,geo.vector)			
			}
	if(nodeDist[i+1]%in%nodes && clado.events$node[which(round(clado.events$node_ht,8)==round(nodeDist[i+1],8))]%in%nn){nodecount=nodecount+1}
	}			
	
	mat[which(mat[,3]==0),3]<-match(mat[which(mat[,3]==0),2],phylo$tip.label)
	
	maps.list=list()
	
	for(k in 1:length(phylo$edge.length)){
		
		#identify branch from edge matrix
		#lookup the name of this branch in the 'mat' matrix compiled above
		#lookup which nat elements have the name of this branch
		#write a vector of the nodeDiff values named with the ranges for each of these elements

		
		lf<-phylo$edge[k,1]
		ri<-phylo$edge[k,2]
		br<-mat[which(mat[,1]==lf & mat[,3]==ri),2]
		natis<-which(sapply(nat,function(x)br%in%x[,1]))
		out.vec<-nodeDiff[natis]
		name.vec<-vector()
		for(n in 1:length(out.vec)){
			name.vec<-c(name.vec,nat[[natis[n]]][which(nat[[natis[n]]][,1]==br),2])	
		}	
		names(out.vec)<-name.vec
		
		#sum adjacent elements with the same name
		
		out.vec.simple<-vector()
		counter=1
		for(i in 1: length(out.vec)){
			if(i == 1 || i == (counter+1)){
				while((length(out.vec)>counter) && (names(out.vec[i])==names(out.vec[counter+1]))){
					counter=counter+1
					}
				hold<-sum(out.vec[i:counter])
				names(hold)<-names(out.vec[i])
				out.vec.simple<-c(out.vec.simple,hold)	
				}
		}
		
		
		maps.list[[k]]<-out.vec.simple
	
	}

	mapped.edge<-matrix(nrow=dim(phylo$edge)[1],ncol=length(unique(names(unlist(maps.list)))))
	colnames(mapped.edge)<-unique(names(unlist(maps.list)))
	
	for(k in 1:dim(phylo$edge)[1]){
		for(j in 1:dim(mapped.edge)[2]){
			hold<-which(names(maps.list[[k]])==colnames(mapped.edge)[j])
			mapped.edge[k,j]<-ifelse(length(hold)==0,0,sum(maps.list[[k]][hold]))
		}
	}
	
	out<-list(edge=phylo$edge,edge.length=phylo$edge.length,tip.label=phylo$tip.label,Nnode=subclade.phylo$Nnode,maps=maps.list,mapped.edge=mapped.edge,Q="NA",logL="NA")
	class(out)<-c("simmap","phylo")
	
	if(identical(anc.phylo,subclade.phylo)){
		class.object=list(class.object=nat,times=nodeDist,spans=nodeDiff)
	} else {
		class.object=CreateClassObject(out,rnd)
	}
	
	if(!return.mat){
	return(list(geo.simmap=out,class.object=class.object))
	}else{
	return(list(geo.simmap=out,class.object=class.object,mat=mat))
	}}