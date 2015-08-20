CreateGeoObject<-function(phylo,maps.object=NULL){
	if(any(grepl("___",phylo$tip.label))){stop("script will not work with '___' in tip labels; remove extra underscores")}
	paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
	paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
	nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
	root <- length(phylo$tip.label) + 1
	heights<-nodeHeights(phylo)
	for (i in 1:dim(phylo$edge)[1]){
		nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
	}
	nodeDist<-c(nodeDist,max(heights))
	nodeDiff<-diff(nodeDist)
	flag=0
	if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
		node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = phylo$Nnode))
		node.order<-node.order+length(phylo$tip.label)
		old.edge<-phylo$edge
		old.phylo<-phylo
		phylo$edge[,1]<-node.order
		for(j in 1:length(phylo$edge[,2])){
			if(phylo$edge[j,2]>length(phylo$tip.label)){
				#match number order in old edge
				#lookup value in new edge
				#replace with value
				phylo$edge[j,2]<-phylo$edge[,1][match(phylo$edge[j,2],old.edge[,1])]
				}
			}
		nodeDist<-vector()
		for (i in 1:dim(phylo$edge)[1]){
			nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
			}
		nodeDist<-c(nodeDist,max(heights))
		nodeDiff<-diff(nodeDist)
		flag=1
	}
	mat<-matrix(nrow=0, ncol=3)
	counter_three_letters <- 0
	for(i in 1:phylo$Nnode){
		other<-phylo$edge[phylo$edge[,1]==i+length(phylo$tip.label), 2]
		for(b in other){
			int<-matrix(ncol=3)
			int[1]<-i+length(phylo$tip.label)
			if(b>length(phylo$tip.label)){
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
	nat<-list()
	for(i in 1:length(nodeDiff)){
		if(i==1){
		nat[[i]]<-list(mat[mat[,1]==(length(phylo$tip.label)+i),2])} else {
		IN<-vector()
		P<-mat[as.numeric(mat[,1])<=(length(phylo$tip.label)+i),c(2,3)]
		IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(length(phylo$tip.label)+i),1])
		nat[[i]]<-list(IN)
		}
		}	
	for(i in 2:length(nodeDiff)){			##THIS LOOP checks for an error
		if(length(unlist(nat[[i]]))!=(length(unlist(nat[[i-1]]))+1)){
			print(paste("ERROR at node",i+length(phylo$tip.label)))
			}	
		}
	#now replace 0s in mat with numbers
		for(i in 1:length(mat[,1])){
		if(mat[i,3]==0){
		mat[i,3]<-as.character(match(mat[i,2],phylo$tip.label))
	}}	
	mat<-as.data.frame(mat)
	if(is.null(maps.object)){
		colnames(mat)<-c("node.number","descendant.branch","next.node")
			geography.matrix<-list()
			for(i in 1:phylo$Nnode){ 
				var.list<-unlist(nat[[i]])
				len=length(var.list)
				int.mat<-matrix(nrow=len,ncol=len)
				rownames(int.mat)<-var.list
				colnames(int.mat)<-var.list
				geography.matrix[[i]]<-int.mat
			}
		}
		
	if(!is.null(maps.object)){
		if(is.list(maps.object)){
			if(flag==0){
				map<-.edge.mapped(phylo,maps.object)
				map.v<-map[match(mat[,3],map[,2]),3]
			}
			if(flag==1){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
				map<-.edge.mapped(old.phylo,maps.object) 
				map<-cbind(phylo$edge,map[,3])
				map.v<-map[match(mat[,3],map[,2]),3]
			}		
			mat<-data.frame(mat,map.v)
		}
		if(is.matrix(maps.object)){
			if(flag==0){
				map<-maps.object
				map.v<-map[match(mat[,3],map[,2]),3]
				}
			if(flag==1){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
				map<-cbind(phylo$edge,maps.object[,3])
				map.v<-map[match(mat[,3],map[,2]),3]
			}
			mat<-data.frame(mat,map.v)
		}
		colnames(mat)<-c("node.number","descendant.branch","next.node","region")
		geography.matrix<-list()
		for(i in 1:phylo$Nnode){ 
			var.list<-unlist(nat[[i]])
			len=length(var.list)
			int.mat<-matrix(nrow=len,ncol=len)
			rownames(int.mat)<-var.list
			colnames(int.mat)<-var.list
			diag(int.mat)<-1
			for(j in 1:length(var.list)){
				for(k in 1:length(var.list)){
					if((lower.tri(int.mat,diag=TRUE)[j,k]==TRUE)){
						int.mat[j,k]<-ifelse(mat[match(var.list[j],mat[,2]),4]==mat[match(var.list[k],mat[,2]),4],1,0)
					}
				}
			}
			int.mat[upper.tri(int.mat)==TRUE]<-t(int.mat)[upper.tri(t(int.mat))==TRUE]
			geography.matrix[[i]]<-int.mat
		}
	}						
		return(list(node.descendants=mat,geography.object=geography.matrix))
	}
	
