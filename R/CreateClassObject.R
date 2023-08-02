CreateClassObject<-function(map,rnd=5,return.mat=FALSE){
  if(any(grepl("___",map$tip.label))){stop("script will not work with '___' in tip labels; remove extra underscores")}
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
  nodeDist<-vector(mode = "numeric", length = map$Nnode)
  totlen<-length(map$tip.label)
  root <-totlen  + 1
  heights<-nodeHeights(map)
  for (i in 1:dim(map$edge)[1]){
    nodeDist[[map$edge[i, 1] - totlen]] <- heights[i]
  }
  nodeDist<-c(nodeDist,max(heights))
  nodeDiff<-diff(nodeDist)
  flag=0
  if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
    node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = map$Nnode))
    node.order<-node.order+totlen
    old.edge<-map$edge
    old.map<-map
    map$edge[,1]<-node.order
    for(j in 1:length(map$edge[,2])){
      if(map$edge[j,2]>totlen){
        #match number order in old edge
        #lookup value in new edge
        #replace with value
        map$edge[j,2]<-map$edge[,1][match(map$edge[j,2],old.edge[,1])]
      }
    }
    nodeDist<-vector()
    for (i in 1:dim(map$edge)[1]){
      nodeDist[[map$edge[i, 1] - totlen]] <- heights[i]
    }
    nodeDist<-c(nodeDist,max(heights))
    nodeDiff<-diff(nodeDist)
    flag=1
  }
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:map$Nnode){
    other<-map$edge[map$edge[,1]==i+totlen, 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+totlen
      if(b>totlen){
        counter_three_letters <- counter_three_letters + 1
        int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
        int[3]<-b
      } else {
        int[2]<-map$tip.label[[b]]
        int[3]<-0 
      }
      mat<-rbind(mat,int)
    }
  }
  if(geiger::is.phylo(map)){
    for(i in 1:length(mat[,1])){
      if(mat[i,3]==0){
        mat[i,3]<-as.character(match(mat[i,2],map$tip.label))
      }}	
    maps.object<-map$maps
    brchange<-vapply(maps.object,function(x) any(length(names(x))>1),1)
    if(sum(brchange)!=0){
      newtimes<-vector()
      for(i in 1:length(brchange)){
        if(brchange[i]!=0){
          #look up branch in stochastic map
          intlen<-length(maps.object[[i]])
          inttime<-vector()
          for(j in 1:(intlen-1)){
            inttime<-c(inttime,maps.object[[i]][j])
            nt<-as.numeric(heights[i,1]+sum(inttime))
            newtimes<-c(newtimes,nt)
          }
        }
      }
      old.Dist<-nodeDist
      old.Diff<-nodeDiff
      nodeDist<-sort(c(nodeDist,newtimes))
      nodeDiff<-diff(nodeDist)
    }
	if(any(nodeDiff<= (2*(10^-rnd)))){stop("CreateClassObject.R: potential rounding error, two time bins very similar")}
	nat<-list()
	nodecount=1
	for(i in 1:length(nodeDiff)){
		if(i==1){
			hold.m<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
			int<-dim(hold.m)[1]
			for(m in 1:int){
				hold.m[m,2]<-names(maps.object[[which(map$edge[,2]==as.numeric(hold.m[m,2]))]])[1]
			}
			nat[[i]]<-hold.m
			} else {
			if(nodeDist[i]%in%old.Dist){ #a new species appears, this is a node
				P<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
				hold.m<-rbind(P[as.numeric(P[,2])<=totlen,],P[as.numeric(P[,2])>(totlen+nodecount),])
				int<-dim(hold.m)[1]
				for(m in 1:int){
					iden<-which(map$edge[,2]==as.numeric(hold.m[m,2]))
					if(brchange[iden]==0 || !(hold.m[m,1]%in%nat[[i-1]][,1])){
						hold.m[m,2]<-names(maps.object[[iden]])[1]
					} else{
						num=1
						#while(round(nodeDist[i+1]-old.Dist[map$edge[iden,1]-totlen],6)>round(sum(maps.object[[iden]][1:num]),6)){
						while((round(nodeDist[i+1]-old.Dist[map$edge[iden,1]-totlen],rnd)-round(sum(maps.object[[iden]][1:num]),rnd))>(2*(10^-rnd))){
						num=num+1}
						hold.m[m,2]<-names(maps.object[[iden]])[num]
					}
				}
				nat[[i]]<-hold.m
			} else{ # this is a branch change
				P<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
				hold.m<-rbind(P[as.numeric(P[,2])<=totlen,],P[as.numeric(P[,2])>(totlen+nodecount),])
				int<-dim(hold.m)[1]
				for(m in 1:int){
					iden<-which(map$edge[,2]==as.numeric(hold.m[m,2]))
					if(brchange[iden]==0 | !(hold.m[m,1]%in%nat[[i-1]][,1])){
						hold.m[m,2]<-names(maps.object[[iden]])[1]
					} else{					
						num=1
						#while(round(nodeDist[i+1]-old.Dist[map$edge[iden,1]-totlen],5)>round(sum(maps.object[[iden]][1:num]),5)){
						while((round(nodeDist[i+1]-old.Dist[map$edge[iden,1]-totlen],rnd)-round(sum(maps.object[[iden]][1:num]),rnd))>(2*(10^-rnd))){
						num=num+1}
						hold.m[m,2]<-names(maps.object[[iden]])[num]
					}
				}
			nat[[i]]<-hold.m
			}
		}
	if(nodeDist[i+1]%in%old.Dist){nodecount=nodecount+1}
	}	
    }
    

  if(return.mat){
	mat[which(mat[,3]==0),3]<-match(mat[which(mat[,3]==0),2],map$tip.label)
	return(list(class.object=nat,times=round(nodeDist,8),spans=nodeDiff,mat=mat,phylo=map))
  }else{
  return(list(class.object=nat,times=round(nodeDist,8),spans=nodeDiff))
  }
}