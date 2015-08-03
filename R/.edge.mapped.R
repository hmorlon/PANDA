.edge.mapped<-function(phylo,maps.object,round=TRUE){
	hold<-vector()
	if(round==TRUE){
	for(i in 1:length(maps.object)){
		if(length(names(maps.object[[i]]))==1){hold<-c(hold,names(maps.object[[i]]))}
		else{hold<-c(hold,names(which(maps.object[[i]]==max(maps.object[[i]]))))}		
	}}
	if(round==FALSE){
	stop("The current version does not support changes in biogeography along branches")
	}
	mat<-cbind(phylo$edge,hold,deparse.level=0)
	return(mat)
	}
