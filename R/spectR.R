#plot spectral density
spectR <- function(phylo,method=c("standard")){
		
	if(method=="standard"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=F)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=1)
	
	#get eigengap
		abs(diff(m))->gaps
			as.matrix(gaps)->gapMat
				c(1:length(gapMat))->modalities
			cbind(modalities,gapMat)->gapMatCol
		subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))->eigenGap
					
		res<-list(eigenvalues=e$values, eigengap=eigenGap[,1])
	}
	
	if(method=="normal1"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=T)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=0)

		#get eigengap
		abs(diff(m))->gaps
			as.matrix(gaps)->gapMat
				c(1:length(gapMat))->modalities
			cbind(modalities,gapMat)->gapMatCol
		subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))->eigenGap
				
		res<-list(eigenvalues=e$values, eigengap=eigenGap[,1])
			
			}

	if(method=="normal2"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=F)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=0)

		#get eigengap
		abs(diff(m))->gaps
			as.matrix(gaps)->gapMat
				c(1:length(gapMat))->modalities
			cbind(modalities,gapMat)->gapMatCol
		subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))->eigenGap
		
		res<-list(eigenvalues=e$values, eigengap=eigenGap[,1])

	}
	class(res)	<- "spectR"
	return(res)					
}


