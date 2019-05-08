

#get Jensen-Shannon divergence	
trait_JSDcluster <- function(phylo,mat,plot=F){
	dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x*log(x/y))
	JSD <- function(x,y) sqrt(0.5*KLD(x,(x+y)/2)+0.5*KLD(y,(x+y)/2))
		matrixColSize <- length(colnames(inMatrix))
		matrixRowSize <- length(rownames(inMatrix))
		colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize) 
	  	inMatrix = apply(inMatrix,1:2,
  		function(x) ifelse (x==0,pseudocount,x))
	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
 }
#gaussian kernel convolution	 
dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
  if(has.na <- any(is.na(x))) {
    na.omit(x)->x
    if(length(x) == 0)
        stop('too infinite.')
  }
  	kernelG<-function(x, mean=0, sd=1) 
		dnorm(x, mean = mean, sd = sd)
  x<-log(x)		
  sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}

#take square-root of Jensen-Shannon divergence
JSDist <- function(x,y) sqrt(dist.JSD(x,y))


x<-lapply(1:dim(mat)[2],function(l){
	trait_spectR(phylo,mat[,l])$eigenvalues
	} )

d<-c()
Ds<-c()
	for(i in 1:length(x)){
			d[[i]]<-dens(x[[i]])$x
			Ds<-as.data.frame(cbind(Ds,d[[i]]))
			}
			colnames(Ds)<-colnames(mat)
	
	#compute divergence matrix
	if(min(Ds)<0){JSD<-as.matrix(JSDist(Ds-min(Ds)+1e-9))}
	else{JSD<-as.matrix(JSDist(Ds))}	
	
	#cluster on k-medoids
	clustersMedoid <- pamk(JSD)
	
	#write table for divergence matrix
	write.table(JSD,file='JSD_divergenceMatrix.txt')
	
	#write table for clusters with silhouette widths
	write.table(clustersMedoid[[1]]$silinfo$widths[,c(1,3)],file='JSD_divergenceMatrix_clusters.txt')
	
	if(plot==T){
	#plot heatmap
	heatmap(JSD,symm=T)
		dev.new()
	#plot hierarchical clustering with bootstrap support
	clustersHierarchy <- pvclust(JSD,r=seq(0.5,1.5,0.2))
	plot(clustersHierarchy,cex=0.3)
	pvrect(clustersHierarchy,alpha=0.9)
	return(clustersMedoid[[1]]$clusinfo)
	}
	else{return(clustersMedoid[[1]]$clusinfo)}
}	
	
