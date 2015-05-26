#get BIC values for tree matrix and random matrix for 1:t modalities
BICompare <- function(phylo,t,method=c("gaussian","poisson")){
	options(warn=-1)
	#get lLk, BIC
	kmeansBIC <- function(fit){
		m <- ncol(fit$centers)
		n <- length(fit$cluster)
		k <- nrow(fit$centers)
		D <- fit$tot.withinss
	return(data.frame(BIC = D + log(n)*m*k))
	}
		phyloM <- as.matrix(dist.nodes(phylo))	
		rDP <- c()
		c()->r
		c()->p	
		if(method=="gaussian"){
			rDP = matrix(rnorm(length(phyloM),median(phyloM),							3*sd(phyloM)))	

			kmeansBIC(kmeans(rDP,t,algorithm="Hartigan-Wong"))->r
			kmeansBIC(kmeans(phyloM,t,algorithm="Hartigan-Wong"))->p
	}
		if(method=="poisson"){
			rDP = matrix(rpois(length(phyloM),median(phyloM)))

			kmeansBIC(kmeans(rDP,t,algorithm="Lloyd"))->r
			kmeansBIC(kmeans(phyloM,t,algorithm="Lloyd"))->p
	}			
rp<-cbind(p,r)
colnames(rp)<-c("tree BIC","control BIC")
return(rp)

}
