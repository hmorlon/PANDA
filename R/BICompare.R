#get BIC values for tree matrix and random matrix for 1:t modalities
BICompare <- function(phylo,t,method=c("gaussian","poisson")){
	#get lLk, BIC
	kmeansBIC <- function(fit){
		m <- ncol(fit$centers)
		n <- length(fit$cluster)
		k <- nrow(fit$centers)
		D <- fit$tot.withinss
	return(data.frame(BIC = D + log(n)*m*k))
	}
		phyloM <- as.matrix(dist.nodes(phylo))	
		DP <- c()
		if(method=="gaussian"){
			rDP = matrix(rnorm(length(phyloM),median(phyloM),							3*sd(phyloM)))#gaussian
				}
		matrix()->r
		matrix()->p				
		for(i in c(1:t)){	
			kmeansBIC(kmeans(rDP,t,algorithm="Hartigan-Wong"))->r[i]
			kmeansBIC(kmeans(phyloM,t,algorithm="Hartigan-Wong"))->p[i]
			}

		if(method=="gaussian"){
			rDP = matrix(rpois(length(phyloM),median(phyloM)))#poisson
				}
		for(i in c(1:t)){	
			kmeansBIC(kmeans(rDP,t,algorithm="Lloyd"))->r[i]
			kmeansBIC(kmeans(phyloM,t,algorithm="Lloyd"))->p[i]
			}
				
rp<-cbind(p,r)
colnames(rp)<-c("tree BIC","control BIC")
return(rp)
}
