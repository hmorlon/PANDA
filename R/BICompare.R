#get BIC values for tree matrix and random matrix for t modalities
#returns assignation of nodes to the t clusters
#plot tree with t clusters

BICompare <- function(phylo,t,meth=c("ultrametric")){
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
		c()->q	
	##if tree is non-ultrametric	
		if(meth=="non-ultrametric"){	
			rDP = as.matrix(
					dist.nodes(
				rtree(n=length(phylo$tip.label),
			br=runif(100,min=min(phylo$edge.length),
		max=max(phylo$edge.length)))))
			}	
	##if tree is ultrametric	
		if(meth=="ultrametric"){
			rDP = as.matrix(
					dist.nodes(
				rcoal(n=length(phylo$tip.label),
					br=runif(100,0,max(branching.times(phylo)/20)))))
			}	
	##get BICs for tree and control
		kmeansBIC(kmeans(rDP,t,algorithm="Hartigan-Wong"))->r
			kmeans(phyloM,t,algorithm="Hartigan-Wong")->q
				kmeansBIC(q)->p
rp<-cbind(p,r)
colnames(rp)<-c("tree BIC","control BIC")
res<-list("BIC_test"=rp,"clusters"=q$cluster,"BSS/TSS"=q$betweenss/q$totss)

class(res)	<- "BICompare"
return(res)

}
