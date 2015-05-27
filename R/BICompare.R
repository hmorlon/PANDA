#get BIC values for tree matrix and random matrix for t modalities
#returns assignation of nodes to the t clusters
#plot tree with t clusters

BICompare <- function(phylo,t){
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
			rDP = matrix(rnorm(length(phyloM),median(phyloM),3*sd(phyloM)))
			kmeansBIC(kmeans(rDP,t,algorithm="Hartigan-Wong"))->r
			kmeans(phyloM,t,algorithm="Hartigan-Wong")->q
			kmeansBIC(q)->p
rp<-cbind(p,r)
colnames(rp)<-c("tree BIC","control BIC")
res<-list("BIC_test"=rp,"clusters"=q$cluster)

col_edge<-rainbow(t)[res[[2]][phylo$edge[,2]]]
col_tip<-rainbow(t)[res[[2]][1:length(phylo$tip.label)]]
plot(phylo,edge.color=col_edge,tip.color=col_tip,type="fan")

return(res)


}
