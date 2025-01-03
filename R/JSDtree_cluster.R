JSDtree_cluster <- function(JSDtree,alpha=0.9,draw=TRUE)
{


#cluster JSD matrix on medoids
maxCluster<-dim(JSDtree)[1]-1
clustersMedoid <- pamk(JSDtree,krange=1:maxCluster)
clustersMedoidSupport <- pam(JSDtree,clustersMedoid$nc)

if(draw == TRUE){
#plot heatmap
heatmap(JSDtree,symm=TRUE)

#plot hierarchical clustering with bootstrap support
	dev.new()
	clustersHierarchy <- pvclust(JSDtree)
	plot(clustersHierarchy,cex=0.3)
	pvrect(clustersHierarchy,alpha=alpha)
}
else{}

#print clustersMedoid
res	<- list(clusters=clustersMedoid$nc, cluster_assignments=clustersMedoid[[1]][[3]],cluster_support=clustersMedoidSupport[[7]])
return(res)

  }

