JSDtree_cluster <- function(JSDtree,alpha=0.8)
{

#plot heatmap
heatmap(JSDtree,symm=T)

#plot hierarchical clustering with bootstrap support
quartz()
	clustersHierarchy <- pvclust(JSDtree)
	plot(clustersHierarchy,cex=0.3)
	pvrect(clustersHierarchy,alpha=alpha)
	
#cluster JSD matrix on medoids
clustersMedoid <- pamk(JSDtree)
clustersMedoidSupport <- pam(JSDtree,clustersMedoid$nc)


#print clustersMedoid
res	<- list(clusters=clustersMedoid$nc, cluster_assignments=clustersMedoid[[1]][[3]],cluster_support=clustersMedoidSupport[[7]])
return(res)

  }

