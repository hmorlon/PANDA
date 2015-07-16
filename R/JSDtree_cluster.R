JSDtree_cluster <- function(JSDtree,alpha)
{
  if (!inherits(JSDtree, "JSDtree"))
      stop("object \"JSDtree\" is not of class \"JSDtree\"")

#cluster JSD matrix on medoids
clustersMedoid <- pamk(JSDtree$JSD)
	clustersMedoidSupport <- pam(JSDtree$JSD,clustersMedoid$nc)

#plot heatmap
heatmap(JSDtree$JSD,symm=T)

#plot hierarchical clustering with bootstrap support
quartz()
	clustersHierarchy <- pvclust(JSDtree$JSD)
	plot(clustersHierarchy,cex=0.5)
		pvrect(clustersHierarchy,alpha=alpha)

#print clustersMedoid
res	<- list(clusters=clustersMedoid$nc, cluster_assignments=clustersMedoid[[1]][[3]],cluster_support=clustersMedoidSupport[7])
return(res)

  }

