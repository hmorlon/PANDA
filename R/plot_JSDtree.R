plot_JSDtree <- function(JSDtree,alpha)
{
  if (!inherits(JSDtree, "JSDtree"))
      stop("object \"JSDtree\" is not of class \"JSDtree\"")

#plot heatmap
heatmap(JSDtree$JSD,symm=T)

#plot hierarchical clustering with bootstrap support
quartz()
	clustersHierarchy <- pvclust(JSDtree$JSD)
	plot(clustersHierarchy,cex=0.5)
		pvrect(clustersHierarchy,alpha=alpha)

  }

