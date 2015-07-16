\name{JSDtree_cluster}
\alias{JSDtree_cluster}

\title{
Heatmap comparison, hierarchical clustering, and k-medoid clustering of phylogenies
}
\description{
Plots a hierarchical cluster with bootstrap support, plots a heatmap, and prints k-medoid cluster assignments and support based on the Jensen-Shannon distance between more than 10 trees.
}
\usage{
JSDtree_cluster(JSDtree,alpha)
}

\arguments{
  \item{JSDtree}{
an object of class 'JSDtree', output of the 'JSDtree' function
}
 \item{alpha}{
   	the confidence value for demarcating clusters in the hierarchical clustering plot
	}
	}


\value{
	\item{hierarchical plot}{a plot of hierarchical cluster with approximately unbiased (AU) and boostrap probability (BP) percentage support}
	\item{heatmap}{a heatmap plot of Jensen-Shannon distances}
	\item{clusters}{the optimal number of clusters around medoids (see pamk documentation)}
	\item{cluster_assignments}{assignments of trees to clusters} 
	\item{cluster_support}{cluster to which each tree belongs, the most similar cluster, and the silhouette width of the observation; average silhouette width for each cluster; average silhouette width across all clusters (see silhouette documentation)}
}

\author{E Lewitus}

\seealso{
\code{\link{JSDtree}}
}

\examples{
trees<-TESS::sim.globalBiDe.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
res<-JSDtree(trees)
JSDtree_cluster(res,alpha=0.9)
}
