\name{plot_JSDtree}
\alias{plot_JSDtree}

\title{
Heatmap comparison of phylogenies
}
\description{
Plots a hierarchical cluster with bootstrap support and heatmap based on the Jensen-Shannon distance between more than 10 trees.
}
\usage{
plot_JSDtree(JSDtree,alpha)
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
a plot of heatmap with hierarchical cluster and plot of hierarchical cluster with approximately unbiased (AU) and boostrap probability (BP) percentage support.
}

\author{E Lewitus}

\seealso{
\code{\link{JSDtree}}
}

\examples{
trees<-TESS::sim.globalBiDe.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
res<-JSDtree(trees)
#plot_JSDtree(res,alpha=0.8)
}
