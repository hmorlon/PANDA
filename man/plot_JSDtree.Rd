\name{plot_JSDtree}
\alias{plot_JSDtree}

\title{
Heatmap comparison of phylogenies
}
\description{
Plots a hierarchical cluster and heatmap based on the Jensen-Shannon distance between more than 10 trees.
}
\usage{
plot_JSDtree(JSDtree)
}

\arguments{
  \item{JSDtree}{
an object of class 'JSDtree', output of the 'JSDtree' function
}
	}


\value{
a plot of heatmap with hierarchical cluster
}

\author{E Lewitus}

\seealso{
\code{\link{JSDtree}}
}

\examples{
trees<-TESS::sim.globalBiDe.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
res<-JSDtree(trees,alpha=0.8)
#plot_JSDtree(res)
}
