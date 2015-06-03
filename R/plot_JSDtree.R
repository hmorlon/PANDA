plot_JSDtree <- function(JSDtree)
{
  if (!inherits(JSDtree, "JSDtree"))
      stop("object \"JSDtree\" is not of class \"JSDtree\"")

heatmap(JSDtree,symm=T)
  
  }

