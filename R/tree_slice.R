tree_slice <-
function (tree, slice = NULL) {
  tree <- reorder(tree)
  H <- nodeHeights(tree)
  edges <- which(H[, 2] > slice & H[, 1] < slice)
  nodes <- tree$edge[edges, 2]
  trees <- list()
  class(trees) <- "multiPhylo"
  for (i in 1:length(nodes)) {
    if (nodes[i] > Ntip(tree)) {
      trees[[i]] <- extract.clade(tree, node = nodes[i])
      trees[[i]]$root.edge <- H[which(tree$edge[, 2] == 
                                        nodes[i]), 2] - slice
    }
    else {
      z <- list(edge = matrix(c(2, 1), 1, 2), edge.length = H[which(tree$edge[, 
                                                                              2] == nodes[i]), 2] - slice, tip.label = tree$tip.label[nodes[i]], 
                Nnode = 1L)
      class(z) <- "phylo"
      trees[[i]] <- z
    }
  }
  return(trees)
}
