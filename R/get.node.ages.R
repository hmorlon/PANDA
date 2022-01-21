get.node.ages <-  function(nodes, phylo. = phylo){
  ALL_nodes_ages <- as.data.frame(apply(data.frame(nodesID=names(branching.times(phylo)),ages=branching.times(phylo)), 2, as.numeric))
  nodes_ages_selected <- sort(ALL_nodes_ages$ages[ALL_nodes_ages$nodesID %in% nodes])
  return(nodes_ages_selected)
}