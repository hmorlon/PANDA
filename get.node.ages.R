get.node.ages <-  function(nodes, phy. = phy){
  ALL_nodes_ages <- as.data.frame(apply(data.frame(nodesID=names(branching.times(phy)),ages=branching.times(phy)), 2, as.numeric))
  nodes_ages_selected <- sort(ALL_nodes_ages$ages[ALL_nodes_ages$nodesID %in% nodes])
  return(nodes_ages_selected)
}