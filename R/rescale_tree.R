rescale_tree <-
function(tree){
  tree$edge.length <- tree$edge.length/sum(tree$edge.length)
  return(tree)}
