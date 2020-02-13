tree_change_given_switches <-
function(tree,switches){ 
  maxlen <- max(node.depth.edgelength(tree))
  loop=0
  ksi <- ncol(switches)
  n <- Ntip(tree)
  z <- switches[3,]
  output <- tree
  for (index in 1:ksi) { 
    i <- z[index]
    node_length <- node.depth.edgelength(output)
    i <- i - (maxlen - max(node_length)) 
    b <- switches[1,index]
    bs <-  switches[2,index]
    output <- tree_switch(output,i,b,bs,node_length)
    loop <- loop +1 }
  return(output)}
