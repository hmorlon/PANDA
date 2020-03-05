tree_change <-
function(tree,ksi,maxlen){ 
  loop=0
  n <- Ntip(tree)
  while (loop<ksi){
    loop=0
    switches <- matrix(0,nrow=3,ncol=ksi)
    b = sample(x=seq(1:(2*(n-1))), size=ksi, prob=tree$edge.length, replace = T) 
    l_a1=node.depth.edgelength(tree)[tree$edge[b,1]]
    l_b1=node.depth.edgelength(tree)[tree$edge[b,2]]
    z <- sort(runif(ksi, l_a1, l_b1))
    switches[3,] <- z 
    output <- tree
    for (index in 1:ksi) { 
      i <- z[index]
      node_length <- node.depth.edgelength(output)
      i <- i - (maxlen - max(node_length)) 
      branches <- intersect(which(node_length[output$edge[,1]]<i),which(node_length[output$edge[,2]]>i))
      br <- length(branches)
      b <-  branches[sample(br, size=1)]
      bs <- branches[branches!=b][sample(br-1, size=1)] 
      if (is.integer(bs)){ if (length(bs)==1){
        switches[1,index] <- b
        switches[2,index] <- bs 
        output <- tree_switch(output,i,b,bs,node_length)
        loop <- loop +1 }}}}
  return(list(output,switches))}
