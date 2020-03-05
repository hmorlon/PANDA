tree_change_temporal <-
function(tree,ksi,maxlen,mean,sd){ 
  if (sd<0) stop(print("Please provide a realistic standart deviation for the temporal distribution of switches"))
  if (all(rnorm(10,mean,sd)<0)|all(rnorm(10,mean,sd)>maxlen)) stop(print("Please provide a realistic mean for the temporal distribution of switches"))
  
  loop=0
  n <- Ntip(tree)
  while (loop<ksi){
    loop=0
    switches <- matrix(0,nrow=3,ncol=ksi)
    z <- c()
    while(length(z)<ksi){position <- rnorm(1,mean,sd)      
    if ((position>0)&(position<maxlen)) {z[length(z)+1] <- position}}
    z <- sort(z)
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
        loop <- loop +1 }}}
  }
  return(list(output,switches))}
