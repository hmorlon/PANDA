tree_switch <-
function(tree,z1,b,bs,node_length){
  r <- Ntip(tree)+1
  a1 <- tree$edge[b,1]
  b1 <- tree$edge[b,2]
  a2 <- tree$edge[bs,1]
  b2 <- tree$edge[bs,2]
  n_tree <- tree
  # ii) (a2 - b2) -> (a2' - b2)
  n_tree$edge.length[which(tree$edge[,2]==b2)] <- node_length[b2] - z1
  # iv) (a1 - b1) -> (a2' - b1)
  n_tree$edge[b,1] <- a2
  n_tree$edge.length[b] <- node_length[b1] - z1
  if ((a1==a2)&(a2!=r)){ # v) (a0 - a2) -> (a0 - a2')
    index <- which(tree$edge[,2]==a2)
    n_tree$edge.length[index] <- z1 - node_length[tree$edge[index,1]]}
  if ((a1!=a2)&(a2!=r)){ # iii) (a0 - a2) -> (a1 - a2') and (a2 - b3) -> (a0 - b3)
    index <- which(tree$edge[,2]==a2)
    n_tree$edge[index,1] <- a1
    n_tree$edge[index,2] <- a2
    n_tree$edge.length[index] <- z1 - node_length[a1]
    b3 <- setdiff(tree$edge[which(tree$edge[,1]==a2),2],b2)
    index <- which(tree$edge[,2]==b3) #bt
    a0 <- tree$edge[which(tree$edge[,2]==a2),1]
    n_tree$edge[index,1] <- a0
    n_tree$edge.length[index] <- node_length[b3]-node_length[a0] 
  }
  if((a1!=a2)&(a2==r)){ # (a2 - b3) -> (a1 - a2')  and Re-rooting in b3 (index a2 <-> b3)
    b3 <- setdiff(tree$edge[which(tree$edge[,1]==a2),2],b2)
    index <- which(tree$edge[,2]==b3) #bt
    n_tree$edge[index,1] <- a1
    n_tree$edge[index,2] <- a2
    n_tree$edge.length[index] <- z1-node_length[a1]
    tree <- n_tree
    tree$edge[which(n_tree$edge[,1]==b3),1] <- a2
    tree$edge[which(n_tree$edge[,2]==a2),2] <- b3
    tree$edge[which(n_tree$edge[,1]==a2),1] <- b3
    n_tree <- tree
  }
  return(n_tree)
}
