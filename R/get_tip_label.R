get_tip_label <-
function(tree,branch){
  tip <- tree$edge[branch,2]
  while(!all(tip<(Ntip(tree)+1))){
    new_tip <- c()
    for(i in tip){if (i>Ntip(tree)) {new_tip <- c(new_tip,tree$edge[which(tree$edge[,1]==i),2]) 
    }else{ new_tip <- c(new_tip,i)} }
    tip <- new_tip}
  return(tree$tip.label[tip])
}
