sim_duplications_tree <-
function(host_tree, delta){
  
  host_tree$tip.label <- paste0(host_tree$tip.label, "-seq1")
  
  list_duplicated_tree <- list()
  list_tips <- c()
  list_positions <- c()
  
  count_duplication <- 1
  
  
  j=0 # which tree
  
  while(j<count_duplication){
    
    if (j==0) {tree <- host_tree
    root_age <- 0 # can change this if the duplication of the tree is possible at the MRCA
    }else{ tree <- list_duplicated_tree[[j]]
    root_age <- list_positions[j]
    }
    
    if (class(tree)=="character"){  # singleton
      
      branch_length <- list_positions[j]
      
      time <- rexp(n=1, rate = delta)
      
      while(time<branch_length){
        
        count_duplication <- count_duplication + 1
        
        extracted_tree <- tree
        list_tips <- rbind(list_tips, c(tree,tree))
        
        list_positions <- c(list_positions, branch_length-time) 
        
        extracted_tree <- gsub(pattern = "-seq.*", replacement = paste0("-seq",count_duplication), extracted_tree)
        
        list_duplicated_tree <- c(list_duplicated_tree, list(extracted_tree))
        
        branch_length <- branch_length-time
        time <- rexp(n=1, rate = delta)
      }
      
    }else{
      
      if (root_age>0){
        # duplicate again the tree? 
        
        branch_length <- list_positions[j]
        
        time <- rexp(n=1, rate = delta)
        
        while(time<branch_length){
          
          count_duplication <- count_duplication + 1
          
          extracted_tree <- tree
          list_tips <- rbind(list_tips, extracted_tree$tip.label[c(1,Ntip(extracted_tree))])
          
          list_positions <- c(list_positions, branch_length-time) 
          
          extracted_tree$tip.label <- gsub(pattern = "-seq.*", replacement = paste0("-seq",count_duplication), extracted_tree$tip.label)
          
          list_duplicated_tree <- c(list_duplicated_tree, list(extracted_tree))
          
          branch_length <- branch_length-time
          time <- rexp(n=1, rate = delta)
        }
        
      }
      
      for (i in 1:nrow(tree$edge)){ # for each branch
        
        branch_length <- tree$edge.length[i]
        
        time <- rexp(n=1, rate = delta)
        
        while(time<branch_length){
          # make duplication
          
          count_duplication <- count_duplication + 1
          
          node <- tree$edge[i,2]
          
          if (node>Ntip(tree)){
            extracted_tree <- extract.clade(tree, node=node)
            list_tips <- rbind(list_tips, extracted_tree$tip.label[c(1,Ntip(extracted_tree))])
            
            extracted_tree$tip.label <- gsub(pattern = "-seq.*", replacement = paste0("-seq",count_duplication), extracted_tree$tip.label)
            
          }else{
            extracted_tree <- tree$tip.label[node]
            extracted_tree <- gsub(pattern = "-seq.*", replacement = paste0("-seq",count_duplication), extracted_tree)
            list_tips <- rbind(list_tips, c(tree$tip.label[node],tree$tip.label[node]))
          }
          
          list_positions <- c(list_positions, branch_length-time) 

          list_duplicated_tree <- c(list_duplicated_tree, list(extracted_tree))
          
          branch_length <- branch_length-time
          time <- rexp(n=1, rate = delta)
          
        }
      }
    }
    j=j+1
  }
  
  if(count_duplication>1){
    for (i in 1:(count_duplication-1)){
      
      duplicated_tree <- list_duplicated_tree[[i]]
      
      position <- list_positions[i]
      tips <- list_tips[i,]

      if (class(duplicated_tree)=="character"){
        tree <- pbtree(n=2)
        tree$tip.label <- c(duplicated_tree, "remove")
        tree$edge.length <- c(position,position)
        duplicated_tree <- tree
        
        host_tree <- bind.tree(x=host_tree, y=duplicated_tree, where = which(host_tree$tip.label==unique(tips)), position = position, interactive = FALSE)
        host_tree <- drop.tip(host_tree, tip="remove")
        
      }else{
        duplicated_tree$root.edge <- position
        host_tree <- bind.tree(x=host_tree, y=duplicated_tree, where = getMRCA(host_tree, tip = tips), position = position, interactive = FALSE)
      }
    }}
  
  return(c(list(host_tree), count_duplication-1))
}
