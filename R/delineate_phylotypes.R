delineate_phylotypes <-
function(tree, thresh, sequences){	
  
  if (!inherits(tree, "phylo")) {stop(print("object \"phy\" is not of class \"phylo\".\n"))}
  if (is.rooted(tree)) {stop(print("Please provide a rooted phylogeny"))}
  if (all(tree$tip.label %in% rownames(sequences))) {stop(print("Please provide a nucleotidic alignment with sequence names matching the tip labels of the phylogenetic tree"))}
  
  N.tip<-Ntip(tree)
  bins<-matrix(0,nrow=N.tip,ncol=2)
  rownames(bins)<-tree$tip.label
  
  if(N.tip==1){
    print("single")
    sequence<-tree$tip.label
    bins[sequence,1]<-1
  } else {
    it_OTU <- 0 
    for (i in (N.tip+1):(N.tip+Nnode(tree))){
    extracted_tree <- extract.clade(tree, node = i)
      if (Ntip(extracted_tree)<Ntip(tree)/8){
        if (all(bins[extracted_tree$tip.label,1]==0)){
        
          # Pi estimate
          mean_dist <- pi_estimator(sequences[which(rownames(sequences) %in% extracted_tree$tip.label),])
        
          if (mean_dist<(100-thresh)/100){
            it_OTU <- it_OTU + 1
            print(it_OTU)
            bins[extracted_tree$tip.label,1] <- it_OTU
            bins[extracted_tree$tip.label,2] <- tree$node.label[i-N.tip]
          }
        }
      }
    }
    }
  return(bins)
}
