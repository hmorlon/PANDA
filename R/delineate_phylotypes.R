delineate_phylotypes <-
function(tree, thresh=97, sequences, method="pi"){	
  
  if (!inherits(tree, "phylo")) {stop(print("object \"phy\" is not of class \"phylo\".\n"))}
  if (is.rooted(tree)) {stop(print("Please provide a rooted phylogeny"))}
  if (!all(tree$tip.label %in% rownames(sequences))) {stop(print("Please provide a nucleotidic alignment with sequence names matching the tip labels of the phylogenetic tree"))}
  
  if (method %in% c("pi","theta", "mean")){stop(print("Please provide a method among 'pi' or 'theta'"))}
  
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
        
          if (method=="pi") mean_dist <- pi_estimator(sequences[which(rownames(sequences) %in% extracted_tree$tip.label),])
          if (method=="theta") mean_dist <- theta_estimator(sequences[which(rownames(sequences) %in% extracted_tree$tip.label),])
          if (method=="mean"){
            extracted_sequences <- as.DNAbin(sequences[which(rownames(sequences) %in% extracted_tree$tip.label),])
            distances_matrix <- dist.dna(x=extracted_sequences, pairwise.deletion=T, as.matrix=F, model="raw") 
            mean_dist <- mean(as.vector(distances_matrix), na.rm=TRUE)
          }
          
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
