delineate_phylotypes <-
function(tree, thresh=97, sequences, method="pi"){	
  
  if (!inherits(tree, "phylo")) {stop(print("object \"phy\" is not of class \"phylo\".\n"))}
  if (!is.rooted(tree)) {stop(print("Please provide a rooted phylogeny"))}
  if (!all(tree$tip.label %in% rownames(sequences))) {stop(print("Please provide a nucleotidic alignment with sequence names matching the tip labels of the phylogenetic tree"))}
  
  if (!method %in% c("pi","theta", "mean")){stop(print("Please provide a method among 'pi' or 'theta'"))}
  
  N.tip<-Ntip(tree)
  bins<-matrix(0,nrow=N.tip,ncol=2)
  rownames(bins) <- tree$tip.label
  bins[,2] <- tree$tip.label
  
  it_OTU <- 0 
  if (N.tip==1){
    print("single")
    sequence<-tree$tip.label
    bins[sequence,1]<-1
  } else {
    for (i in (N.tip+1):(N.tip+Nnode(tree))){
    extracted_tree <- extract.clade(tree, node = i)
        if (all(bins[extracted_tree$tip.label,1]==0)){
        
          if (method=="pi") mean_dist <- pi_estimator(sequences[which(rownames(sequences) %in% extracted_tree$tip.label),])
          if (method=="theta") mean_dist <- theta_estimator(sequences[which(rownames(sequences) %in% extracted_tree$tip.label),])
          if (method=="mean"){
            extracted_sequences <- as.DNAbin(sequences[which(rownames(sequences) %in% extracted_tree$tip.label),])
            distances_matrix <- dist.dna(x=extracted_sequences, pairwise.deletion=TRUE, as.matrix=FALSE, model="raw") 
            mean_dist <- mean(as.vector(distances_matrix), na.rm=TRUE)
          }
          
          if (mean_dist<(100-thresh)/100){
            it_OTU <- it_OTU + 1
            bins[extracted_tree$tip.label,1] <- it_OTU
            
            extracted_sequences <- sequences[which(rownames(sequences) %in% extracted_tree$tip.label),]
            bins[extracted_tree$tip.label,2] <- rownames(extracted_sequences)[which.max(sapply(1:nrow(extracted_sequences), function(k) length(which(!extracted_sequences[k,] %in% c("-", "N", "n")))))]
            
          }
        }
      }
  }
  bins <- data.frame(bins, stringsAsFactors = FALSE)
  colnames(bins) <- c("phylotype", "representative_sequence")
  
  print(paste0("Number of phylotypes (including singletons): ", it_OTU+length(which(bins$phylotype=="0"))))
  print(paste0("Number of phylotypes (excluding singletons): ", it_OTU))
  print(paste0("Number of singletons: ", length(which(bins$phylotype=="0"))))
  return(bins)
}
