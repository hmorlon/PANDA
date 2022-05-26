mantel_test_nbpartners <-
function(network, tree_A, tree_B = NULL, method="Jaccard_binary", nperm=1000, correlation="Pearson"){
  
  if (!correlation %in% c("Pearson", "Spearman")) {stop("\"correlation\" must be among 'Pearson' or 'Spearman'.")}
  if (!is.numeric(nperm)) {stop("Please provide a numeric number of permutations (\"nperm\").")}
  
  if (!inherits(tree_A, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!is.null(tree_B)) {if (!inherits(tree_B, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}}
  
  if (is.null(method)) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', and 'UniFrac_unweighted'.")}
  if (method %in% c("GUniFrac", "UniFrac_unweighted", "PBLM", "PBLM_binary")) {if (is.null(tree_B)) stop("Please provide a phylogenetic tree \"tree_B\" for guild B.")}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "Bray-Curtis", "GUniFrac", "UniFrac_unweighted")) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', and 'UniFrac_unweighted'.")}
  
  if (nrow(network)<2){stop("Please provide a \"network\" with at least 2 species in clade B.")}
  if (ncol(network)<2){stop("Please provide a \"network\" with at least 2 species in clade A.")}
  
  
  # Only keep species with at least 1 interaction
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  
  # A in columns and B in rows
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  names(nb_A) <- "nb_A"
  names(nb_B) <- "nb_B"
  
  # Check names
  if (all(is.null(colnames(network)))|all(is.null(rownames(network)))) {stop("Please provide a \"network\" with row names and columns names matching the species names.")}
  
  if (!all(colnames(network) %in% tree_A$tip.label)){stop("Please provide a \"tree_A\" for all the species in clade A (the columns of the intercation network).")}
  if (only_A==FALSE) { if (!all(rownames(network) %in% tree_B$tip.label)){stop("Please provide a \"tree_B\" for all the species in clade B (the rows of the intercation network).")}}
  
  tree_A <- ape::drop.tip(tree_A,tip=tree_A$tip.label[which(!tree_A$tip.label %in% colnames(network))])
  if (only_A==FALSE) { tree_B <- ape::drop.tip(tree_B,tip=tree_B$tip.label[which(!tree_B$tip.label %in% rownames(network))])}
  
  
  if (!is.rooted(tree_A)){tree_A <- phytools::midpoint.root(tree_A) }
  if (only_A==FALSE) { if (!is.rooted(tree_B)){tree_A <- phytools::midpoint.root(tree_B) }}
  
  if (only_A==TRUE) { 
    network <- network[1:nrow(network),tree_A$tip.label]
  } else {
    network <- network[tree_B$tip.label,tree_A$tip.label]
  }
  
  
  
  
  compute_eco_dist <- function(network){
    # binary Jaccard distances
    if (method=="Jaccard_binary"){
      jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
      eco_A <- jaccard_A
    }
    
    # quantitative Jaccard distances
    if (method=="Jaccard_weighted"){
      jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
      eco_A <- jaccard_A
    }
    
    
    # Bray-Curtis dissimilarity 
    if (method=="Bray-Curtis"){
      bray_A <- as.matrix(vegan::vegdist(t(network), "bray", binary=F))
      eco_A <- bray_A
    }
    
    # Unifrac (generalized UniFrac, with alpha=0.5)
    if (method=="GUniFrac"){
      unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
      index=1
      eco_A <- unifrac_A$unifracs[,,index]
    }
    
    # Unifrac (unweighted UniFrac)
    if (method=="UniFrac_unweighted"){
      unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
      index=2
      eco_A <- unifrac_A$unifracs[,,index]
    }
    return(eco_A)
  }
  
  eco_A <- compute_eco_dist(network)
  
  # Perform Mantel test:
  
  # cophenetic distances
  cophe_A <- ape::cophenetic.phylo(tree_A)

  nb_A <- ncol(network)
  nb_B <- nrow(network)
  
  results <- c(as.integer(nb_A), as.integer(nb_B), NA, NA, NA, NA, NA, NA)
  names(results) <-  c("nb_A","nb_B","mantel_cor_A","pvalue_upper_A","pvalue_lower_A", "mantel_cor_B", "pvalue_upper_B", "pvalue_lower_B")
  
  if (length(unique(as.vector(cophe_A)))<3) {
    print("The phylogenetic distance matrix is composed of only 2 different values (because of polytomies?).")
    return(results)}
  if (length(unique(as.vector(eco_A)))<3) {
    print("The ecological distance matrix is composed of only 1 value (identical patterns of interactions across species?).")
    return(results)}
  
  
  if (correlation=="Pearson") {correlation="pearson"}
  if (correlation=="Spearman") {correlation="spearman"}
  
  original_correlation <- cor(as.vector(as.dist(eco_A)), as.vector(as.dist(cophe_A)), use = "everything", method = correlation)
  
  
  # Make randomizations

  random_correlation <- vector(mode = "numeric", length = nperm)
  vector_cophe_A <- as.vector(as.dist(cophe_A))
  
  for (i in 1:nperm){
    rand_network <- network
    
    for (k in 1:nb_A){
      rand_network[,k] <- sample(rand_network[,k])
    }
    
    eco_A_rand <- compute_eco_dist(rand_network)
    
    random_correlation[i] <- cor(as.vector(as.dist(eco_A_rand)), vector_cophe_A, use = "everything", method = correlation)
  
  }
  
  results <- c(original_correlation, min(c(length(which(c(random_correlation,original_correlation)>=original_correlation))/nperm,1)), min(c(length(which(c(random_correlation,original_correlation)<=original_correlation))/nperm, 1)))

  return(results)
}
