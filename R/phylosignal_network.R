phylosignal_network <-
function(network, tree_A, tree_B=NULL, method = "Jaccard_weighted", nperm = 10000, correlation = "Pearson", only_A = FALSE){

  if (is.null(tree_B)) {only_A <- TRUE} 
  
  if (!inherits(tree_A, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!is.null(tree_B)) {if (!inherits(tree_B, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}}
  
  if (is.null(method)) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'GUniFrac', 'UniFrac_unweighted', 'PBLM', 'PBLM_binary', and 'degree'.")}
  if (method %in% c("GUniFrac", "UniFrac_unweighted", "PBLM", "PBLM_binary")) {if (is.null(tree_B)) stop("Please provide a phylogenetic tree \"tree_B\" for guild B.")}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "GUniFrac", "UniFrac_unweighted", "PBLM", "PBLM_binary", "degree")) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'GUniFrac', 'UniFrac_unweighted', 'PBLM', 'PBLM_binary', and 'degree'.")}
  
  if (!correlation %in% c("Pearson", "Spearman", "Kendall")) {stop("Please pick a \"correlation\" among Pearson, Spearman, and Kendall.")}
  
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
  
  tree_A <- drop.tip(tree_A,tip=tree_A$tip.label[which(!tree_A$tip.label %in% colnames(network))])
  if (only_A==FALSE) { tree_B <- drop.tip(tree_B,tip=tree_B$tip.label[which(!tree_B$tip.label %in% rownames(network))])}
  

  if (!is.rooted(tree_A)){tree_A <- midpoint.root(tree_A) }
  if (only_A==FALSE) { if (!is.rooted(tree_B)){tree_A <- midpoint.root(tree_B) }}
  
  if (only_A==TRUE) { 
    network <- network[1:nrow(network),tree_A$tip.label]
    } else {
    network <- network[tree_B$tip.label,tree_A$tip.label]
    }
  
  # binary Jaccard distances
  if (method=="Jaccard_binary"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
    if (only_A==FALSE) jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=T))
    eco_A <- jaccard_A
    if (only_A==FALSE) eco_B <- jaccard_B
  }
  
  # quantitative Jaccard distances
  if (method=="Jaccard_weighted"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
    if (only_A==FALSE) jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=F))
    eco_A <- jaccard_A
    if (only_A==FALSE) eco_B <- jaccard_B
  }
  
  # Unifrac (generalized UniFrac, with alpha=0.5)
  if (method=="GUniFrac"){
    if (only_A==FALSE) unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=2
    eco_A <- unifrac_A$unifracs[,,index]
    if (only_A==FALSE) eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Unifrac (unweighted UniFrac)
  if (method=="UniFrac_unweighted"){
    if (only_A==FALSE) unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=4
    eco_A <- unifrac_A$unifracs[,,index]
    if (only_A==FALSE) eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Degree
  if (method=="degree"){
    network_binary <- network
    network_binary[network_binary>0] <- 1
    
    eco_A <- as.matrix(dist(colSums(network_binary)))
    if (only_A==FALSE) eco_B <- as.matrix(dist(rowSums(network_binary)))
  }
  
  # PBLM (non binary)
  if ((method=="PBLM")&(only_A==FALSE)){
    model_pblm <- R.utils::withTimeout(pblm(assocs=network, tree1=tree_B, tree2=tree_A, bootstrap=F, nreps=0), timeout = 60*60*24, onTimeout = "silent")
    
    if (!is.null(model_pblm)) {
      results <- c(round(nb_A), round(nb_B), model_pblm$signal.strength$estimate[2], model_pblm$signal.strength$estimate[1], model_pblm$MSE )
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))
    }else{
      results <- c(round(nb_A), round(nb_B), NA, NA, NA, NA, NA, NA)
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))}
    }
  
  # PBLM binary
  if ((method=="PBLM_binary")&(only_A==FALSE)){
    network_binary <- network
    network_binary[network_binary>0] <- 1
    
    model_pblm <- R.utils::withTimeout(pblm(assocs=network_binary, tree1=tree_B, tree2=tree_A, bootstrap=F, nreps=0), timeout = 60*60*24, onTimeout = "silent")
    
    if (!is.null(model_pblm)) {
      results <- c(round(nb_A), round(nb_B), model_pblm$signal.strength$estimate[2], model_pblm$signal.strength$estimate[1], model_pblm$MSE )
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))
    }else{
      results <- c(round(nb_A), round(nb_B), NA, NA, NA, NA, NA, NA)
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))}
  }
  
  # Mantel tests
  if (!method %in% c("PBLM_binary","PBLM")){
    
    # cophenetic distances
    cophe_A <- cophenetic.phylo(tree_A)
    if (only_A==FALSE) cophe_B <- cophenetic.phylo(tree_B)

    results <- c(round(nb_A), round(nb_B), NA, NA, NA, NA, NA, NA)
    names(results) <-  c("nb_A","nb_B","mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")

    if (length(unique(as.vector(cophe_A)))<3) {
      print("The phylogenetic distance matrix of guild A is composed of only 2 different values (because of polytomies?).")
      return(results)}
    if (only_A==FALSE) {if (length(unique(as.vector(cophe_B)))<3) {
      print("The phylogenetic distance matrix of guild B only is composed of only 2 different values (because of polytomies?).")
      return(results)}}
    if (length(unique(as.vector(eco_A)))<3) {
      print("The ecological distance matrix of guild A is composed of only 1 value (identical patterns of interactions across species?).")
      return(results)}
    if (only_A==FALSE) {if (length(unique(as.vector(eco_B)))<3) {
      print("The ecological distance matrix of guild B is composed of only 1 value (identical patterns of interactions across species?).")
      return(results)}}

    
    if (correlation=="Pearson"){
      mantel_A <- RPANDA::mantel_test(as.dist(eco_A) ~ as.dist(cophe_A),  nperm = nperm, correlation="Pearson")
      if (only_A==FALSE) mantel_B <- RPANDA::mantel_test(as.dist(eco_B) ~ as.dist(cophe_B),  nperm = nperm, correlation="Pearson")
    }
    
    if (correlation=="Spearman"){
      mantel_A <- RPANDA::mantel_test(as.dist(eco_A) ~ as.dist(cophe_A),  nperm = nperm, correlation="Spearman")
      if (only_A==FALSE) mantel_B <- RPANDA::mantel_test(as.dist(eco_B) ~ as.dist(cophe_B),  nperm = nperm, correlation="Spearman")
    }
    
    if (correlation=="Kendall"){
      mantel_A <- RPANDA::mantel_test(as.dist(eco_A) ~ as.dist(cophe_A),  nperm = nperm, correlation="Kendall")
      if (only_A==FALSE) mantel_B <- RPANDA::mantel_test(as.dist(eco_B) ~ as.dist(cophe_B),  nperm = nperm, correlation="Kendall")
    }
    
    if (only_A==TRUE) mantel_B <- c(NA, NA, NA)
    
    results <- c(round(nb_A), round(nb_B), mantel_A[1], mantel_A[2], mantel_A[3], mantel_B[1], mantel_B[2], mantel_B[3])
    names(results) <-  c("nb_A","nb_B","mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")
    return(results)
  }
}
