phylosignal_sub_network <-
function(network, tree_A, tree_B=NULL, method = "Jaccard_weighted", 
                                    nperm = 1000, correlation = "Pearson", minimum=10, degree=FALSE, permutation ="shuffle"){
  
  host_tree <- tree_A
  symbiont_tree <- tree_B
  
  if (!inherits(host_tree, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!is.null(tree_B)) {if (!inherits(symbiont_tree, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "Bray-Curtis", "GUniFrac", "UniFrac_unweighted")) {stop("Please provide a \"method\" to compute phylogenetic signals.")}
  
  if (all(is.null(colnames(network)))|all(is.null(rownames(network)))) {stop("Please provide a network with row names and columns names matching the species names.")}
  
  if (method %in% c("GUniFrac", "UniFrac_unweighted")) {if (is.null(tree_B)) stop("Please provide a phylogenetic tree \"tree_B\" for guild B.")}
  
  if (!correlation %in% c("Pearson", "Spearman", "Kendall")) {stop("Please pick a \"correlation\" among Pearson, Spearman, and Kendall.")}
  
  if (nrow(network)<2){stop("Please provide a \"network\" with at least 2 species in clade B.")}
  if (ncol(network)<2){stop("Please provide a \"network\" with at least 2 species in clade A.")}
  
  if (minimum<2){stop("The minimal number of descending species (\"minimum\") must be with at least of 2 (or even larger!).")}
  
  if (!permutation %in% c("shuffle","nbpartners")) {stop("Please provide a type of \"permutation\" among 'shuffle' and 'nbpartners'.")}

  # only keep species having at least one interaction
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  host_tree <- ape::drop.tip(host_tree, tip=host_tree$tip.label[!host_tree$tip.label %in% colnames(network)])
  if (!is.null(symbiont_tree)){
    symbiont_tree <- ape::drop.tip(symbiont_tree, tip=symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(network)])
    network <- network[symbiont_tree$tip.label,host_tree$tip.label]
  }else{
    network <- network[,host_tree$tip.label]
  }
  
  # check a second time (in case of a missing species) 
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  host_tree <- ape::drop.tip(host_tree, tip=host_tree$tip.label[!host_tree$tip.label %in% colnames(network)])
  if (!is.null(symbiont_tree)){
    symbiont_tree <- ape::drop.tip(symbiont_tree, tip=symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(network)])
    network <- network[symbiont_tree$tip.label,host_tree$tip.label]
  }else{
    network <- network[,host_tree$tip.label]
  }
  
  
  nb_sub_clades <- 0
  results_sub_clades <- c()
  for (i in sort(unique(host_tree$edge[,1]))){  # include root and can be non binary
    sub_host_tree <- ape::extract.clade(host_tree, i)
    if (Ntip(sub_host_tree)>=minimum){
      sub_network <- network[,sub_host_tree$tip.label]
      sub_network <- sub_network[which(rowSums(sub_network)>0),,drop=FALSE]
      if (!is.null(symbiont_tree)){
        sub_symbiont_tree <- ape::drop.tip(symbiont_tree, tip= symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(sub_network)])
      }else{sub_symbiont_tree <- NULL}
      
      if (nrow(sub_network)>1){
        nb_sub_clades <- nb_sub_clades+1
        mantel_test <- RPANDA::phylosignal_network(sub_network, sub_host_tree, sub_symbiont_tree, method = method, nperm = nperm, correlation = correlation, permutation = permutation)
        
        if (degree==TRUE){
          mantel_degree <- rep("NA", 5)
          tryCatch({
            mantel_degree <- RPANDA::phylosignal_network(sub_network, sub_host_tree, sub_symbiont_tree, method = "degree", nperm = nperm, correlation = correlation)
          }, error=function(e){cat("clade ",i,": ", conditionMessage(e), "\n")})
          results_sub_clades <- rbind(results_sub_clades, c(i, mantel_test[1:5],NA,NA, mantel_degree[3:5] ))
        }else{
          results_sub_clades <- rbind(results_sub_clades, c(i, mantel_test[1:5],NA,NA))
        }
        
      }
    }
  }
  if (degree==TRUE){
    colnames(results_sub_clades) <- c("node", "nb_A", "nb_B", "mantel_cor", "pvalue_upper", "pvalue_lower", "pvalue_upper_corrected","pvalue_lower_corrected", "degree_mantel_cor", "degree_pvalue_upper", "degree_pvalue_lower") 
  }else{
    colnames(results_sub_clades) <- c("node", "nb_A", "nb_B", "mantel_cor", "pvalue_upper", "pvalue_lower", "pvalue_upper_corrected","pvalue_lower_corrected") 
  }
  results_sub_clades <- data.frame(results_sub_clades, stringsAsFactors = FALSE)
  results_sub_clades$nb_A <- as.integer(as.numeric(results_sub_clades$nb_A))
  results_sub_clades$nb_B <- as.integer(as.numeric(results_sub_clades$nb_B))
  results_sub_clades$mantel_cor <- as.numeric(results_sub_clades$mantel_cor)
  results_sub_clades$pvalue_upper <- as.numeric(results_sub_clades$pvalue_upper)
  results_sub_clades$pvalue_lower <- as.numeric(results_sub_clades$pvalue_lower)
  
  results_sub_clades$pvalue_upper_corrected <- results_sub_clades$pvalue_upper*nb_sub_clades
  results_sub_clades$pvalue_lower_corrected <- results_sub_clades$pvalue_lower*nb_sub_clades
  results_sub_clades$pvalue_upper_corrected[results_sub_clades$pvalue_upper_corrected>1] <- 1
  results_sub_clades$pvalue_lower_corrected[results_sub_clades$pvalue_lower_corrected>1] <- 1
  
  if (degree==TRUE){
    results_sub_clades$degree_mantel_cor <- as.numeric(results_sub_clades$degree_mantel_cor)
    results_sub_clades$degree_pvalue_upper <- as.numeric(results_sub_clades$degree_pvalue_upper)
    results_sub_clades$degree_pvalue_lower <- as.numeric(results_sub_clades$degree_pvalue_lower)
    
    results_sub_clades$degree_pvalue_upper_corrected <- results_sub_clades$degree_pvalue_upper*nb_sub_clades
    results_sub_clades$degree_pvalue_lower_corrected <- results_sub_clades$degree_pvalue_lower*nb_sub_clades
    results_sub_clades$degree_pvalue_upper_corrected[results_sub_clades$degree_pvalue_upper_corrected>1] <- 1
    results_sub_clades$degree_pvalue_lower_corrected[results_sub_clades$degree_pvalue_lower_corrected>1] <- 1
  }
  
  return(results_sub_clades)
}
