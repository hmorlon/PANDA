phylosignal_sub_network <-
function(network, tree_A, tree_B, method = "Jaccard_weighted", 
                                    nperm = 1000, correlation = "Pearson", minimum=10, degree=F){
  
  host_tree <- tree_A
  symbiont_tree <- tree_B
  
  if (!inherits(host_tree, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!inherits(symbiont_tree, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "GUniFrac", "UniFrac_unweighted")) {stop("Please provide a \"method\" to compute phylogenetic signals.")}
  
  if (all(is.null(colnames(network)))|all(is.null(rownames(network)))) {stop("Please provide a network with row names and columns names matching the species names.")}
  
  if (!correlation %in% c("Pearson", "Spearman", "Kendall")) {stop("Please pick a \"correlation\" among Pearson, Spearman, and Kendall.")}
  
  if (nrow(network)<2){stop("Please provide a \"network\" with at least 2 species in clade B.")}
  if (ncol(network)<2){stop("Please provide a \"network\" with at least 2 species in clade A.")}
  
  host_tree <- drop.tip(host_tree, tip=host_tree$tip.label[!host_tree$tip.label %in% colnames(network)])
  symbiont_tree <- drop.tip(symbiont_tree, tip=symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(network)])
  
  network <- network[symbiont_tree$tip.label,host_tree$tip.label]
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  
  set.seed(1)
  
  nb_sub_clades <- 0
  results_sub_clades <- c()
  for (i in sort(unique(host_tree$edge[,1]))){  # include root  and can be non binary
    sub_host_tree <- extract.clade(host_tree, i)
    if (Ntip(sub_host_tree)>=minimum){
      sub_network <- network[,sub_host_tree$tip.label]
      sub_network <- sub_network[which(rowSums(sub_network)>0),,drop=F]
      sub_symbiont_tree <- drop.tip(symbiont_tree, tip= symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(sub_network)])
      
      if (nrow(sub_network)>1){
        nb_sub_clades <- nb_sub_clades+1
        mantel_test <- phylosignal_network(sub_network, sub_host_tree, sub_symbiont_tree, method = method, nperm = nperm, correlation = correlation)
        
        if (degree==TRUE){
        mantel_degree <- rep("NA", 5)
        tryCatch({
          mantel_degree <- phylosignal_network(sub_network, sub_host_tree, sub_symbiont_tree, method = "degree", nperm = nperm, correlation = correlation)
                  }, error=function(e){cat("clade ",i,": ", conditionMessage(e), "\n")})
        results_sub_clades <- rbind(results_sub_clades, c(i, mantel_test[1:5],NA,NA, mantel_degree[3:5] ))
        }else{
          results_sub_clades <- rbind(results_sub_clades, c(i, mantel_test[1:5],NA,NA))
        }
        
      }
    }
  }
  if (degree==TRUE){
    colnames(results_sub_clades) <- c("node", "nb_A", "nb_B", "mantel_cor", "pvalue_high", "pvalue_low", "pvalue_high_corrected","pvalue_low_corrected", "degree_mantel_cor", "degree_pvalue_high", "degree_pvalue_low") 
    }else{
      colnames(results_sub_clades) <- c("node", "nb_A", "nb_B", "mantel_cor", "pvalue_high", "pvalue_low", "pvalue_high_corrected","pvalue_low_corrected") 
    }
  results_sub_clades <- data.frame(results_sub_clades, stringsAsFactors = F)
  results_sub_clades$nb_A <- round(as.numeric(results_sub_clades$nb_A))
  results_sub_clades$nb_B <- round(as.numeric(results_sub_clades$nb_B))
  results_sub_clades$mantel_cor <- as.numeric(results_sub_clades$mantel_cor)
  results_sub_clades$pvalue_high <- as.numeric(results_sub_clades$pvalue_high)
  results_sub_clades$pvalue_low <- as.numeric(results_sub_clades$pvalue_low)
  
  results_sub_clades$pvalue_high_corrected <- results_sub_clades$pvalue_high*nb_sub_clades
  results_sub_clades$pvalue_low_corrected <- results_sub_clades$pvalue_low*nb_sub_clades
  results_sub_clades$pvalue_high_corrected[results_sub_clades$pvalue_high_corrected>1] <- 1
  results_sub_clades$pvalue_low_corrected[results_sub_clades$pvalue_low_corrected>1] <- 1
  
  if (degree==TRUE){
    results_sub_clades$degree_mantel_cor <- as.numeric(results_sub_clades$degree_mantel_cor)
    results_sub_clades$degree_pvalue_high <- as.numeric(results_sub_clades$degree_pvalue_high)
    results_sub_clades$degree_pvalue_low <- as.numeric(results_sub_clades$degree_pvalue_low)
    
  results_sub_clades$degree_pvalue_high_corrected <- results_sub_clades$degree_pvalue_high*nb_sub_clades
  results_sub_clades$degree_pvalue_low_corrected <- results_sub_clades$degree_pvalue_low*nb_sub_clades
  results_sub_clades$degree_pvalue_high_corrected[results_sub_clades$degree_pvalue_high_corrected>1] <- 1
  results_sub_clades$degree_pvalue_low_corrected[results_sub_clades$degree_pvalue_low_corrected>1] <- 1
  }
  
  return(results_sub_clades)
}
