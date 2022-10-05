pi_per_site <-
function(sequences){
  alleles <- table(sequences)
  alleles <- alleles[!names(alleles) %in% c("-","N","n")]
  n_alleles <- sum(alleles)
  return(1-sum(alleles*(alleles-1))/(n_alleles*(n_alleles-1)))
}
