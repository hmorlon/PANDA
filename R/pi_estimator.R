pi_estimator <-
function(sequences) {
  sequences <- sequences[,sapply(1:ncol(sequences), function(i) length(which(!sequences[,i] %in% c("-","N","n")))>1),drop=F]
  pi_sites <- sapply(1:ncol(sequences), function(i) pi_per_site(sequences[,i]))
  pi_sites[is.na(pi_sites)] <- 0
  sum(pi_sites)/ncol(sequences)
}
