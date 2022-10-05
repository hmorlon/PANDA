theta_estimator <-
function (sequences) {
  sequences <- sequences[,sapply(1:ncol(sequences), function(i) length(which(!sequences[,i] %in% c("-","N","n")))>1),drop=F]
  S <- length(seg_sites_fasta(sequences))
  sum_harmonic <- 0
  for (j in 1:ncol(sequences)){if (length(which(sequences[,j] %in% c("-","N","n")))<nrow(sequences)-1){
    for (i in 1:(length(which(!sequences[,j] %in% c("-","N","n")))-1)){
      sum_harmonic <- sum_harmonic + 1/i}}}
  return(S/sum_harmonic)
}
