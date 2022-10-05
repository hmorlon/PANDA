seg_sites_fasta <-
function (x) {
  S <- which(sapply(1:ncol(x), function(i) length(which(!unique(x[,i]) %in% c("-","N","n"))))>1)
  return(S)
}
