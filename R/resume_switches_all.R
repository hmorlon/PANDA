resume_switches_all <-
function(index,name,nb_tree,host_tree=NULL){
  dir.create(file.path("infer_PHS/"), showWarnings = FALSE)
  results <- read.table(paste("results/results_",name,"_",index,".txt",sep=""),header=T)
  results <- results[order(results$ksi),]
  results <- results[which(is.finite(results$minloglik)),]
  est_ksi <- results$ksi[which.min(results$minloglik)]
  if (is.null(host_tree)) {
    host_tree <- read.tree(paste0("host_tree_",name,".tre"))
    host_tree <- ladderize(host_tree)
    host_tree$edge.length <- host_tree$edge.length/sum(host_tree$edge.length)
  }
  resume_switches(name,index,est_ksi,nb_tree,host_tree)
}
