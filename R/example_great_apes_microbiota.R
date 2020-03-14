example_great_apes_microbiota <-
function(name,path=getwd()){
  data(great_apes_microbiota, envir = environment())
  write.tree(great_apes_microbiota[[1]],file = paste(path,"/host_tree_",name,".tre",sep=""))
  write.dna(great_apes_microbiota[[2]],paste(path,"/alignment_",name,"_OTU0001.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(great_apes_microbiota[[2]]))
  write.dna(great_apes_microbiota[[3]],paste(path,"/alignment_",name,"_OTU0002.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(great_apes_microbiota[[3]]))
  write.dna(great_apes_microbiota[[4]],paste(path,"/alignment_",name,"_OTU0003.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(great_apes_microbiota[[4]]))
}
