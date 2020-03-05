example_great_apes_microbiota <-
function(name,path=getwd()){ 
  data(great_apes_microbiota)
  write.tree(host_tree,file = paste(path,"/host_tree_",name,".tre",sep=""))
  write.dna(alignment_OTU0001,paste(path,"/alignment_",name,"_OTU0001.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(alignment_OTU0001))
  write.dna(alignment_OTU0002,paste(path,"/alignment_",name,"_OTU0002.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(alignment_OTU0002))
  write.dna(alignment_OTU0003,paste(path,"/alignment_",name,"_OTU0003.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(alignment_OTU0003))
}
