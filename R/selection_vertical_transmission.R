selection_vertical_transmission <-
function(name,index){
  table <- read.table(paste("results/results_",name,"_",index,".txt",sep=""),header=T)
  df <- 1
  ratio <- 2*(table$minloglik[table$ksi==0]-min(table$minloglik))[1]
  proba <- 1-pchisq(ratio, df=df)
  results <- rbind(c("-2log(Likelihood_ratio)","df","p-value"),c(ratio, df, proba))
  write.table(results, paste("results/model_selection_vertical_",name,"_",index,".txt",sep=""), row.names=F,col.names=F,quote = F,sep="\t")
}
