fit_HOME <-
function(index,name,nb_tree=10000,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),nb_cores=1,raref=FALSE,...){ 
  print(paste("Index: ",index,sep=""))
  if (!file.exists(paste("data/data_model_",name,"_",index,".RData",sep=""))) stop(print("Please start by running the previous steps of HOME (prepare_data_HOME...)"))
  
  if (file.exists(paste("data/alignment_variant_",name,"_",index,".fas",sep=""))){
    load(paste("data/data_model_",name,"_",index,".RData",sep=""))
    if (N_variant>0){
      
      print("test")
      Sys.sleep(10)
      
      #output <- optimize(f=inference_vertical_transmission,lower=0.0001,upper=50,tol=0.05,name=name,index=index,sequences=variant_sequences) 
      
      print("test_end")
      
      results <- cbind(0,output$minimum,output$objective)
      colnames(results) <- c("ksi","mu","minloglik")
      write.table(results, paste("results/results_",name,"_",index,".txt",sep=""), append=FALSE, col.names=T, row.names=F,quote = F,sep="\t")
      
      inference_ksi <- function(ksi){
        print(paste("Number of switch(es): ",ksi,sep=""))
        load(file=paste("simulated_trees/simulated_trees_",name,"_",ksi,".RData",sep=""))
        
        if (n!=Ntip(list_tree[[1]])) {missing_symbiont <- setdiff(list_tree[[1]]$tip.label,rownames(variant_sequences))
        for(i in 1:nb_tree){for (missing in missing_symbiont){list_tree[[i]] <- drop.tip(list_tree[[i]], tip=missing)}}}
        
        #output <- optimize(f=inference_switches, lower=0.0001, upper=50,tol=0.05,ksi=ksi, randomize=F, name=name,index=index,sequences=variant_sequences,nb_tree=nb_tree,list_tree=list_tree,eig_val=eig_val, eig_vect=eig_vect, ivp=ivp, propinv=propinv)
        
        mu <- output$minimum
        minloglik <- output$objective
        
        if (is.finite(minloglik)){
          results <- cbind(ksi,mu,minloglik)
          colnames(results) <- c("ksi","mu","minloglik")
          write.table(results, paste("results/results_",name,"_",index,".txt",sep=""),append=TRUE, col.names=FALSE, row.names=F,quote = F,sep="\t")}
        
        if (raref==T){ #####   rarefactions
          list <- seq(1,nb_tree,nb_tree/20)
          list_loglik <- vector("numeric",length(list))
          for (i in 1:length(list)) {
            
            output <- optimize(f=inference_switches, lower=0.0001, upper=50,tol=0.05,ksi=ksi,randomize=T, name=name,index=index,sequences=variant_sequences,nb_tree=list[i],list_tree=list_tree[c(1:list[i])],eig_val=eig_val, eig_vect=eig_vect, ivp=ivp, propinv=propinv)
            
            list_loglik[i] <- output$objective}
          list_loglik<- cbind(list,as.vector(list_loglik))
          write.table(list_loglik, paste("rarefactions/loglik_rarefaction_",name,"_",index,"_",ksi,".txt",sep=""),col.names=F,row.names = F)
        } 
      }
      
      output <- mclapply(lambda,inference_ksi,mc.cores=nb_cores)
      
      ####  Model selection (vertical transmission) to save time 
      output <- selection_vertical_transmission(name,index)
    }}}
