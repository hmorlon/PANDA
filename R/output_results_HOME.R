output_results_HOME <-
function(iter,name,name_index,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),nb_tree=10000,empirical,randomize,raref,nb_random=10,figure=FALSE,...){ 

  index <- name_index[iter]
  if (empirical==F){load(file=paste("data/simulation_data_",name,"_",index,".RData",sep=""))}
  if (!file.exists(paste("data/data_model_",name,"_",index,".RData",sep=""))) stop(print("Please start by running the previous steps of HOME (fit_HOME...)"))
  load(paste("data/data_model_",name,"_",index,".RData",sep=""))
  if(!exists("path")) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  setwd(path)
  
  if (N_variant>0&!is.na(N_variant)){
    
    if (empirical==F){ksi <- simul[iter]
    indep <- (simul[iter]=="indep")}
    
    transparent_theme <- theme(panel.grid = element_line("black"),
                               axis.line = element_line("black"),
                               panel.grid.major=element_line("black"),
                               panel.background = element_rect(fill = "transparent",colour = NA),
                               plot.background = element_rect(fill = "transparent",colour = NA))
    transparent_theme_y_only <- theme(panel.grid = element_blank(),
                                      axis.line = element_line("black"),
                                      panel.grid.major.y=element_line("black"),
                                      panel.background = element_rect(fill = "transparent",colour = NA),
                                      plot.background = element_rect(fill = "transparent",colour = NA),
                                      legend.key=element_blank(),legend.background=element_blank())
    transparent_theme_no <- theme(panel.grid = element_blank(),
                                  axis.line = element_line("black"),
                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                  plot.background = element_rect(fill = "transparent",colour = NA),
                                  legend.key=element_blank(),legend.background=element_blank())
    
    # Plot trees (host tree)
    if (empirical==F){
      if (figure==TRUE){
        pdf(paste("figures/host_tree_",name,"_",index,".pdf",sep=""),width=14,height=12) 
        tree <- host_tree
        for(i in 1:n){tree$tip.label[i] <- paste("    ",tree$tip.label[i],sep="")}
        plot(tree,edge.width=4)
        nodes <- as.character(1:(2*n-2))
        for (i in 1:length(nodes)){if (nchar(nodes[length(nodes)])>nchar(nodes[i])){nodes[i] <- paste(paste(rep("0",nchar(nodes[length(nodes)])-nchar(nodes[i])),collapse=""),nodes[i],sep="")}}
        edgelabels(nodes,frame="circle",col="#78281f",bg="#f39c12")
        add.scale.bar()
        dev.off()}
    }
    
    #Plot lik landscape
    results <- read.table(paste("results/results_",name,"_",index,".txt",sep=""),header=T)
    results <- results[order(results$ksi),]
    results <- results[which(is.finite(results$minloglik)),]
    est_ksi <- results$ksi[which.min(results$minloglik)]
    
    if (nrow(results)>1){
      if (figure==TRUE){
        if (empirical==F){
          p <- ggplot(results, aes(x=ksi))+
            geom_line(alpha=0.8, aes(y=minloglik), colour="#d35400",size=1.5)+
            geom_point(aes(x=ksi,y=minloglik), color="#d35400")+
            geom_hline(yintercept = simulated_likelihood,color="#154360",linetype="dashed", size=1)+  #### SIMULATIONS
            labs(x = "Number of switches", y="-log(Likelihood)")+transparent_theme+scale_x_continuous(breaks=seq(0,max(max(lambda),10),5))
        }else{p <- ggplot(results, aes(x=ksi))+
          geom_line(alpha=0.8, aes(y=minloglik), colour="#d35400",size=1.5)+
          geom_point(aes(x=ksi,y=minloglik), color="#d35400")+
          labs(x = "Number of switches", y="-log(Likelihood)")+transparent_theme+scale_x_continuous(breaks=seq(0,max(max(lambda),10),5))
        }
        pdf(paste("figures/profil_model_switches_",name,"_",index,".pdf",sep=""),width=5,height=4) 
        print(p)
        dev.off()}
      
      if (figure==TRUE){
        p <- ggplot(results, aes(x=ksi))+ 
          geom_line(alpha=0.8, aes(y=mu), colour="#d35400",size=1.5)+
          geom_point(aes(x=ksi,y=mu), color="#d35400")+scale_x_continuous(breaks=seq(0,max(max(lambda),10),5))+
          labs(x = "Number of switches", y="Substitution rate")+transparent_theme
        pdf(paste("figures/profil_model_switches_mu_",name,"_",index,".pdf",sep=""),width=5,height=4) 
        print(p)
        dev.off()}
      
      # Plot loglikelihood distibution
      if (est_ksi!=0){
        results_ll <- read.table(paste("results/optim_ll_",name,"_",index,"_",est_ksi,".txt",sep=""),header=F)
        if (figure==TRUE){
          if (empirical==F){p <- ggplot(results_ll, aes(x=""))+
            geom_violin(alpha=0.5, aes(y=V1),fill="#d35400", colour="#d35400",size=1.5,trim=T)+
            geom_hline(yintercept = simulated_likelihood,color="#154360",linetype="dashed", size=1)+
            labs(x = "Simulated trees (nb switches)", y="-log(Likelihood)")+transparent_theme#_y_only#+scale_x_continuous(breaks=c())
          }else{p <- ggplot(results_ll, aes(x=""))+
            geom_violin(alpha=0.5, aes(y=V1),fill="#d35400", colour="#d35400",size=1.5,trim=T)+
            labs(x = "Simulated trees (nb switches)", y="-log(Likelihood)")+transparent_theme#_y_only#+scale_x_continuous(breaks=c())
          }
          print("test1")
          pdf(paste("figures/profil_ll_",name,"_",index,"_",est_ksi,".pdf",sep=""),width=5,height=4)
          print(p)
          dev.off()
          print("test2")
          }
      }
      
      #Plot vertical transmission
      results_vertical <- read.table(paste("results/model_selection_vertical_",name,"_",index,".txt",sep=""),header=T)
      if (is.finite(results_vertical[1,1])){ 
        x <- seq(0.1,max((results_vertical[1,1]+5),10),0.01)
        chi_square <- data.frame(cbind(x,dchisq(x, df=1)))
        colnames(chi_square) <- c("x","proba")
        p <- ggplot(chi_square, aes(x=x))+ 
          geom_line(alpha=0.8, aes(y=proba), colour="#333333",size=1)+
          geom_vline(xintercept=results_vertical[1,1],linetype="dashed",size=1.5,colour="#d35400")+
          geom_vline(xintercept=qchisq(0.95,df=1),size=0.8,colour="#154360")+
          geom_vline(xintercept=qchisq(0.99,df=1),size=0.8,colour="#5499c7")+
          labs(x = "-2log(vertical null model / model switches)", y="Probability density")+transparent_theme_no
        if (figure==TRUE){
          pdf(paste("figures/pure_vertical_transmission_",name,"_",index,".pdf",sep=""),width=5,height=4) 
          print(p)
          dev.off()}
      }
      
      #Plot independent evolution
      est_ksi <- results$ksi[which.min(results$minloglik)]
      est_mu <- results$mu[which.min(results$minloglik)]
      if(randomize==T){
        results_randomize <- read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T)
        if (nrow(results_randomize)!=nb_random){print(paste0("Problem with number of replicates in randomizations: ", index))}
        ksi_randomize <- length(which(results_randomize$ksi<=est_ksi))/length(results_randomize$ksi)
        mu_randomize <- length(which(round(results_randomize$mu,digits=10)<=round(est_mu,digits=10)))/length(results_randomize$mu)
        write.table(rbind(c("Distribution_test","p-value"),c("Empirical_ranking_(ksi_distribution)",round(c(ksi_randomize),4)),c("Empirical_ranking_(mu_distribution)",round(c(mu_randomize),4))), paste("results/model_selection_independent_stats_",name,"_",index,".txt",sep=""), row.names=F,col.names=F,quote = F,sep="\t")
        results_ll <- results[which.min(results$minloglik),]
        results_ll$replicate <- 0
        results_randomize <- rbind(results_randomize,results_ll)
        if (figure==TRUE){
          p <- ggplot(results_randomize, aes(x=ksi))+
            geom_point(aes(x=ksi,y=mu), shape=21, fill="#154360", color="#333333",size=3.5)+transparent_theme+
            labs(x = "Estimated number of switches", y="Estimated substitution rate")+scale_x_continuous(breaks=seq(0,max(max(lambda),10),5))+
            geom_point(aes(x=est_ksi,y=est_mu), shape=21, fill="#d35400",color="#333333", size=3)
          pdf(paste("figures/independent_evolution_randomize_",name,"_",index,".pdf",sep=""),width=5,height=4)
          print(p)
          dev.off()}
      }
      
      #Rarefactions
      if (raref==T){
        table <- seq(1,nb_tree,nb_tree/20)
        table <- cbind(table,rep(results$minloglik[which(results$ksi==0)],20))
        for (ksi in lambda){
          table <- cbind(table,as.numeric(read.table(paste(path,"/rarefactions/loglik_rarefaction_",name,"_",index,"_",ksi,".txt",sep=""))[,2]))
        }
        table <- rbind(table,c(nb_tree,results$minloglik))
        table <- data.frame(table)
        colnames(table) <- c("table","0",lambda)
        table <- melt(table, id.vars=c("table"))
        colnames(table) <- c("table","ksi","value")
        p <- ggplot(table,aes(x=table,col=ksi))+geom_line(aes(y=value),size=1.25) +labs(x = "Number of simulated trees", y="-log(Likelihood)")+transparent_theme
        pdf(paste("figures/rarefactions_",name,"_",index,".pdf",sep=""),width=8,height=5) 
        print(p)
        dev.off()
      }
      
      ###### File resume ##
      if(randomize==T){
        write.table(cbind("name",index),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)
        write.table(cbind("n-host",n),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("N_variant",N_variant),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_mu",est_mu),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_ksi",est_ksi),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind(selected_model,Q[1,3]),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("strict-vertical",as.numeric(results_vertical[1,3])),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("indep-ksi",ksi_randomize),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("indep-mu",mu_randomize),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
      }else{
        write.table(cbind("name",index),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)
        write.table(cbind("n-host",n),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("N_variant",N_variant),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_mu",est_mu),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_ksi",est_ksi),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind(selected_model,Q[1,3]),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("strict-vertical",as.numeric(results_vertical[1,3])),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
      }
      
      if (figure==TRUE){
        HTMLStart(outdir=paste(path,"/figures/",sep=""), file=paste("report_",name,"_",index,sep=""), extension="html", echo=FALSE, HTMLframe=T)
        HTML.title(paste("Results",name,index,sep=" "), HR=1)
        HTML.title("Description of the data", HR=2)
        if (empirical==F){
          HTML(rbind(c("Number of symbiont-host association:",n),
                     c("Simulated substitution rate:",simulated_mu), 
                     c("Number of simulated switches:",ksi),   ## ONLY SIMULATIONS
                     c("Seed for simulations:",seed),   ## ONLY SIMULATIONS
                     c("Sequences length:",N),
                     c("Probability variant sites:",proportion_variant),
                     c("Number of invariant sites:",N_invariant),
                     c("Number of strictly variant sites:",N_variant)))
        }else {
          HTML(rbind(c("Number of symbiont-host association:",n),
                     c("Sequences length:",N),
                     c("Number of invariant sites:",N_invariant),
                     c("Number of strictly variant sites:",N_variant)))
        }
        
        if (nrow(results)>1){
          HTML.title("Summary of the most likely scenario:", HR=2)
          table <- read.table(paste("results/results_",name,"_",index,".txt",sep=""),header=T)
          if ((empirical==F)&(randomize==F)){
            HTML(rbind(c("Substitution model:",selected_model),
                       c("Transition/Transversion ratio:",round(Q[1,3],4)),
                       c("Estimated substitution rate:",round(est_mu,4)),
                       c("Estimated number of switches:",round(est_ksi,4)),
                       c("Associated likelihood:",round(table$minloglik[which.min(table$minloglik)],4)),
                       c("p-value: Strict vertical transmission:",round(as.numeric(results_vertical[1,3]),4))))}
          if ((empirical==F)&(randomize==T)){
            HTML(rbind(c("Substitution model:",selected_model),
                       c("Transition/Transversion ratio:",round(Q[1,3],4)),
                       c("Estimated substitution rate:",round(est_mu,4)),
                       c("Estimated number of switches:",round(est_ksi,4)),
                       c("Associated likelihood:",round(table$minloglik[which.min(table$minloglik)],4)),
                       c("p-value: Strict vertical transmission:",round(as.numeric(results_vertical[1,3]),4)),
                       c("p-value: Independent evolutions (nb switches):",round(ksi_randomize,4)),
                       c("p-value: Independent evolutions (subs. rate)",round(mu_randomize,4))))}
          if ((empirical==F)&(randomize==F)){
            HTML(rbind(c("Simulated substitution rate:",round(simulated_mu,4)),
                       c("Original likelihood:",round(simulated_likelihood,4)),
                       c("Substitution model:",selected_model),
                       c("Transition/Transversion ratio:",round(Q[1,3],4)),
                       c("Estimated substitution rate:",round(est_mu,4)),
                       c("Estimated number of switches:",round(est_ksi,4)),
                       c("Associated likelihood:",round(table$minloglik[which.min(table$minloglik)],4)),
                       c("p-value: Strict vertical transmission:",round(as.numeric(results_vertical[1,3]),4))))}
          if ((empirical==T)&(randomize==T)){
            HTML(rbind(c("Substitution model:",selected_model),
                       c("Transition/Transversion ratio:",round(Q[1,3],4)),
                       c("Estimated substitution rate:",round(est_mu,4)),
                       c("Estimated number of switches:",round(est_ksi,4)),
                       c("Associated likelihood:",round(table$minloglik[which.min(table$minloglik)],4)),
                       c("p-value: Strict vertical transmission:",round(as.numeric(results_vertical[1,3]),4)),
                       c("p-value: Independent evolutions (nb switches):",round(ksi_randomize,4)),
                       c("p-value: Independent evolutions (subs. rate)",round(mu_randomize,4))))}
          if ((empirical==T)&(randomize==F)){
            HTML(rbind(c("Substitution model:",selected_model),
                       c("Transition/Transversion ratio:",round(Q[1,3],4)),
                       c("Estimated substitution rate:",round(est_mu,4)),
                       c("Estimated number of switches:",round(est_ksi,4)),
                       c("Associated likelihood:",round(table$minloglik[which.min(table$minloglik)],4)),
                       c("p-value: Strict vertical transmission:",round(as.numeric(results_vertical[1,3]),4))))}
          
          table <- read.table(paste("results/model_selection_vertical_",name,"_",index,".txt",sep=""),header=T)
          strict_vertical <- as.numeric(table[3])
          if(strict_vertical>=0.05){HTML(rbind("CONCLUSION:","STRICT VERTICAL TRANSMISSION"," "))}
          if (randomize==T){if(ksi_randomize>0.05){HTML(rbind("CONCLUSION:","INDEPENDENT EVOLUTIONS (NUMBER OF SWITCHES)"," "))}
            if((ksi_randomize>=0.05)&(ksi_randomize<0.15)){HTML(rbind("INCREASE THE NUMBER OF RANDOMIZATIONS (NUMBER OF SWITCHES)"))}
            if(mu_randomize>0.05){HTML(rbind("CONCLUSION:","INDEPENDENT EVOLUTIONS (SUBSTITUTION RATE)"," "))}
            if((mu_randomize>=0.05)&(mu_randomize<0.15)){HTML(rbind("INCREASE THE NUMBER OF RANDOMIZATIONS (SUBSTITUTION RATE)"))}}
          
          HTML.title("Host-switches inference:", HR=2)
          table <- read.table(paste("results/results_",name,"_",index,".txt",sep=""),header=T)
          HTML(rbind(c("Ksi","Mu","-log(Likelihood)"),c(round(as.numeric(table[which.min(table$minloglik),1:3]),4))))
          HTML.title("Profil of minus log likelihood as a function of the number of switches.",HR=4)
          HTMLInsertGraph(paste("profil_model_switches_",name,"_",index,".pdf",sep=""),WidthHTML=400, Caption="", GraphBorder=0,Align="center")
          HTML.title("Estimated substitution rate as a function of the number of switches.",HR=4)
          HTMLInsertGraph(paste("profil_model_switches_mu_",name,"_",index,".pdf",sep=""),WidthHTML=400, Caption="", GraphBorder=0,Align="center")
          HTML.title("Strict vertical transmission model (rejected if p-value < 0.05):", HR=2)
          table <- read.table(paste("results/model_selection_vertical_",name,"_",index,".txt",sep=""),header=T)
          HTML(rbind(c("Ratio log(Likelihood)","df","p-value"),round(as.numeric(table[1,1:3]),4)))
          HTML.title("Results of the likelihood ratio test. The grey curve correponds to the Chi2 distribution with df=1. The dark blue line (resp. light) stands for the 0.05 (resp. 0.01) p-value threshold and the dashed orange line is the observed LRT ratio.",HR=4)
          HTMLInsertGraph(paste("pure_vertical_transmission_",name,"_",index,".pdf",sep=""),WidthHTML=400, Caption="", GraphBorder=0,Align="center")
          
          if(randomize==T){
            HTML.title("Host-symbiont independent evolutions (rejected if p-values < 0.05):", HR=2)
            HTML.title("Representation of the estimated numbers of switches and the estimated substitution rates for the randomized alignments (in blue) and the empirical alignment (in orange). Independent evolutions can be rejected if the orange dot stands alone in the bottom left corner.",HR=4)
            HTMLInsertGraph(paste("independent_evolution_randomize_",name,"_",index,".pdf",sep=""),WidthHTML=400, Caption="", GraphBorder=0,Align="center")
            HTML(rbind(c("Distribution test","p-value"),c("Empirical ranking (ksi distribution)",round(c(ksi_randomize),4)),c("Empirical ranking (mu distribution)",round(c(mu_randomize),4))),header=T)
          }
          
          HTML.title("Estimated substitution model", HR=2)
          HTML.title("Estimated rate matrix:", HR=3)
          HTML(cbind(c(" ","A","C","G","T"),rbind(c("A","C","G","T"),round(Q,2))))
          HTML.title("Nucleotide frequencies:", HR=3)
          HTML(rbind(c("A","C","G","T"),round(propinv,2)))
          
          if (empirical==F){if (indep==F){if (as.numeric(ksi)>0){
            HTML.title("Estimated switches", HR=2)
            HTML.title("Simulated switches:",HR=3)
            table <- read.table(paste("simulations/simulated_switches_",name,"_",index,".txt",sep=""),header=F)
            table <- cbind(c("Branch origin","Branch arrival","Absolute position"),table)
            colnames(table) <- c("Simulated switch(es):",as.factor(c(1:(ncol(table)-1))))
            HTML(table,row.names=FALSE)}}
            
            if (file.exists(paste(path,"/figures/host_tree_switches_",name,"_",index,".pdf",sep=""))){HTMLInsertGraph(paste("host_tree_switches_",name,"_",index,".pdf",sep=""),WidthHTML=600, Caption="", GraphBorder=0,Align="center")
            }else{HTMLInsertGraph(paste("host_tree_",name,"_",index,".pdf",sep=""),WidthHTML=600, Caption="", GraphBorder=0,Align="center")}
            
            if (raref==T){
              HTML.title("Rarefactions:", HR=2)
              HTML.title(c("-log(Likelihood) as a function of the number of simulated trees"),align="center",HR=4)
              HTMLInsertGraph(paste("rarefactions_",name,"_",index,".pdf",sep=""),WidthHTML=600, Caption="", GraphBorder=0,Align="center")
            } 
            HTML.title("More detailed table:", HR=3)
            HTML.title("Host switches inference:")
            results <- read.table(paste("results/results_",name,"_",index,".txt",sep=""),header=T)
            HTML(round(results[order(results$ksi),],4),row.names=FALSE)
            
            if(randomize==T){
              HTML.title("Independent evolution results:", HR=3)
              results <- read.table(paste("results/model_selection_independent_",name,"_",index,".txt",sep=""),header=T)
              HTML(round(results,4),row.names=FALSE)
            }
          }}
        
        HTMLStop()
        file.remove(file=paste(path,"/figures/","report_",name,"_",index,".html",sep=""))
        file.remove(file=paste(path,"/figures/","report_",name,"_",index,"_menu.html",sep=""))
        
      }}}else{ ### NO VARTIATION WITHIN ALIGNMENT
        write.table(cbind("name",index),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)
        write.table(cbind("n-host",n),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("N_variant",N_variant),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_mu",NA),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_ksi",NA),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind(NA,NA),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("strict-vertical",NA),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("indep-ksi",NA),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("indep-mu",NA),paste("figures/results_randomize_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        
        write.table(cbind("name",index),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=F)
        write.table(cbind("n-host",n),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("N_variant",N_variant),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_mu",NA),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("est_ksi",NA),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind(NA,NA),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
        write.table(cbind("strict-vertical",NA),paste("figures/results_",name,"_",index,".csv",sep=""),col.names=F,row.names=F,sep=";",quote=F,append=TRUE)
      }
}
