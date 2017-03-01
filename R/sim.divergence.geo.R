sim.divergence.geo<-function(phylo,pars, Nsegments=2500, plot=FALSE, geo.object){
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
  nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
  root <- length(phylo$tip.label) + 1
  heights<-phytools::nodeHeights(phylo)
  totlen<-max(heights)
  len<-root-1
  nodeDist<-c(as.numeric(sort(max(ape::branching.times(phylo))-ape::branching.times(phylo))),totlen)	
  nodeDiff<-diff(nodeDist)
  if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
  old.labels<-as.numeric(names(sort(ape::branching.times(phylo),decreasing=TRUE)))
  old.edge<-phylo$edge
  if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
  	checkmat<-cbind(old.labels,seq(root,len+phylo$Nnode))
  	old.edge<-phylo$edge
  	for(j in 1:phylo$Nnode){phylo$edge[which(old.edge==checkmat[j,1])]<-checkmat[j,2]}
  	}
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:phylo$Nnode){
    other<-phylo$edge[phylo$edge[,1]==i+len, 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+len
      if(b>len){
        counter_three_letters <- counter_three_letters + 1
        int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
        int[3]<-b
      } else {
        int[2]<-phylo$tip.label[[b]]
        int[3]<-0 ##NOTE :I am considering tips to be "0" here and use this below
      }
      mat<-rbind(mat,int)
    }
  }		
  nat<-list()
  branchesPresent <- rep(NA, length(nodeDiff)) 
  for(i in 1:length(nodeDiff)){
    if(i==1){
        nat[[i]]<-mat[mat[,1]==(len+i),2]
      } else {
        IN<-vector()
        P<-mat[as.numeric(mat[,1])<=(len+i),c(2,3)]
        IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(len+i),1])
        nat[[i]]<-IN
      }
    branchesPresent[i] = length(nat[[i]])
  }	
    
  masterbranch<-list()
  segmentSize <- rep(NA, phylo$Nnode)
  mappings <- list()
  segsize = sum(nodeDiff)/Nsegments
  newDist<-geo.object$times
  newDiff<-geo.object$spans
  geography.object<-geo.object$geography.object

  for(i in 1:phylo$Nnode){ ##for each node interval
    
    nati<-nat[[i]]
    mappings[[i]] = rep(NA, branchesPresent[i])
    
    if(i==1){
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
    }else{      
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
      for (j in 1:branchesPresent[i]){
        if(length(which(nati[j] == pastNati)) > 0){ 
          mappings[[i]][j] = which(nati[j] == pastNati)
        }else{
          mappings[[i]][j] = which(mat[mat[,3]==mat[mat[,2]==nati[j],1],2] == pastNati)
        }
      }     
    }   
    pastNati <- nati
  }
  
  ## Looping over the parameters
  
  out <- matrix( nrow = nrow(pars), ncol = len)

  for (p in 1:nrow(pars)){
    sig2 = pars[p,1]
    max = pars[p,2]
    alpha = pars[p,3]
    root.value = pars[p,4]
    psi=pars[p,5]
    theta=pars[p,6]
    
    timecount=1      
    for(i in 1:phylo$Nnode){   
        traitMat <- matrix(nrow = branchesPresent[i], ncol = segmentSize[i]+1)
        
        if (i == 1){
          traitMat[,1] = root.value
        }else{
          traitMat[,1] = masterbranch[[i-1]][mappings[[i]],(segmentSize[i-1]+1)] #added +1 here since segmentSize[i-1] is penultimate column, not the last one
        }
        
        tempInd<- 1:branchesPresent[i] # hack to have fast selection of not k, seemed to be faster than a call to which()
        
        for(k in 1:segmentSize[i]){        
          	for(j in 1:branchesPresent[i]){
	          	elms<-which(geography.object[[timecount]][match(unlist(nat[[i]])[j],rownames(geography.object[[timecount]])),]==1)##these are elements that are equal to 1 (sympatric lineages);if length(elms)=0 the repulsion component just falls out
	            diffMeOthers <- traitMat[j,k] - traitMat[elms[which(elms!=j)],k]
	            sign1<-sign(diffMeOthers) ##Alternative: sign1<-ifelse(diffMeOthers>0,1,-1), but diff seems to be faster
	            temp1 = max *segsize #old version:temp1 = 1/(length(elms)) * max *segsize 
      		    temp2 = sqrt(sig2*segsize)
	            if(any(sign1==0)){
	            	no.eq<-which(traitMat[,k] ==traitMat[j,k])
					if(j==min(no.eq)){
						sign1<-sign(rnorm(1))
	            		traitMat[j,k+1]<-traitMat[j,k] +psi*(theta-traitMat[j,k])*segsize + temp1 * sum(sign1*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2)
					} else{
						minn<-min(no.eq)
						sign2<-sign(traitMat[minn,k+1]-traitMat[minn,k])
	            		traitMat[j,k+1]<-traitMat[j,k] +psi*(theta-traitMat[j,k])*segsize + temp1 * sum(-sign2*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2)
					}			
	            } else{ # this is necessary because without this, if lineages have identical trait values, they do not have repulsion (instead repulsion should be strongest when lineages are equal). Also, ifelse approach provides one sign value for all identical matches, which could be fine but this randomly chooses sign for each 0.
	            	traitMat[j,k+1]<-traitMat[j,k] + psi*(theta-traitMat[j,k])*segsize + temp1 * sum(sign1*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2)
              	}
				if((k!=segmentSize[i])&&(j==branchesPresent[i])){timecount= ifelse(round(((k*segsize)+nodeDist[i]),8)>=round(newDist[timecount+1],8),timecount+1,timecount)}
					###loop for last segment size (to preserve exact branch lengths)
					if(k==segmentSize[i] && nodeDiff[i]%%segsize!=0){
						segsizeT= nodeDiff[i]%%segsize
	            		diffMeOthers <- traitMat[j,k] - traitMat[elms[which(elms!=j)],k]
						temp1B =  max *segsizeT #old version 1/length(elms) * max *segsizeT 
	        			temp2B = sqrt(sig2*segsizeT)
	            		sign1<-sign(diffMeOthers) # Alternative: sign1<-ifelse(diffMeOthers>0,1,-1), but diff seems to be faster
	            		if(any(sign1==0)){
	            			no.eq<-which(traitMat[,k] ==traitMat[j,k])
							if(j==min(no.eq)){
								sign1<-sign(rnorm(1))
	            				traitMat[j,k+1]<-traitMat[j,k] + psi*(theta-traitMat[j,k])*segsizeT +temp1B * sum(sign1*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2B)
							} else{
								minn<-min(no.eq)
								sign2<-sign(traitMat[minn,k+1]-traitMat[minn,k])
	            				traitMat[j,k+1]<-traitMat[j,k] + psi*(theta-traitMat[j,k])*segsizeT +temp1B * sum(-sign2*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2B)
							}			
	            		} else{ # this is necessary because without this, if lineages have identical trait values, they do not have repulsion (instead repulsion should be strongest when lineages are equal). Also, ifelse approach provides one sign value for all identical matches, which could be fine but this randomly chooses sign for each 0.
	            		traitMat[j,k+1]<-traitMat[j,k] + psi*(theta-traitMat[j,k])*segsizeT + temp1B * sum(sign1*exp(-alpha*(diffMeOthers)^2)) +rnorm(1,0,temp2B)
						}
					 	if(j==branchesPresent[i]){timecount= ifelse(round((((k-1)*segsize)+nodeDist[i]+segsizeT),8)>=round(newDist[timecount+1],8),timecount+1,timecount)}
						}           
	        	}
    	}
    	
        masterbranch[[i]] = traitMat
      }
    out[p,] = as.vector(masterbranch[[phylo$Nnode]][, tail(segmentSize,n=1)+1])      
    if(p==1){colnames(out)<-unlist(nat[[i]])}			
  }

  if(plot==TRUE){
    print("plotting last simulated dataset")
    M=seq(0,sum(nodeDiff),length=sum(segmentSize))
    O=list()
    for(i in 1:length(nodeDiff)){
    O[[i]]<-data.frame(seq(nodeDist[[i]],nodeDist[[i+1]],length=(segmentSize[[i]]+1)),as.data.frame(t(masterbranch[[i]])))
    }
    t.plot<-plot(M,1:length(M),col="white", ylim=c(range(sapply(masterbranch,range))), xlab="time", ylab="Value")
    for(i in 1:length(nodeDiff)){
    	for(j in 1:(ncol(O[[i]])-1)){
    		lines(O[[i]][,1],O[[i]][,j+1])
    	}	
    }
  }

  return(out)
}

