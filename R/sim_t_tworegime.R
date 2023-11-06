sim_t_tworegime<-function(regime.map,pars,root.value,Nsegments=2500,model=c("MC","DDexp","DDlin","EB"),verbose=TRUE,rnd=6){

phylo<-regime.map
#return error if non-ultrametric tree
if(!is.ultrametric(phylo)){stop("phylo object must be ultrametric")}
if(is.na(match(model,c("MC","DDexp","DDlin","EB")))){stop("model not specified correctly, must be 'MC','DDexp', 'DDlin',or 'EB'")}

if(is.null(regime.map)){
	stop("provide a regime.map (see ?sim_t_tworegime())")
} else {
	class.object<-try(CreateClassObject(regime.map))
	if(inherits(class.object, "try-error")){
		class.object<-CreateClassObject(regime.map,rnd=6)
		}
				
	SMatrix<-.CreateSMatrix(class.object)
	if(verbose){
	cat(paste0("regime 1 is ",SMatrix$S1,"; regime 2 is ",SMatrix$S2))
	}
	if(length(pars)!=3){stop("pars must be a vector with a value for sig2 and either S1 and S2 for MC model, r1 and r2 for DDexp model, b1 and b2 for DDlin model, or r1 and r2 for EB model")}
	S1=pars[2]
	S2=pars[3]
	sig2=pars[1]
	
	smat = SMatrix$S.matrix

}


paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS

root.value<-root.value

##define a few useful variables
##nodeDist, the distance between each node and the root
##nodeDiff, the distance between each node and the previous node (times for integration)

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
segsize = sum(nodeDiff)/2500
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



if(model=="EB"){
	if((sig2/(sig2*exp(S1*totlen)))>=1e3){warning("r1 parameter leads to sig2 at tips that is 1000x smaller than root value, consider changing")}
	if((sig2/(sig2*exp(S2*totlen)))>=1e3){warning("r2 parameter leads to sig2 at tips that is 1000x smaller than root value, consider changing")}

	if(!all(colnames(smat[[length(smat)]])==nat[[length(nat)]])){stop("order of S.matrix species is incorrect; perhaps this is a sorted matrix?")}	
    newDist<-SMatrix$times
    newDiff<-SMatrix$spans
	S1.first.time=newDist[min(which(lapply(smat,function(x) sum(x[1,])>0) ==TRUE))]
	S2.first.time=newDist[min(which(lapply(smat,function(x) sum(x[2,])>0) ==TRUE))]
	
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
				current.time = newDist[timecount]+(k*segsize)
	          	
	          	for(j in 1:branchesPresent[i]){
	            
       		 		tempS1 =  ifelse(smat[[timecount]][1,j]>0,exp(S1*(current.time-S1.first.time)),0)
       		 		tempS2 =  ifelse(smat[[timecount]][2,j]>0,exp(S2*(current.time-S2.first.time)),0)
            		traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(sig2*((tempS1+tempS2)/sum(smat[[timecount]][,j]))*segsize))
				
				
				if((k!=segmentSize[i])&&(j==branchesPresent[i])){timecount= ifelse(round(((k*segsize)+nodeDist[i]),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
				###loop for last segment size (to preserve exact branch lengths)
				if(k==segmentSize[i] && nodeDiff[i]%%segsize!=0){
					segsizeT= nodeDiff[i]%%segsize
					current.time = newDist[timecount]+((k-1)*segsize)+ segsizeT
       		 		
       		 		tempS1 =  ifelse(smat[[timecount]][1,j]>0,exp(S1*(current.time-S1.first.time)),0)
       		 		tempS2 =  ifelse(smat[[timecount]][2,j]>0,exp(S2*(current.time-S2.first.time)),0)
	            	traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(sig2*((tempS1+tempS2)/sum(smat[[timecount]][,j]))*segsizeT))
			
				 	if(j==branchesPresent[i]){timecount= ifelse(round((((k-1)*segsize)+nodeDist[i]+segsizeT),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
					}
	        	}	        	
			}
        masterbranch[[i]] = traitMat
      }
    out = as.vector(masterbranch[[phylo$Nnode]][, tail(segmentSize,n=1)+1])
    names(out)<-unlist(nat[[i]])
}


if(model=="MC"){

	if(!all(colnames(smat[[length(smat)]])==nat[[length(nat)]])){stop("order of S.matrix species is incorrect; perhaps this is a sorted matrix?")}	
  	newDist<-SMatrix$times
  	newDiff<-SMatrix$spans
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
           
           elmsS1<- if(smat[[timecount]][1,j]==1){which(smat[[timecount]][1,]==1)}else{logical(0)}#these elements are sympatric and in S1
           elmsS2<- if(smat[[timecount]][2,j]==1){which(smat[[timecount]][2,]==1)}else{logical(0)}#these elements are sympatric and in S2
     		 	tempS1 =  ifelse(length(elmsS1)>0,S1*(mean(traitMat[elmsS1,k])-traitMat[j,k])*segsize,0)
     		 	tempS2 =  ifelse(length(elmsS2)>0,S2*(mean(traitMat[elmsS2,k])-traitMat[j,k])*segsize,0)
    		    temp2 = sqrt(sig2*segsize)
          	traitMat[j,k+1]<-traitMat[j,k] +(tempS1+tempS2)/sum(smat[[timecount]][,j]) +rnorm(1,0,temp2)
			if((k!=segmentSize[i])&&(j==branchesPresent[i])){timecount= ifelse(round(((k*segsize)+nodeDist[i]),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
		###loop for last segment size (to preserve exact branch lengths)
		if(k==segmentSize[i] && nodeDiff[i]%%segsize!=0){
			segsizeT= nodeDiff[i]%%segsize
     		 		tempS1B =  ifelse(length(elmsS1)>0,S1*(mean(traitMat[elmsS1,k])-traitMat[j,k])*segsizeT,0)
     		 		tempS2B =  ifelse(length(elmsS2)>0,S2*(mean(traitMat[elmsS2,k])-traitMat[j,k])*segsizeT,0)
       		temp2B = sqrt(sig2*segsizeT)
           	traitMat[j,k+1]<-traitMat[j,k] +(tempS1B+tempS2B)/sum(smat[[timecount]][,j]) +rnorm(1,0,temp2B)
		
		 	if(j==branchesPresent[i]){timecount= ifelse(round((((k-1)*segsize)+nodeDist[i]+segsizeT),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
			}
       	}	        	
	}
      masterbranch[[i]] = traitMat
    }
  out = as.vector(masterbranch[[phylo$Nnode]][, tail(segmentSize,n=1)+1])
  names(out)<-unlist(nat[[i]])
}

if(model%in%c("DDexp","DDlin")){
	if(model=="DDexp"){
		if((sig2/(sig2*exp(S1*rowSums(smat[[length(smat)]])[1])))>=1e3){warning("r1 parameter leads to sig2 at tips that is 1000x smaller than root value, consider changing")}
		if((sig2/(sig2*exp(S2*rowSums(smat[[length(smat)]])[2])))>=1e3){warning("r2 parameter leads to sig2 at tips that is 1000x smaller than root value, consider changing")}
	}
	if(model=="DDlin"){
		test1<-sig2+(S1*rowSums(smat[[length(smat)]])[1]) 
		if(test1<0){warning("b1 parameter leads to negative sig2 values")}
		test2<-sig2+(S2*rowSums(smat[[length(smat)]])[2]) 
		if(test2<0){warning("b2 parameter leads to negative sig2 values")}
	}	

	if(!all(colnames(smat[[length(smat)]])==nat[[length(nat)]])){stop("order of S.matrix species is incorrect; perhaps this is a sorted matrix?")}	
    newDist<-SMatrix$times
    newDiff<-SMatrix$spans
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
	            
	            if(model=="DDexp"){
       		 	tempS1 =  ifelse(smat[[timecount]][1,j]>0,exp(S1*sum(smat[[timecount]][1,])),0)
       		 	tempS2 =  ifelse(smat[[timecount]][2,j]>0,exp(S2*sum(smat[[timecount]][2,])),0)
            	traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(sig2*((tempS1+tempS2)/sum(smat[[timecount]][,j]))*segsize))
				} else{
				tempS1 = ifelse(smat[[timecount]][1,j]>0,(S1*sum(smat[[timecount]][1,])),0)
				tempS2 = ifelse(smat[[timecount]][2,j]>0,(S2*sum(smat[[timecount]][2,])),0)
				hold= (((sig2+tempS1)*smat[[timecount]][1,j])+((sig2+tempS2)*smat[[timecount]][2,j]))/sum(smat[[timecount]][,j])
				if(hold<=0){stop("rates lead to negative rates")}
            	traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(hold*segsize))
				}
				
				
				
				if((k!=segmentSize[i])&&(j==branchesPresent[i])){timecount= ifelse(round(((k*segsize)+nodeDist[i]),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
				###loop for last segment size (to preserve exact branch lengths)
				if(k==segmentSize[i] && nodeDiff[i]%%segsize!=0){
					segsizeT= nodeDiff[i]%%segsize
					
					if(model=="DDexp"){
	            	traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(sig2*((tempS1+tempS2)/sum(smat[[timecount]][,j]))*segsizeT))
					} else{
            		traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(((((sig2+tempS1)*smat[[timecount]][1,j])+((sig2+tempS2)*smat[[timecount]][2,j]))/sum(smat[[timecount]][,j]))*segsizeT))
					}
				 	if(j==branchesPresent[i]){timecount= ifelse(round((((k-1)*segsize)+nodeDist[i]+segsizeT),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
					}
	        	}	        	
			}
        masterbranch[[i]] = traitMat
      }
    out = as.vector(masterbranch[[phylo$Nnode]][, tail(segmentSize,n=1)+1])
    names(out)<-unlist(nat[[i]])

}

return(out)

}
