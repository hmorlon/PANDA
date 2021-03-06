##################################################
#    Bank of Classic 1D Phenotypic Models
##################################################
#r.object here refers to a list (equal to the length of the geo.object, if incorporating biogeography) of matrices with two rows (one for each competitive regime) and the same number of columns as there are species

.createModel_DDlin_multi_geo_ME <- function(tree,geo.object,r.object){    

        comment <- "Diversity dependent exponential model with biogeography and two slope regimes and measurement error."
        paramsNames <- c("m0", "logsigma0", "r1", "r2","lognuisance")
        params0 <- c(0,log(1),0,0,log(1))


        periodizing <- periodizeOneTree_multigeo(tree,geo.object) 
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) )        
        
        aAGamma <- function(i, params){
        	rmat<-r.object$S.matrix[[i]]
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(rep(0, length(vectorU)))
            
        	R1.mat <- geo.object$geography.object[[i]]*outer(rmat[1,],rmat[1,])
			#diag(R1.mat)<-1
			R2.mat <- geo.object$geography.object[[i]]*outer(rmat[2,],rmat[2,])
			#diag(R2.mat)<-1
			nij1<-colSums(R1.mat)
			nij2<-colSums(R2.mat)

            #matrixGamma <- function(t) return(sqrt((exp(params[2])^2)+(params[3]*colSums(geo.object$geography.object[[i]])))*diag(vectorU))
            matrixGamma <- function(t) return(sqrt((exp(params[2])^2)+1/colSums(rmat)*((params[3]*nij1)+(params[4]*nij2)))*diag(vectorU))
			#matrixGamma <- function(t) return(((exp(params[2])+diag(nij1)*params[3])*diag(diag(R1.mat))/colSums(rmat))+((exp(params[2])+diag(nij2)*params[4])*diag(diag(R2.mat))/colSums(rmat)))
            
            matrixA <- diag(0, length(vectorU))
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(((exp(params[2])^2)+(params[3]*max(vapply(geo.object$geography.object,function(x)max(rowSums(x)),1))) > 0) &&((exp(params[2])^2)+(params[4]*max(vapply(geo.object$geography.object,function(x)max(rowSums(x)),1))) > 0))
        
        model <- new(Class="PhenotypicADiag", name="DDlinM", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)


    return(model)
}



##################################################
#    Describe the periods on a 'phylo' tree
##################################################

getStartingTimes <- function(tree){
    # Returns a vector giving the starting time for each branch of a tree in format "phylo"
    
    nBranch = length(tree$edge.length)
    starting_times <- rep(0, times=nBranch)
    
    # we add progressively for each branch the length of all parent branches in the vector "starting_times"
    for(n1 in 1:nBranch){
        n2 <- n1 + 1
        while(n2 <= nBranch){
            if(tree$edge[n2,1]==tree$edge[n1,2]){
                starting_times[n2] <- starting_times[n1] + tree$edge.length[n1]
            }
            n2 <- n2+1
        }
    }
    
    return(starting_times)
}


isATip <- function(tree, branch_number){
    return(!(tree$edge[branch_number,2] %in% tree$edge[,1]))
}

periodizeOneTree_multigeo <- function(tree,geo.object){
    # Returns 3 vectors giving 
    # 1) the periods of the tree, 
    # 2) the starting times of all branches in the tree 
    # 3) the death time of all branches in the tree
    hold<-nodeHeights(tree)
    startingTimes <- hold[,1]
    endTimes <- hold[,2]
    all_time_events <- sort(c(startingTimes, endTimes))
    
    
    nodetimes=max(branching.times(tree))-sort(branching.times(tree),decreasing=TRUE)
	extv<-vapply(geo.object$geography.object,function(x)dim(x)[2],1)
	outv<-c(1)
	for(i in 2:length(extv)){
		if(extv[i]!=extv[i-1]){
			outv<-c(outv,i)
		}}
	
	chg.times=which(!1:length(geo.object$times)%in%c(outv,length(geo.object$times)))
	periods=sort(c(geo.object$times[chg.times],unique(startingTimes),max(endTimes)))
    return(list(periods=periods, startingTimes=startingTimes, endTimes=endTimes))
}

endOfPeriods <- function(periodizing, tree){
    # Returns the list of branching or dying lineages at the beginning of each period : copy
    # Together with the list of places where the new lineage is inserted (or zero if a lineage dies) : paste
    # And the number of lineages on the focal period : nLineages
    # The rule is : at each branching point, the first of the two new branches is assigned its mother label, and the new one takes the last label (n, where n is the number of lineages at that time)
    
    nBranch <- length(periodizing$startingTimes)
    nPeriods <- length(periodizing$periods)
    
    numbersCopy <- rep(0, times=nPeriods)
    numbersPaste <- rep(0, times=nPeriods)
    numbersLineages <- rep(0, times=nPeriods)

    # We initialize the labeling of branches in the tree
    labelingLineages <- rep(0, times=nBranch)
    initialBranches <- periodizing$startingTimes[periodizing$startingTimes==0]
    if(length(initialBranches) == 1){
        labelingLineages[1] <- 1
        n <- 1
    }else{
        labelingLineages[periodizing$startingTimes==0] <- c(1,2)
        n <- 2
    }
    numbersLineages[1] <- n
    numbersCopy[1] <- 1
    numbersPaste[1] <- 2
    
    for(i in 2:nPeriods){
        tau_i <- periodizing$periods[i]
        newBranches <- which(tau_i == periodizing$startingTimes)
        # If tau_i is a birth time on the tree
        if(length(newBranches) == 2){
            n <- n+1
            labelingLineages[newBranches[1]] <- labelingLineages[newBranches[1]-1]
            labelingLineages[newBranches[2]] <- n
            numbersCopy[i] <- labelingLineages[newBranches[1]-1]
            numbersPaste[i] <- n
        # Else, tau_i is only a death time of one or many terminal branches.
        }else{
            numbersCopy[i] <- 0
 			#deadBranches <- which(tau_i == periodizing$endTimes)
            #numbersCopy[i] <- labelingLineages[ deadBranches[1] ]
            numbersPaste[i] <- 0
        }
        numbersLineages[i] <- n
    }

    permutationLabels <- labelingLineages[!(periodizing$endTimes %in% periodizing$startingTimes)]
    labeling <- tree$tip.label[order(permutationLabels)]
    
    return(list(copy=numbersCopy, paste=numbersPaste, nLineages=numbersLineages, labeling=labeling))
}


getLivingLineages <- function(i, eventEndOfPeriods){
    
    livingLineages <- rep(1, times=eventEndOfPeriods$nLineages[i])
    deads <- eventEndOfPeriods$copy[1:i][eventEndOfPeriods$paste[1:i] == 0]
    livingLineages[deads] <- 0
    
    return(livingLineages)
}
    