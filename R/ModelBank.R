##################################################
#    Bank of Classic 1D Phenotypic Models
##################################################

createModel <- function(tree, keyword){
    
    if(keyword == "BM"){

        comment <- "Brownian Motion model with linear drift.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = d dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "d", "sigma")
        params0 <- c(0,0,0,1)

        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(params[3]*vectorU)
            matrixGamma <- function(t) return(params[4]*diag(vectorU))
            matrixA <- diag(0, sizeU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[2]>=0 && params[4]>=0)
        
        model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=getMatrixOfCoalescenceTimes(tree))

    }else if(keyword == "BMbis"){

        comment <- "Brownian Motion model with linear drift.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = d dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "d", "sigma")
        params0 <- c(0,0,0,1)

        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(params[3]*vectorU)
            matrixGamma <- function(t) return(params[4]*diag(vectorU))
            matrixA <- diag(0, sizeU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[2]>=0 && params[4]>=0)
        
        model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "BM_from0"){

        comment <- "Brownian Motion model with linear drift.\nStarts with two lineages having the same value X_0 = (0,0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = d dt + sigma dW_t"
        paramsNames <- c("d", "sigma")
        params0 <- c(0,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list( mean=c(0), var=matrix(c(0)) ) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(params[1]*vectorU)
            matrixB <- function(t) return(params[2]*diag(vectorU))
            matrixA <- diag(0, sizeU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixB))
        }
        
        constraints <- function(params) return(params[2]>=0)
        
        model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, comment=comment, matrixCoalescenceTimes=getMatrixOfCoalescenceTimes(tree))

    }else if(keyword == "BM_from0_driftless"){

        comment <- "Brownian Motion model without drift.\nStarts with two lineages having the same value X_0 = (0,0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma dW_t"
        paramsNames <- c("sigma")
        params0 <- c(1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list( mean=c(0), var=matrix(c(0)) ) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(0*vectorU)
            matrixGamma <- function(t) return(params[1]*diag(vectorU))
            matrixA <- diag(0, sizeU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[1]>=0)
        
        model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=getMatrixOfCoalescenceTimes(tree))

    }else if(keyword == "OU"){

        comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "psi", "theta", "sigma")
        params0 <- c(0,0,1,0,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(params[3]*params[4]*vectorU)
            matrixGamma <- function(t) return(params[5]*diag(vectorU))
            matrixA <- params[3]*diag(vectorU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[5]>=0 && params[3]!=0)
        
        model <- new(Class="PhenotypicOU", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=getMatrixOfCoalescenceTimes(tree))

    }else if(keyword == "OUbis"){

        comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "psi", "theta", "sigma")
        params0 <- c(0,0,1,0,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(params[3]*params[4]*vectorU)
            matrixGamma <- function(t) return(params[5]*diag(vectorU))
            matrixA <- params[3]*diag(vectorU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[5]>=0 && params[3]!=0)
        
        model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "OUter"){

        comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "psi", "theta", "sigma")
        params0 <- c(0,0,1,0,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(params[3]*params[4]*vectorU)
            matrixGamma <- function(t) return(params[5]*diag(vectorU))
            matrixA <- params[3]*diag(vectorU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[5]>=0 && params[3]!=0)
        
        model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "OU_from0"){

        comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 = (0,0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
        paramsNames <- c("psi", "theta", "sigma")
        params0 <- c(0.01,0,1)
        # This model requires the following list of parameters, in this order : psi, theta, sigma
        # One trait in each lineage
        # starts with two lineages with the same value X_0 ~ Normal(0,0)
        # dX_t = psi(theta- X_t) dt + sigma dW_t
        # Tree is an object of type "phylo" (cf. ape)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(0), var=matrix(c(0))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages  
            vectorU <- description$livingLineages
            vectorA <- function(t) return(params[1]*params[2]*vectorU)
            matrixGamma <- function(t) return(params[3]*diag(vectorU))
            matrixA <- params[1]*diag(vectorU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[3]>=0 && params[1]!=0)

        model <- new(Class="PhenotypicOU", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=getMatrixOfCoalescenceTimes(tree))

    }else if(keyword == "EB"){

        comment <- "Early-Burst model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(-1/2rt) dW_t"
        paramsNames <- c("m0", "v0", "sigma0", "r")
        params0 <- c(0,0,100,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages
            vectorU <- description$livingLineages 
            vectorA <- function(t) return(rep(0, sizeU))
            matrixGamma <- function(t) return(params[3]*exp(-params[4]*t)*diag(vectorU))
            matrixA <- diag(0, sizeU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[3]>0 && params[4]<=1 && params[4] >= -1)
        
        model <- new(Class="PhenotypicEB", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=getMatrixOfCoalescenceTimes(tree))

    }else if(keyword == "EBbis"){

        comment <- "Early-Burst model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(-1/2rt) dW_t"
        paramsNames <- c("m0", "v0", "sigma0", "r")
        params0 <- c(0,0,100,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            description <- describeOnePeriod(i, periodizing, tree)
            sizeU <- description$nLineages
            vectorU <- description$livingLineages 
            vectorA <- function(t) return(rep(0, sizeU))
            matrixGamma <- function(t) return(params[3]*exp(-params[4]*t)*diag(vectorU))
            matrixA <- diag(0, sizeU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[3]>0 && params[4]<=1 && params[4] >= -1)
        
        model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "PM"){

        comment <- "Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression."
        paramsNames <- c("m0", "v0", "theta", "psi", "S", "sigma")
        params0 <- c(0,0,0,0.2,0.5,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
            
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(params[3]*params[4]*vectorU)
            matrixGamma <- function(t) return(params[6]*diag(vectorU))
            matrixA <- (params[4]+params[5])*diag(vectorU) - (params[5]/sum(vectorU)) * outer(vectorU,vectorU) 
              
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[6]>=0)
        
        model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "PMbis"){

        comment <- "Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression."
        paramsNames <- c("m0", "v0", "theta", "psi", "S", "sigma")
        params0 <- c(0,0,0,0.2,0.5,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
            
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(params[3]*params[4]*vectorU)
            matrixGamma <- function(t) return(params[6]*diag(vectorU))
            matrixA <- (params[4]+params[5])*diag(vectorU) - (params[5]/sum(vectorU)) * outer(vectorU,vectorU) 
              
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[6]>=0)
        
        model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "PM_IMACS"){

        comment <- "Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression."
        paramsNames <- c("m0", "v0", "theta", "psi", "S", "sigma")
        params0 <- c(0,0,0,0.2,0.5,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
            
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
              
            return(list(v=vectorU))
        }

        constraints <- function(params) return(params[2]>=0 && params[6]>=0)
        
        model <- new(Class="PhenotypicADiagIMACS", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "PM_OUless"){

        comment <- "Simplified Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression, without the OU term."
        paramsNames <- c("m0", "v0", "S", "sigma")
        params0 <- c(0,0,0.5,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)

        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
            
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(0*vectorU)
            matrixGamma <- function(t) return(params[4]*diag(vectorU))
            matrixA <- params[3]*diag(vectorU) - (params[3]/sum(vectorU)) * outer(vectorU,vectorU) 
              
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[4]>=0)
        
        model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else{
        stop("Keyword does not correspond to any model in the model bank")
    }


    return(model)
}



##################################################
#    Describe the periods on a 'phylo' tree
##################################################

getMatrixOfCoalescenceTimes <- function(tree){
    startingTimes <- getStartingTimes(tree)
    # we get the matrix giving the number of the mrca node to each pair of tips
    matrixMRCA <- mrca(tree)
    matrixCoalescenceTimes <- diag(0, length(matrixMRCA[,1]))

    # We fill the matrix of coalescence times
    for(i in 1:length(matrixMRCA[,1])){
        for(j in i:length(matrixMRCA[1,])){
            # if the mrca is an internal node, we take its date from the root
            if( matrixMRCA[i,j] %in% tree$edge[,1]){
                index <- which(tree$edge[,1]==matrixMRCA[i,j])[1]
                matrixCoalescenceTimes[i,j] <- startingTimes[index]
                matrixCoalescenceTimes[j,i] <- startingTimes[index]
            # if the mrca is a tip, we take its death time
            }else{
                index <- which(tree$edge[,2]==matrixMRCA[i,j])
                matrixCoalescenceTimes[i,j] <- startingTimes[index]+tree$edge.length[index]
                matrixCoalescenceTimes[j,i] <- startingTimes[index]+tree$edge.length[index]
            }
        }
    }
    return(matrixCoalescenceTimes)
}

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

removeIdenticalEntries <- function(vector){
    # Remove identical entries on a vector in ascending order
    new_vector <- c(vector[1])
    for(i in 1:(length(vector)-1)){
        if(vector[i] < vector[i+1]-1e-5){
            new_vector <- c(new_vector, vector[i+1])
        }
    }
    return(new_vector)
}

isATip <- function(tree, branch_number){
    return(!(tree$edge[branch_number,2] %in% tree$edge[,1]))
}

periodizeOneTree <- function(tree){
    # Returns 3 vectors giving 
    # 1) the periods of the tree, 
    # 2) the starting times of all branches in the tree 
    # 3) the death time of all branches in the tree
    
    startingTimes <- getStartingTimes(tree)
    endTimes <- startingTimes + tree$edge.length
    all_time_events <- sort(c(startingTimes, endTimes))
    periods <- removeIdenticalEntries(all_time_events)
    
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
            deadBranches <- which(tau_i == periodizing$endTimes)
            numbersCopy[i] <- labelingLineages[ deadBranches[1] ]
            numbersPaste[i] <- 0
        }
        numbersLineages[i] <- n
    }

    permutationLabels <- labelingLineages[!(periodizing$endTimes %in% periodizing$startingTimes)]
    labeling <- tree$tip.label[order(permutationLabels)]
    
    return(list(copy=numbersCopy, paste=numbersPaste, nLineages=numbersLineages, labeling=labeling))
}

############################################
# Old stuff with 'ever changing' labeling, slower
# The following should be deleted at some point
############################################

describeOnePeriod <- function(i, periodizing, tree){
    
    nBranch <- length(periodizing$startingTimes)
    numbering <- rep(0, times=nBranch)
    IDLineage <- 1
    numberCopied <- 0
    livingLineages <- c()
        
    for(k in 1:nBranch){

        # If the branch is running, we assign it a number and we put a 1 in the livingLineages vector
        if( periodizing$startingTimes[k] <= periodizing$periods[i] && periodizing$endTimes[k] > periodizing$periods[i] ){
            numbering[k] <- IDLineage
            IDLineage <- IDLineage + 1
            livingLineages <- c(livingLineages, 1)
            
            # If two branches start at the beginning of the period, we copy the ID of the first
            if(periodizing$startingTimes[k] == periodizing$periods[i] && numberCopied == 0){
                numberCopied <- numbering[k]
                numberPasted <- numbering[k] + 1
            }
            
        # If this is a tip already dead on the considered period, we assign it a number and we put a zero in livingLineages
        }else if( periodizing$startingTimes[k] <= periodizing$periods[i] && isATip(tree,k) ){
            numbering[k] <- IDLineage
            IDLineage <- IDLineage + 1
            livingLineages <- c(livingLineages, 0)
            
            # If a tip branch is stopping at the beginning of the period
            if(periodizing$endTimes[k] == periodizing$periods[i] ){
                # the branch number is collected
                numberCopied <- numbering[k]
                numberPasted <- 0
            }
        }
        
    }
    
    return(list(copy=numberCopied, paste=numberPasted, nLineages=IDLineage-1, livingLineages=livingLineages))
}

endOfPeriods2 <- function(periodizing, tree){
    # Returns the list of branching or dying lineages at the beginning of each period : copy
    # Together with the list of places where the new lineage is inserted (or zero if a lineage dies) : paste
    # And the number of lineages on the focal period : nLineages
    
    nBranch <- length(periodizing$startingTimes)
    nPeriods <- length(periodizing$periods)
    
    numbersCopy <- rep(0, times=nPeriods)
    numbersPaste <- rep(0, times=nPeriods)
    numbersLineages <- rep(0, times=nPeriods)
    
    for(i in 1:nPeriods){
        specificPeriod <- describeOnePeriod(i, periodizing, tree)
        numbersCopy[i] <- specificPeriod$copy
        numbersPaste[i] <- specificPeriod$paste
        numbersLineages[i] <- specificPeriod$nLineages
    }
    
    return(list(copy=numbersCopy, paste=numbersPaste, nLineages=numbersLineages))
}
