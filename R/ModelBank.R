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
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(params[3]*vectorU)
            matrixGamma <- function(t) return(params[4]*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[2]>=0 && params[4]>=0)
        
        model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))

    }else if(keyword == "BMbis"){

        comment <- "Brownian Motion model with linear drift.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = d dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "d", "sigma")
        params0 <- c(0,0,0,1)

        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(params[3]*vectorU)
            matrixGamma <- function(t) return(params[4]*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
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
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(params[1]*vectorU)
            matrixGamma <- function(t) return(params[2]*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[2]>=0)
        
        model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))

    }else if(keyword == "BM_from0_driftless"){

        comment <- "Brownian Motion model without drift.\nStarts with two lineages having the same value X_0 = (0,0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma dW_t"
        paramsNames <- c("sigma")
        params0 <- c(1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list( mean=c(0), var=matrix(c(0)) ) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(0*vectorU)
            matrixGamma <- function(t) return(params[1]*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[1]>=0)
        
        model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))

    }else if(keyword == "OU"){

        comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "psi", "theta", "sigma")
        params0 <- c(0,0,1,0,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(params[3]*params[4]*vectorU)
            matrixGamma <- function(t) return(params[5]*diag(vectorU))
            matrixA <- params[3]*diag(vectorU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[5]>=0 && params[3]!=0)
        
        model <- new(Class="PhenotypicOU", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))

    }else if(keyword == "OUbis"){

        comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
        paramsNames <- c("m0", "v0", "psi", "theta", "sigma")
        params0 <- c(0,0,1,0,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
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
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
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
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(params[1]*params[2]*vectorU)
            matrixGamma <- function(t) return(params[3]*diag(vectorU))
            matrixA <- params[1]*diag(vectorU)
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }
        
        constraints <- function(params) return(params[3]>=0 && params[1]!=0)

        model <- new(Class="PhenotypicOU", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))

    }else if(keyword == "ACDC"){

        comment <- "ACcelerating or DeCelerating rate of evolution.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(rt) dW_t"
        paramsNames <- c("m0", "v0", "sigma0", "r")
        params0 <- c(0,0,100,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(rep(0, length(vectorU)))
            matrixGamma <- function(t) return(params[3]*exp(params[4]*t)*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[3]>0)
        
        model <- new(Class="PhenotypicACDC", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))

    }else if(keyword == "ACDCbis"){

        comment <- "ACcelerating or DeCelerating rate of evolution.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(rt) dW_t"
        paramsNames <- c("m0", "v0", "sigma0", "r")
        params0 <- c(0,0,100,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(rep(0, length(vectorU)))
            matrixGamma <- function(t) return(params[3]*exp(params[4]*t)*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[3]>0)
        
        model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "DD"){

        comment <- "Diversity-dependent model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(r n_t) dW_t"
        paramsNames <- c("m0", "v0", "r", "sigma0")
        params0 <- c(0,0,100,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(rep(0, length(vectorU)))
            matrixGamma <- function(t) return(params[4]*exp(params[4]*sum(vectorU))*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[4]>0)
        
        model <- new(Class="PhenotypicDD", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, comment=comment, matrixCoalescenceJ=getMatrixCoalescenceJ(tree, periodizing$periods), nLivingLineages=eventEndOfPeriods$nLivingLineages)

    }else if(keyword == "DDbis"){

        comment <- "Diversity-dependent model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(r n_t) dW_t"
        paramsNames <- c("m0", "v0", "r", "sigma0")
        params0 <- c(0,0,100,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
        
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(rep(0, length(vectorU)))
            matrixGamma <- function(t) return(params[4]*exp(params[3]*sum(vectorU))*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[2]>=0 && params[4]>0)
        
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
              
            return(list(u=vectorU, OU=TRUE))
        }

        constraints <- function(params) return(params[2]>=0 && params[6]>=0)
        
        model <- new(Class="PhenotypicPM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

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
        
        model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "PMter"){

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

    }else if(keyword == "PM_OUless"){

        comment <- "Simplified Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression, without the OU term."
        paramsNames <- c("m0", "v0", "S", "sigma")
        params0 <- c(0,0,0.5,1)
        
        periodizing <- periodizeOneTree(tree)
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)

        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
            
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
              
            return(list(u=vectorU, OU=FALSE))
        }

        constraints <- function(params) return(params[2]>=0 && params[4]>=0)
        
        model <- new(Class="PhenotypicPM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)

    }else if(keyword == "PM_OUlessbis"){

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

getMatrixCoalescenceJ <- function(tree, periods){
    # The entry (k,l) of the matrix is the index j such that tau_j = t_{k,l}
    matrixCoalescenceTimes <- findMRCA(tree, type="height")
    n <- length(matrixCoalescenceTimes[,1])
    matrixCoalescenceJ <- diag(0, n)
    for(k in 1:n){
        for(l in 1:n){
            matrixCoalescenceJ[k,l] <- which(periods == matrixCoalescenceTimes[k,l])
        }
    }

    return(matrixCoalescenceJ)
}

isATip <- function(tree, branch_number){
    return(!(tree$edge[branch_number,2] %in% tree$edge[,1]))
}

periodizeOneTree <- function(tree){
    # Returns 3 vectors giving 
    # 1) the periods of the tree, 
    # 2) the starting times of all branches in the tree 
    # 3) the death time of all branches in the tree
    
    nodeheight <- nodeHeights(tree)
    startingTimes <- nodeheight[,1]
    endTimes <- nodeheight[,2]
    all_time_events <- sort(c(startingTimes, endTimes))
    # the following removes identical entries in the vector
    periods <- unique(all_time_events)
    
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
    numbersLivingLineages <- rep(0, times=nPeriods)

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
    numbersLivingLineages[1] <- n
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
            numbersLivingLineages[i] <- numbersLivingLineages[i-1]+1
        # Else, tau_i is only a death time of one or many terminal branches.
        }else{
            deadBranches <- which(tau_i == periodizing$endTimes)
            numbersCopy[i] <- labelingLineages[ deadBranches[1] ]
            numbersPaste[i] <- 0
            numbersLivingLineages[i] <- numbersLivingLineages[i-1]-1
        }
        numbersLineages[i] <- n
    }

    permutationLabels <- labelingLineages[!(periodizing$endTimes %in% periodizing$startingTimes)]
    labeling <- tree$tip.label[order(permutationLabels)]
    
    return(list(copy=numbersCopy, paste=numbersPaste, nLineages=numbersLineages, labeling=labeling, nLivingLineages=numbersLivingLineages))
}

getLivingLineages <- function(i, eventEndOfPeriods){
    
    livingLineages <- rep(1, times=eventEndOfPeriods$nLineages[i])
    deads <- eventEndOfPeriods$copy[1:i][eventEndOfPeriods$paste[1:i] == 0]
    livingLineages[deads] <- 0
    
    return(livingLineages)
}

