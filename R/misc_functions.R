#code from Nick Matzke's code on BioGeoBEARS wiki
.events_txt_list_into_events_table<-function(events_txt_list, trtable=NULL, recalc_abs_ages=TRUE)
	{
	
	if (is.null(events_txt_list))
		{
		errortxt = paste("\nWARNING in events_txt_list_into_events_table(): your events_txt_list has NO events!\n\nThis means your tree has NO d/e/a events across the whole tree.\nThis is *expected* e.g. if you inferred d=e=0 under DEC+J. Input a list of '' or NA to avoid this error.\n\n", sep="")
		cat(errortxt)
		errortxt2 = paste("events_txt_list_into_events_table() is returning NULL which will might cause issues later.\n\n", sep="")
		cat(errortxt2)
		return(NULL)
		}
	
	# Convert NAs to "none"
	events_txt_list[is.na(events_txt_list)] = "none"
	
	# Remove lines with no events or NA:
	noneTF = events_txt_list == "none"
	keepTF = (noneTF == FALSE)
	events_txt_list = events_txt_list[keepTF]
	
	
	# If no anagenetic events, return NULL
	if (length(events_txt_list) == 0)
		{
		events_table = NULL
		return(events_table)
		}


	# Include the trtable, if that is input
	if (length(trtable) > 0)
		{
		trtable_subset = NULL
		}

	
	# Convert the events text back into a table:
	tmptable = NULL
	for (i in 1:length(events_txt_list))
		{
		#print(events_txt_list)
		tmptable_rows = .events_txt_into_events_table(events_txt_list[i])
		rownums_in_trtable = as.numeric(tmptable_rows$nodenum_at_top_of_branch)
		#print(tmptable_rows)
		num_newrows = nrow(tmptable_rows)
		tmptable = rbind(tmptable, tmptable_rows)
		} # END for (i in 1:length(events_txt_list))
	events_table = .dfnums_to_numeric(.adf2(tmptable))
	names(events_table) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")	
	
	return(events_table)
	}
	
### adf2 from CRAN version of BioGeoBEARS

.adf2<-function (x) 
{
    rownames = 1:nrow(x)
    return(as.data.frame(x, row.names = rownames, stringsAsFactors = FALSE))
}

### dfnums_to_numeric from CRAN version of BioGeoBEARS
.dfnums_to_numeric<-	function (dtf, max_NAs = 0.5, printout = FALSE, roundval = NULL) 
{
    dtf_classes = .cls.df(dtf, printout = FALSE)
    dtf_names = names(dtf)
    numcols = ncol(dtf)
    cls_col_list = c()
    for (i in 1:numcols) {
        cls_col = NA
        cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], 
            "')", sep = "")
        eval(parse(text = cmdstr))
        cls_col_list[i] = cls_col
    }
    for (i in 1:numcols) {
        if (cls_col_list[i] == "list") {
            (next)()
        }
        if (cls_col_list[i] != "numeric") {
            newcol = NA
            cmdstr = paste("newcol = as.numeric(as.character(dtf$'", 
                dtf_names[i], "'))", sep = "")
            suppressWarnings(eval(parse(text = cmdstr)))
            if (sum(is.na(newcol)) < (max_NAs * length(newcol))) {
                cmdstr = paste("dtf$'", dtf_names[i], "' = newcol", 
                  sep = "")
                suppressWarnings(eval(parse(text = cmdstr)))
                if (!is.null(roundval)) {
                  cmdstr = paste("dtf$'", dtf_names[i], "' = round(dtf$'", 
                    dtf_names[i], "', digits=roundval)", sep = "")
                  suppressWarnings(eval(parse(text = cmdstr)))
                }
            }
        }
    }
    tmp_classes = .cls.df(dtf)
    dtf_classes$newclasses = tmp_classes[, ncol(tmp_classes)]
    if (printout) {
        cat("\n")
        cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, ") reports: dataframe 'dtf_classes' has ", 
            nrow(dtf_classes), " rows, ", ncol(dtf_classes), 
            " columns.\n", sep = "")
        cat("...names() and classes() of each column below...\n", 
            sep = "")
        cat("\n")
        print(dtf_classes)
    }
    return(dtf)
}


###cls.df from CRAN version of BioGeoBEARS
.cls.df<-function (dtf, printout = FALSE) 
{
    if (inherits(dtf, "matrix")) {
        dtf = as.data.frame(dtf, stringsAsFactors = FALSE)
    }
    dtf_names = names(dtf)
    numcols = ncol(dtf)
    cls_col_list = c()
    for (i in 1:numcols) {
        cls_col = NA
        cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], 
            "')", sep = "")
        eval(parse(text = cmdstr))
        cls_col_list[i] = cls_col
    }
    colnum = 1:numcols
    dtf_classes = cbind(colnum, dtf_names, cls_col_list)
    dtf_classes = data.frame(dtf_classes, row.names = colnum)
    if (printout) {
        cat("\n")
        cat("cls.df(dtf) reports: dataframe 'dtf' has ", nrow(dtf), 
            " rows, ", numcols, " columns.\n", sep = "")
        cat("...names() and classes() of each column below...\n", 
            sep = "")
        cat("\n")
        print(dtf_classes)
        cat("\n")
    }
    return(dtf_classes)
}

#events_txt_into_events_table from Nick Matzke's code on BioGeoBEARS wiki
.events_txt_into_events_table<-function(branch_events_txt)
	{
	words = strsplit(branch_events_txt, split=";")[[1]]
	
	events_table_for_branch = t(sapply(X=words, FUN=.event_txt_into_events_row))
	row.names(events_table_for_branch) = NULL
	events_table_for_branch
	
	events_table_for_branch = .adf2(events_table_for_branch)
	events_table_for_branch
	names(events_table_for_branch) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")
	
	return(events_table_for_branch)
	}

#code from Nick Matzke's code on BioGeoBEARS wiki
.event_txt_into_events_row<-function(word)
	{

	split_key_item<-function(word2)
		{
		output_pair = c("", "")
		words3 = strsplit(word2, split=":")[[1]]
		numwords = length(words3)
		
		output_pair[1:numwords] = words3[1:numwords]
		
		return(output_pair)
		}


	words2 = strsplit(word, split=",")[[1]]
	output = sapply(X=words2, FUN=split_key_item)
	tmprow = matrix(data=output[2,], nrow=1)
	return(tmprow)
	}

.phyOU<-function(phy,alpha){
             if(alpha<=.Machine$double.eps) return(phy) # reduce to BM
             
             times <- branching.times(phy)
             names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             Tmax<-times[1]
             phy2<-phy
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[which(names(times) == phy$edge[i, 1])]
               t1 <- max(times) - age
               t2 <- t1+bl
               phy2$edge.length[i] <- (1/(2*alpha))*exp(-2*alpha * (Tmax-t2)) * (1 - exp(-2 * alpha * t2)) -
                 (1/(2*alpha))*exp(-2*alpha * (Tmax-t1)) * (1 - exp(-2 * alpha * t1))
             }
             phy <- phy2
             return(phy)
           }

.create.function.list.EB1<-function(geo.simmap){
	if(is.null(geo.simmap)){stop('provide geo.simmap')}
	states<-colnames(geo.simmap$mapped.edge)
	funlist<-list()	
	for(i in 1:length(states)){
		eval(parse(text=paste0("funlist[[",i,"]]<-function(x,df,times){;return(x);}")))
		}
			
	return(funlist)

}
           
.create.function.list.EBmulti<-function(regime.simmap){
	if(is.null(regime.simmap)){stop('provide regime.simmap')}
	states<-colnames(regime.simmap$mapped.edge)
	class.object<-try(CreateClassObject(regime.simmap))
	if(inherits(class.object, "try-error")){class.object<-try(CreateClassObject(regime.simmap,rnd=6))}
	if(inherits(class.object, "try-error")){class.object<-CreateClassObject(regime.simmap,rnd=7)}
	
	funlist<-list()	
	for(i in 1:length(states)){
		first.time.bin=min(which(lapply(class.object$class.object,function(x) states[i]%in%x[,2])==TRUE))
		if(first.time.bin==1){ #if state is present at the root
		eval(parse(text=paste0("funlist[[",i,"]]<-function(x,df,times){;return(x);}")))
		} else {
		time.since.root = class.object$times[first.time.bin]
		eval(parse(text=paste0("funlist[[",i,"]]<-function(x,df,times){;return(x-",time.since.root,");}")))
		}
		}
			
	return(funlist)

}

.fit_t_EB<-function(phylo,data,regime.map=NULL,error=NULL, beta=NULL, sigma=NULL, method=c("L-BFGS-B","Nelder-Mead"), upper=list(beta=0,sigma=Inf), lower=-Inf, control=list(maxit=20000), diagnostic=FALSE, echo=FALSE){
	
	
if(is.null(regime.map)){ 	# single slope version 
			
		#convert phylo to simmap object
		hold<-rep("A",length(phylo$tip.label))
		hold[1]<-"B"
		names(hold)<-phylo$tip.label
		smap<-make.simmap(phylo,hold,message=F)
		new.maps<-list()
		for(i in 1:length(phylo$edge.length)){
			new.maps[[i]]<-phylo$edge.length[i]
			names(new.maps[[i]])<-"A"
			}
		new.mapped.edge<- as.matrix(rowSums(smap$mapped.edge))
		colnames(new.mapped.edge)<-"A"	
		smap$maps<-new.maps
		smap$mapped.edge<-new.mapped.edge
		
		#create function
		new_list_function<-.create.function.list.EB1(smap)

		#fit model
		sigma.constraint<-rep(1, dim(smap$mapped.edge)[2])
		beta.constraint<-rep(1, dim(smap$mapped.edge)[2])
		
		up=c(rep(upper$beta,length(beta.constraint)),Inf,0)
				
		out<-.fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=NULL,input.times=NULL,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=up, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	

		
}  else if (!is.null(regime.map)) { # multi-slope version where rates are reset to root state at beginning of new regime

		#first check that regime.map and phylo and data are concordant
		if(!all(as.phylo(phylo)$tip.label == as.phylo(regime.map)$tip.label)) { stop("regime map doesn't match phylogeny")}
		if(length(data) != length(as.phylo(regime.map)$tip.label)) { stop("number of lineages in data and regime map don't match")}
		if(! all (names(data) %in% as.phylo(regime.map)$tip.label)) { stop("names of lineages in data and regime map don't match")}
		if(! all (as.phylo(regime.map)$tip.label %in% names(data)) ) { stop("names of lineages in data and regime map don't match")}
		
		
		new_list_function<-.create.function.list.EBmulti(regime.map)
				
		#fit model
		sigma.constraint<-rep(1, dim(regime.map$mapped.edge)[2])
		beta.constraint<-seq(1,by=1,length.out=dim(regime.map$mapped.edge)[2])
		
		up=c(rep(upper$beta,length(beta.constraint)),Inf,0)
		
		out<-.fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=NULL,class.df=NULL,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=up, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		
}    
		return(out)
}

.fit_t_general <- function(tree, data, fun, class.df, input.times, error=NULL, beta=NULL, sigma=NULL, model=c("exponential","linear"), method=c("L-BFGS-B"), upper=Inf, lower=-Inf, control=list(maxit=20000), diagnostic=TRUE, echo=TRUE, constraint=NULL, return.tree=FALSE) {
  
  if(!inherits(tree,"simmap")==TRUE) stop("For now only simmap-like mapped trees are allowed.","\n")
  old.tree<-tree
  
  # Parameters
  if(!is.null(names(data))) data <- data[tree$tip.label]
  data<-as.matrix(data)
  method=method[1]
  rownames(data)<-tree$tip.label
  model=model[1]
  # Compute node time from the root to the tips
  times<-max(nodeHeights(tree))-nodeHeights(tree)[match(1:tree$Nnode+length(tree$tip),tree$edge[,1]),1]
  
  names(times)<-1:tree$Nnode+length(tree$tip)
  # Set the root to zero
  times<-max(times)-times
  # Max time
  mtot=max(nodeHeights(tree))
  onestate<-ifelse(dim(tree$mapped.edge)[2]==1,TRUE,FALSE) 

  # Number of species
  n=length(tree$tip.label)
  # Number of traits (for future versions)
  k=1
  tree <- reorderSimmap(tree, order="postorder")

    # Number of maps (selective regimes)
    nstates <- dim(old.tree$mapped.edge)[2]
    if(is.null(constraint)){
           number_maps_beta <- number_maps_sigma <- nstates 
           index.user.sigma <- index.user.beta <- 1:number_maps_sigma
        
    }else{ 
        if(is.null(constraint[["sigma"]])) sigConst <- FALSE else sigConst <- TRUE
        if(is.null(constraint[["beta"]])) betConst <- FALSE else betConst <- TRUE
        
        if(sigConst & betConst){ # i) both beta and sigma constrained
            number_maps_beta <- length(unique(constraint$beta[!is.na(constraint$beta)]))
            number_maps_sigma <- length(unique(constraint$sigma))  
            # set the indices
            index.user.sigma <- constraint$sigma
            index.user.beta <- constraint$beta
            
        }else if(!sigConst & betConst){ # ii) only beta constrained
            number_maps_beta <- length(unique(constraint$beta[!is.na(constraint$beta)]))
            number_maps_sigma <- nstates 
            
            index.user.sigma <- 1:number_maps_sigma
            index.user.beta <- constraint$beta
            
        }else if(sigConst & !betConst){ # iii) only sigma constrained
            number_maps_beta <- nstates 
            number_maps_sigma <- length(unique(constraint$sigma))
            
            index.user.sigma <- constraint$sigma
            index.user.beta <- 1:number_maps_beta
        }

    }
    
  
  # Param likelihood contrasts (we are organizing the tree to be postorder for an efficient traversal)
  #ind=reorder(tree,"postorder",index.only=TRUE)
  phy=tree
  #phy$edge.length<-phy$edge.length[ind]
  #phy$edge<-phy$edge[ind,]
  
  # check for simmap like format
  #if(inherits(tree,"simmap")){
  #  phy$mapped.edge<-phy$mapped.edge[ind,]
  #  phy$maps<-phy$maps[ind]
  #}
  #phy <- reorderSimmap(tree, order="pruningwise")
  
  # Random starting value if not provided
  if(is.null(beta)){
    beta=rep(0, number_maps_beta)
  }
  if(is.null(sigma)){
    sigma=rep(sum(pic(data,old.tree)^2)/n, number_maps_sigma)
  }
  
  if(model=="linear"){
    startval=c(beta,log(sigma))
    nbeta=length(beta)
    nsigma=length(sigma)
  }else if(model=="exponential"){
    startval=c(beta,log(sigma))
    nbeta=length(beta)
    nsigma=length(sigma)
  }
  
  # Error estimation?
  if(!is.null(error)){
    ## Index error
    index_error<-sapply(1:n, function(x){ which(phy$edge[,2]==x)})
    startval=c(startval,0.001)
    # to construct a mixed model (refer to l.172 of the code below)
    if(is.numeric(error)){
      error_meas = error^2
    }else{
      error_meas = numeric(n)       
    }
  }
  
  ##--------------Fonction-generale-DD-Env-------------------------------------------##
  
  BranchtransformMAPS<-function(phy,beta,mtot,times,fun,sigma=NULL,model=NULL,errorValue=NULL){
    #Transformations
    tips <- length(phy$tip.label)
    res <- phy
    
    if(model=="exponential"){  
      # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
      f<-function(x, sigma, beta, funInd){sigma*exp(beta*fun[[funInd]](x,df=class.df,times=input.times))}
      
    }else if(model=="linear"){
      # Clim-lin function
      f<-function(x, sigma, beta, funInd){sigma+beta*fun[[funInd]](x,df=class.df,times=input.times)}
     
    }
    
    # Loops over the edges
    for (i in 1:length(phy$edge.length)) {
      
      age <- times[phy$edge[i, 1] - tips] # retrieve the age at the node
      currentmap<-phy$maps[[i]]           # retrieve the corresponding maps on the edge "i"
      indlength<-length(currentmap)       # How many mapping there are?
      tempedge<-numeric(indlength)        # temporary vector for the mappings
      
      # loop pour traverser les "maps"
      for(betaval in 1:indlength){
          
        if(onestate){
            regimenumber=1
        }else{
            regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[betaval])  # retrieve the regimes within maps
        }

        bet<-beta[regimenumber]           # select the corresponding parameter for beta
        sig<-sigma[regimenumber]          # select the corresponding parameter for sigmz
        bl<-currentmap[[betaval]]         # branch length under the current map
        
        int <- try(integrate(f, lower=age, upper=(age + bl), subdivisions=500, rel.tol = .Machine$double.eps^0.05, sigma=sig, beta=bet, funInd=regimenumber), silent=TRUE)
          
           if(inherits(int ,'try-error')){
             warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
             tempbeta <- NA_real_
           } else {
             tempbeta <- int$value
           }
          
        tempedge[betaval] <- tempbeta 
        # on met à jour age parcequ'on va passer au maps suivant pour la lignée i
        # update "age" because we're moving to the next map for lineage i.
        age<-age+bl
      }
      # update branch length
      res$edge.length[i]<-sum(tempedge)
    }
    
    phy<-res
    
    if(!is.null(errorValue)){
      phy$edge.length[index_error]<-phy$edge.length[index_error] + error_meas + errorValue^2
    }
    
    return(phy)
  }
  
  
  ##---------------------------------------------------------------------------------##
  clikCLIM <- function( param, dat, phylo, mtot, times, fun=fun, model, results=FALSE, return.tree=FALSE) {
    
    if(model=="exponential"){
       
        # create vector of parameters
        beta <- numeric(nstates)
        sigma <- numeric(nstates)
        
        # assign values
        beta[] <- c(param[seq_len(nbeta)])[index.user.beta]
        sigma[] <- c(exp(param[nbeta+seq_len(nsigma)]))[index.user.sigma]
        
        # constrain some values to zero (should be only for beta as sigma=0 is undefined)
        beta[is.na(beta)] <- 0
        
      if(!is.null(error)) errorValue <- param[nbeta+nsigma+1] else errorValue <- NULL
      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      if(return.tree){ return(phylo)}
      if(any(is.na(phylo$edge.length)))  return(1000000)
	  if(any(phylo$edge.length<0)) return(1000000)
      LL<-mvLL(phylo,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))
      
    }else{

        # create vector of parameters
        beta <- numeric(nstates)
        sigma <- numeric(nstates)
        
        # assign values
        beta[] <- c(param[seq_len(nbeta)])[index.user.beta]
        sigma[] <- c(exp(param[nbeta+seq_len(nsigma)]))[index.user.sigma]
        
        # constrain some values to zero (should be only for beta as sigma=0 is undefined)
        beta[is.na(beta)] <- 0
        
      if(!is.null(error)) errorValue <- log(param[nbeta+nsigma+1]) else errorValue <- NULL
      
      # test=sigma+(beta*maxN)
      # if(any(test<=0)){
          #	LL<-list()
        #	LL$logl<-Inf
        #	}else{

      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      if(return.tree){ return(phylo)}
      if(any(is.na(phylo$edge.length))) return(1000000) # instead of checking the parameter values as done previously, I return a high-likelihood value when there are NAs in the branch lengths. Note also that returning Inf value doesn't work with L-BFGS-B algorithm
	  if(any(phylo$edge.length<0)) return(1000000)
      LL<-mvLL(phylo,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))
      
    }
    if(is.na(LL$logl) | is.infinite(LL$logl)){return(1000000)}
    if(results==FALSE){
      return(-LL$logl)
    }else{
      return(list(LL=-LL$logl, mu=LL$theta, s2=sigma))
    }
  }
  
  ##------------------------------------Optimization-------------------------------##
  phyloTrans=NULL
  if(method=="L-BFGS-B" | method=="Nelder-Mead"){
    estim<-optim(par=startval,fn=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)},control=control, hessian=TRUE, method=method, lower=lower, upper=upper)
    
  }else if(method=="fixed"){
    estim <- list()
    estim$par <- param <- c(beta,log(sigma))
    #estim$par <- c(beta, sigma)
    estim$value <- clikCLIM(param=estim$par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)
    estim$convergence <- 0
    
    ## Return the tree --- just some modifications to the previous code to allow retrieving the tree
    # create vector of parameters
    beta <- numeric(nstates)
    sigma <- numeric(nstates)
    # assign values
    beta[] <- c(param[seq_len(nbeta)])[index.user.beta]
    sigma[] <- c(exp(param[nbeta+seq_len(nsigma)]))[index.user.sigma]
    # constrain some values to zero (should be only for beta as sigma=0 is undefined)
    beta[is.na(beta)] <- 0
    phyloTrans <- BranchtransformMAPS(phy, beta, mtot, times, fun, sigma, model, errorValue=NULL)
  }
  
  # Results
  # Prepar the tables for the results
  beta <- numeric(nstates)
  sigma <- numeric(nstates)
        
  if(model=="exponential"){
      
     # assign values
     beta[] <- c(estim$par[seq_len(nbeta)])[index.user.beta]
     sigma[] <- c(exp(estim$par[nbeta+seq_len(nsigma)]))[index.user.sigma]
     # constrain some values to zero (should be only for beta as sigma=0 is undefined)
     beta[is.na(beta)] <- 0
      
    resultList<-matrix(c(beta, sigma), ncol=nstates, byrow=T)
    colnames(resultList)<-c(colnames(tree$mapped.edge))
    rownames(resultList)<-c("beta","sigma")
      
  }else{
      
     # assign values
     beta[] <- c(estim$par[seq_len(nbeta)])[index.user.beta]
     sigma[] <- c(estim$par[nbeta+seq_len(nsigma)])[index.user.sigma]
     
     # constrain some values to zero (should be only for beta as sigma=0 is undefined)
     beta[is.na(beta)] <- 0
      
    resultList<-matrix(c(beta, exp(sigma)), ncol=nstates, byrow=T)
    colnames(resultList)<-c(colnames(tree$mapped.edge))
    rownames(resultList)<-c("beta","sigma")
  }
  
  if(!is.null(error)) errorValue <- estim$par[nsigma+nbeta+1]^2 else errorValue <- NULL
  
  # LogLikelihood
  LL<--estim$value
  # parameter (anc + sigma)
  if(model=="exponential"){
    nparam=1+length(estim$par)
  }else{
    nparam=1+length(estim$par) #sigma estimated in optimization
  }
  
  # AIC
  AIC<--2*LL+2*nparam
  # AIC corrected
  AICc<-AIC+((2*nparam*(nparam+1))/(Ntip(phy)-nparam-1)) #Hurvich et Tsai, 1989
  
  #ancestral states estimates
  anc<-clikCLIM(param=estim$par, dat=data, phy, mtot, times, fun=fun, model, results=TRUE)$mu
  
  if(return.tree){
  transformed.phylo.ML<-clikCLIM(param=estim$par, dat=data, phy, mtot, times, fun=fun, model, results=TRUE,return.tree=TRUE)
  }
  
  ##---------------------Diagnostics--------------------------------------------##
  
  if(estim$convergence==0 & diagnostic==TRUE){
    cat("\n","successful convergence of the optimizer","\n")
  }else if(estim$convergence==1 & diagnostic==TRUE){
    cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
  }else if(diagnostic==TRUE){
    cat("\n","convergence of the optimizer has not been reached, try simpler model","\n")
  }
  
  # Hessian eigen decomposition to check the derivatives
  if(method=="L-BFGS-B" | method=="Nelder-Mead"){
    hess<-eigen(estim$hessian)$values
  }else{
    hess<-0
  }
  if(any(hess<0)){
    hess.value<-1
    if(diagnostic==TRUE){
      cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
  }else{
    hess.value<-0
    if(diagnostic==TRUE){
      cat("a reliable solution has been reached","\n")}
  }
  
  ##-------------------Print results--------------------------------------------##
  if(echo==TRUE){
    cat("\n")
    cat("Summary results for the",model," model","\n")
    cat("LogLikelihood:","\t",LL,"\n")
    cat("AIC:","\t",AIC,"\n")
    cat("AICc:","\t",AICc,"\n")
    cat(nparam,"parameters")
    cat("\n")
    cat("Estimated rates matrix","\n")
    print(resultList)
    cat("\n")
    cat("Estimated ancestral state","\n")
    cat(anc)
    cat("\n")
    if(!is.null(error)){
      cat("\n")
      cat("Estimated error","\n")
      cat(errorValue)
      cat("\n") 
    }
  }
  
  if(return.tree){
  results<-list(LogLik=LL, AIC=AIC, AICc=AICc, rates=resultList, anc=anc, convergence=estim$convergence, hess.values=hess.value, error=errorValue, param=estim$par, phyloTrans=transformed.phylo.ML)
  } else{
  results<-list(LogLik=LL, AIC=AIC, AICc=AICc, rates=resultList, anc=anc, convergence=estim$convergence, hess.values=hess.value, error=errorValue, param=estim$par, phyloTrans=phyloTrans)
	}
}

.fit_t_DD<-function(phylo,data,model=c("exponential","linear","both"),geo.map=NULL,subgroup.map=NULL,subgroup=NULL,regime.map=NULL,error=NULL, beta=NULL, sigma=NULL, method=c("Nelder-Mead","L-BFGS-B","BB"), upper=Inf, lower=-Inf, control=list(maxit=20000), diagnostic=FALSE, echo=FALSE){
	
	if(!model%in%c("exponential","linear","both")){ stop("model must be stated as 'exponential' , 'linear', or 'both' ")}
	
if(is.null(geo.map)&&is.null(subgroup.map)&&is.null(regime.map)){ 	# single slope version without BioGeoBEARS biogeography or subgroup pruning
	
		hold<-rep("A",length(phylo$tip.label))
		hold[1]<-"B"
		names(hold)<-phylo$tip.label
		smap<-make.simmap(phylo,hold,message=F)
		new.maps<-list()
		for(i in 1:length(phylo$edge.length)){
			new.maps[[i]]<-phylo$edge.length[i]
			names(new.maps[[i]])<-"A"
			}
		new.mapped.edge<- as.matrix(rowSums(smap$mapped.edge))
		colnames(new.mapped.edge)<-"A"	
		smap$maps<-new.maps
		smap$mapped.edge<-new.mapped.edge
	
		times = as.numeric(sort(max(branching.times(smap))-branching.times(smap)))
		
		class.df<-.return.class.df_sympatric(smap)
		new_list_function<-create.function.list(smap,times=times,df=class.df)

		sigma.constraint<-rep(1, dim(smap$mapped.edge)[2])
		beta.constraint<-rep(1, dim(smap$mapped.edge)[2])

		if(model%in%c("exponential","linear")){
			out<-fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=class.df,input.times=times,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	
		} else {
			out.exp<-fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=class.df,input.times=times,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	
			out.lin<-fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=class.df,input.times=times,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	
			out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}

		
}  else if (is.null(geo.map)&&is.null(subgroup.map)&&!is.null(regime.map)) { # two slope version without BioGeoBEARS biogeography or subgroup pruning

		if(!all(as.phylo(phylo)$tip.label == as.phylo(regime.map)$tip.label)) { stop("regime map doesn't match phylogeny")}
		if(length(data) != length(as.phylo(regime.map)$tip.label)) { stop("number of lineages in data and regime map don't match")}
		if(! all (names(data) %in% as.phylo(regime.map)$tip.label)) { stop("names of lineages in data and regime map don't match")}
		if(! all (as.phylo(regime.map)$tip.label %in% names(data)) ) { stop("names of lineages in data and regime map don't match")}
		
		class.object<-try(CreateClassObject(regime.map))
		if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(regime.map,rnd=6))}
		if(class(class.object)=="try-error"){class.object<-CreateClassObject(regime.map,rnd=7)}

		class.df<-.return.class.df_subgroup(regime.map,class.object)
		new_list_function<-create.function.list(regime.map,times=class.object$times,df=class.df)
				
		sigma.constraint<-rep(1, dim(regime.map$mapped.edge)[2])
		beta.constraint<-seq(1,by=1,length.out=dim(regime.map$mapped.edge)[2])
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp=fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin=fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
		
}  else if (is.null(subgroup.map)&&is.null(regime.map)&&!is.null(geo.map)) { # single slope version with BioGeoBEARS biogeography but no subgroup pruning

		geo.simmap<-geo.map
		hold<-CreateClassObject(geo.simmap)
		geo.class.df<-.return.class.df(geo.simmap,hold)
		class.object=hold
		class.df=geo.class.df
		new_list_function <- create.function.list(geo.simmap=geo.simmap, df=class.df,times=class.object$times)

		sigma.constraint<-rep(1, dim(geo.map$mapped.edge)[2])
		beta.constraint<-rep(1, dim(geo.map$mapped.edge)[2])
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
		
}  else if (is.null(subgroup.map)&&!is.null(regime.map)&&!is.null(geo.map)) {  # two slope version with BioGeoBEARS biogeography but no subgroup pruning
	
		#need to add other cases down the line
		stop("two-slope with biogeography not currently implemented")
	
}  else if (is.null(geo.map)&&is.null(regime.map)&&!is.null(subgroup.map)) { # single slope version with subgroup pruning but no BioGeoBEARS biogeography
				
		if(!all(phylo$tip.label %in% as.phylo(subgroup.map)$tip.label)) { stop("some lineages in phylogeny don't appear in subgroup map")}
		if( is.null(subgroup) || (!subgroup%in%colnames(subgroup.map$mapped.edge))){ stop("specify a subgroup that appears as a mapped regime in subgroup.map")}
		trimclass.subgroup.trimmed<-.trimSimmap(subgroup.map,trim.class=subgroup)
		
		if(length(colnames(trimclass.subgroup.trimmed$mapped.edge))==1){ return("FLAG")}
		
		class.object<-try(CreateClassObject(trimclass.subgroup.trimmed))
		
		if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(trimclass.subgroup.trimmed,rnd=6))}
		if(class(class.object)=="try-error"){class.object<-CreateClassObject(trimclass.subgroup.trimmed,rnd=7)}

		subgroup.class.df<-.return.class.df_subgroup(trimclass.subgroup.trimmed,class.object)
		
		subgroup.class.df[,which(colnames(trimclass.subgroup.trimmed$mapped.edge)!=subgroup)+1]=1 
	
		
		trimclass.subgroup.trimmed.tips<-drop.tip.simmap(trimclass.subgroup.trimmed,trimclass.subgroup.trimmed$tip.label[which(!trimclass.subgroup.trimmed$tip.label%in%names(data))])
		subgroup.class.df.trimmed<-subgroup.class.df[,c(1,match(colnames(trimclass.subgroup.trimmed.tips$mapped.edge),colnames(trimclass.subgroup.trimmed$mapped.edge))+1)]		
	
		subgroup.map.region.root=max(nodeHeights(trimclass.subgroup.trimmed))
		trimclass.subgroup.trimmed.tips.root=max(nodeHeights(trimclass.subgroup.trimmed.tips))

		if(round(subgroup.map.region.root,5)!=round(trimclass.subgroup.trimmed.tips.root,5)){

			trimmed.class.object<-try(CreateClassObject(trimclass.subgroup.trimmed.tips))
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-try(CreateClassObject(trimclass.subgroup.trimmed.tips,rnd=6))}
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-CreateClassObject(trimclass.subgroup.trimmed.tips,rnd=7)}

			shifted.times<-trimmed.class.object$times+(subgroup.map.region.root-trimclass.subgroup.trimmed.tips.root)
			new.subgroup.class.df.trimmed<-subgroup.class.df.trimmed[c(which(round(class.object$times,5)==round(min(shifted.times),5)):dim(subgroup.class.df.trimmed)[1]),]
			
			new.subgroup.class.df.trimmed$interval<-c(1:dim(new.subgroup.class.df.trimmed)[1])
			class.object$times<-class.object$times[which(round(class.object$times,5)>=round(min(shifted.times),5))]-round(subgroup.map.region.root-trimclass.subgroup.trimmed.tips.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
			subgroup.class.df.trimmed<-new.subgroup.class.df.trimmed
		}
	
		class.df=subgroup.class.df.trimmed
		new_list_function<-create.function.list(trimclass.subgroup.trimmed.tips,times=class.object$times,df=class.df)
				
		sigma.constraint<-rep(1, dim(trimclass.subgroup.trimmed.tips$mapped.edge)[2])
		beta.constraint<-rep(NA, dim(trimclass.subgroup.trimmed.tips$mapped.edge)[2])
		beta.constraint[which(colnames(trimclass.subgroup.trimmed.tips$mapped.edge)==subgroup)]<-1
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp=fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin=fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
		
}  else if (is.null(geo.map)&&!is.null(regime.map)&&!is.null(subgroup.map)) {  # two slope version with subgroup pruning but no BioGeoBEARS biogeography

		if(!all(phylo$tip.label %in% as.phylo(subgroup.map)$tip.label)) { stop("some lineages in phylogeny don't appear in subgroup map")}
		if( is.null(subgroup) || (!subgroup%in%colnames(subgroup.map$mapped.edge))){ stop("specify a subgroup that appears as a mapped regime in subgroup.map")}

		class.by.class.object<-try(.CreateClassbyClassObject_mvMORPH(map.guild=subgroup.map,map.regime=regime.map,trim.class=subgroup))
		if(class(class.by.class.object)=="try-error"){class.by.class.object<-try(.CreateClassbyClassObject_mvMORPH(map.guild=subgroup.map,map.regime=regime.map,trim.class=subgroup,rnd=6))}
		if(class(class.by.class.object)=="try-error"){class.by.class.object<-.CreateClassbyClassObject_mvMORPH(map.guild=subgroup.map,map.regime=regime.map,trim.class=subgroup,rnd=7)}
		regime.class.df<-.return.class.df_subgroup(class.by.class.object$regime.simmap,class.by.class.object$regime.class.object)
		regime.class.df[,which(colnames(class.by.class.object$regime.simmap$mapped.edge)=='Z')+1]=1

		regime.simmap.region.trimmed<-drop.tip.simmap(class.by.class.object$regime.simmap,class.by.class.object$regime.simmap$tip.label[which(!class.by.class.object$regime.simmap$tip.label%in%names(data))])
		
		regime.class.df.trimmed<-regime.class.df[,c(1,match(colnames(regime.simmap.region.trimmed$mapped.edge),colnames(class.by.class.object$regime.simmap$mapped.edge))+1)]		

		regime.simmap.region.root=max(nodeHeights(class.by.class.object$regime.simmap))
		regime.simmap.region.trimmed.root=max(nodeHeights(regime.simmap.region.trimmed))

		if(round(regime.simmap.region.root,5)!=round(regime.simmap.region.trimmed.root,5)){
		
			trimmed.class.object<-try(CreateClassObject(regime.simmap.region.trimmed,rnd=5))
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-try(CreateClassObject(regime.simmap.region.trimmed,rnd=6))}
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-CreateClassObject(regime.simmap.region.trimmed,rnd=7)}

			shifted.times<-trimmed.class.object$times+(regime.simmap.region.root-regime.simmap.region.trimmed.root)
			new.regime.class.df.trimmed<-regime.class.df.trimmed[c(which(round(class.by.class.object$regime.class.object$times,5)==round(min(shifted.times),5)):dim(regime.class.df.trimmed)[1]),]
			
			new.regime.class.df.trimmed$interval<-c(1:dim(new.regime.class.df.trimmed)[1])
			class.by.class.object$regime.class.object$times<-class.by.class.object$regime.class.object$times[which(round(class.by.class.object$regime.class.object$times,5)>=round(min(shifted.times),5))]-round(regime.simmap.region.root-regime.simmap.region.trimmed.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the class.by.class.object object, which could lead to errors
			regime.class.df.trimmed<-new.regime.class.df.trimmed
		}
		class.df=regime.class.df.trimmed
		class.object=class.by.class.object$regime.class.object
		
		new_list_function<-create.function.list(regime.simmap.region.trimmed,df=class.df,times=class.object$times)
		
		sigma.constraint<-rep(1, dim(regime.simmap.region.trimmed$mapped.edge)[2])
		beta.constraint<-rep(NA, dim(regime.simmap.region.trimmed$mapped.edge)[2])
		beta.constraint[which(colnames(regime.simmap.region.trimmed$mapped.edge)!="Z")]<-1:(dim(regime.simmap.region.trimmed$mapped.edge)[2]-1)
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
}  else if (is.null(regime.map)&&!is.null(subgroup.map)&&!is.null(geo.map)) { # single slope version with subgroup pruning and BioGeoBEARS biogeography
	
		#need to add other cases down the line (check 20180802_testing.R and CreateGeoByClassObject.R for inspiration)
		stop("single slope version with subgroup pruning and BioGeoBEARS biogeography not yet implemented")
	
	}  else if (!is.null(regime.map)&&!is.null(subgroup.map)&&!is.null(geo.map)) {  # two slope version with subgroup pruning and BioGeoBEARS biogeography
	
		#need to add other cases down the line
		stop("two slope version with subgroup pruning and BioGeoBEARS biogeography not yet implemented")
	}
		return(out)
}



.CreateClassbyClassObject_mvMORPH<-function(map.guild,map.regime,trim.class,rnd=5){

trc=trim.class

##trim tree

new.map<-trimSimmap(map.guild,trc)

##create class.object if not provided

class.object.guild<-CreateClassObject(new.map,rnd=rnd)

#trim regime to match guild

map.regime.trimmed<-drop.tip.simmap(map.regime,map.regime$tip.label[which(!map.regime$tip.label%in%new.map$tip.label)])

class.object.regime<-CreateClassObject(map.regime.trimmed,rnd=rnd,return.mat=TRUE)

##first concatenate geo.object and class.object timings 
	regt<-round(class.object.regime$times,rnd)
	guit<-round(class.object.guild$times,rnd)
	nodeDist<-sort(unique(c(regt,guit)))
	nodeDiff<-diff(nodeDist)
	if(any(nodeDiff<= (2*(10^-rnd)))){stop("potential rounding error, two time bins very similar, try changing rnd digits")}
	

#initialize counter for geo.class.object and class.object
u<-0 #guild class object
y<-0 #regime class object


###NEED TO UPDATE THIS so that only class objects are reconciled, not multiplied against one another
hold.guild<-list()
hold.regime<-list()

for(i in 1:length(nodeDiff)){

	if((nodeDist[i]%in%regt) && (nodeDist[i]%in%guit)){ #if timing is the same for both
		u = u+1
		y = y+1
		gui.int<-class.object.guild$class.object[[u]]
		reg.int<-class.object.regime$class.object[[y]]
		reg.int[which(gui.int[,2]!=trc),2]<-'Z'
		hold.guild[[i]]<-gui.int
		hold.regime[[i]]<-reg.int
	}
	if((nodeDist[i]%in%regt) && (!nodeDist[i]%in%guit)){ #this means that geo.object changes but class object doesn't
		y = y+1
		gui.int<-class.object.guild$class.object[[u]]
		reg.int<-class.object.regime$class.object[[y]]
		reg.int[which(gui.int[,2]!=trc),2]<-'Z'
		hold.guild[[i]]<-gui.int
		hold.regime[[i]]<-reg.int
	}
	if((!nodeDist[i]%in%regt) && (nodeDist[i]%in%guit)){ #this means that class.object changes but geo object doesn't
		u = u+1
		gui.int<-class.object.guild$class.object[[u]]
		reg.int<-class.object.regime$class.object[[y]]
		reg.int[which(gui.int[,2]!=trc),2]<-'Z'
		hold.guild[[i]]<-gui.int
		hold.regime[[i]]<-reg.int
	}
	}


	phylo<-class.object.regime$phylo
	mat<-class.object.regime$mat
	maps.list=list()
	
	for(k in 1:length(phylo$edge.length)){
		
		#identify branch from edge matrix
		#lookup the name of this branch in the 'mat' matrix compiled above
		#lookup which nat elements have the name of this branch
		#write a vector of the nodeDiff values named with the ranges for each of these elements

		
		lf<-phylo$edge[k,1]
		ri<-phylo$edge[k,2]
		br<-mat[which(mat[,1]==lf & mat[,3]==ri),2]
		natis<-which(sapply(hold.regime,function(x)br%in%x[,1]))
		out.vec<-nodeDiff[natis]
		name.vec<-vector()
		for(n in 1:length(out.vec)){
			name.vec<-c(name.vec,hold.regime[[natis[n]]][which(hold.regime[[natis[n]]][,1]==br),2])	
		}	
		names(out.vec)<-name.vec
		
		#sum adjacent elements with the same name
		
		out.vec.simple<-vector()
		counter=1
		for(i in 1: length(out.vec)){
			if(i == 1 || i == (counter+1)){
				while((length(out.vec)>counter) && (names(out.vec[i])==names(out.vec[counter+1]))){
					counter=counter+1
					}
				hold<-sum(out.vec[i:counter])
				names(hold)<-names(out.vec[i])
				out.vec.simple<-c(out.vec.simple,hold)	
				}
		}
		
		
		maps.list[[k]]<-out.vec.simple
	
	}

	mapped.edge<-matrix(nrow=dim(phylo$edge)[1],ncol=length(unique(names(unlist(maps.list)))))
	colnames(mapped.edge)<-unique(names(unlist(maps.list)))
	
	for(k in 1:dim(phylo$edge)[1]){
		for(j in 1:dim(mapped.edge)[2]){
			hold<-which(names(maps.list[[k]])==colnames(mapped.edge)[j])
			mapped.edge[k,j]<-ifelse(length(hold)==0,0,sum(maps.list[[k]][hold]))
		}
	}
	
	outsmap<-list(edge=phylo$edge,edge.length=phylo$edge.length,tip.label=phylo$tip.label,Nnode=phylo$Nnode,maps=maps.list,mapped.edge=mapped.edge,Q="NA",logL="NA")
	class(outsmap)<-c("phylo","simmap")

	
return(list(subgroup.simmap=new.map,subgroup.class.object=list(class.object=hold.guild,times=nodeDist,spans=nodeDiff),regime.simmap=outsmap,regime.class.object=list(class.object=hold.regime,times=nodeDist,spans=nodeDiff)))#new phylo object, #new times, #new spans, #new geo object

}

.return.class.df<-function(simmap,class.object){
	states<-colnames(simmap$mapped.edge)
	for(i in 1:length(states)){
		st.id=paste("c(",paste(which(grepl(paste(strsplit(states[i],split="")[[1]],collapse="|"),states)),collapse=","),")",sep="") #this gives the columns to extract from class.df
		#eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]==states[',i,']))',sep="")))
		eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]%in%states[',st.id,']))',sep="")))
	}
	eval(parse(text=paste('return(data.frame(interval=1:length(d1),',paste('d',1:length(states),sep="",collapse=','),'))',sep="")))
}

.return.class.df_sympatric<-function(simmap){
	states<-colnames(simmap$mapped.edge)
	d1<-2:length(simmap$tip.label)
	eval(parse(text=paste('return(data.frame(interval=1:length(d1),',paste('d',1:length(states),sep="",collapse=','),'))',sep="")))
}


.return.class.df_subgroup<-function(simmap,class.object){
	states<-colnames(simmap$mapped.edge)
	for(i in 1:length(states)){
		st.id=i #this gives the columns to extract from class.df
		#eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]==states[',i,']))',sep="")))
		eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]%in%states[',st.id,']))',sep="")))
	}
	eval(parse(text=paste('return(data.frame(interval=1:length(d1),',paste('d',1:length(states),sep="",collapse=','),'))',sep="")))
}