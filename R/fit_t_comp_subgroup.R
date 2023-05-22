fit_t_comp_subgroup<-function(full.phylo,data,subgroup,subgroup.map,model=c("MC","DDexp","DDlin"),ana.events=NULL,clado.events=NULL,stratified=FALSE,regime.map=NULL,error=NULL,par=NULL,method="Nelder-Mead",bounds=NULL){

	if(is.null(names(data))){stop("data missing taxa names")}
	if(!is.null(dim(data))){stop("data needs to be a single trait")}
	is_tip <- full.phylo$edge[,2] <= length(full.phylo$tip.label)
	if(sum(diff(full.phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp_subgroup cannot be used with ladderized full.phylogenies')}
	if(length(unique(ape::branching.times(full.phylo)))<length(ape::branching.times(full.phylo))){stop("fit_t_comp requires phylogenies where no nodes occur at precisely the same time [see ape::branching.times(phylo)]")}

	if(is.null(bounds[["lower"]]) & is.null(bounds[["upper"]])){
        bounds$lower = -Inf
        bounds$upper = Inf
    }
    
    if(model=="MC" && (is.null(ana.events) || is.null (clado.events))){stop("MC model without biogeography is currently not implemented, please supply ana.events and clado.events")}
    if(model=="MC" && !is.null(regime.map)){stop("two-regime version of MC model with subgroup pruning is currently not implemented")}
    
    if((is.null(ana.events) && !is.null(clado.events))||(!is.null(ana.events) && is.null(clado.events))) {stop("please provide both ana.events and clado.events when fitting models with biogeography (see ?fit_t_comp_subgroup)")}
    if(!is.null(ana.events) && !is.null(regime.map)){ stop("two-regime models with biogeography and subgroup pruning currently not implemented")}
    if(is.null(ana.events)){ #subgroup pruning for DD models without biogeography
    
    	if(is.null(regime.map)){ #single regime subgroup DD models without biogeography
			if(model=="DDexp"){ 
			   	if(is.null(par)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=par[1]
					beta=par[2]
				}

				opt<-.fit_t_DD(phylo=full.phylo,data=data, error= error,model="exponential",par=par,subgroup=subgroup,subgroup.map=subgroup.map,method=method,bounds=bounds)
    			
    			sig2<-opt$rates["sigma",]
				r<-opt$rates["beta",]
				z0<-opt$anc
				if(!is.null(opt$error)){ 
					mserr = opt$error
					npar=4
					}else{
					mserr = NA
					npar=3
					}
				results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = npar, sig2 = as.numeric(sig2), r = as.numeric(r), z0 = as.numeric(z0), convergence = opt$convergence, nuisance=mserr)
				return(results)
			} 
			
			if(model=="DDlin"){
				if(is.null(par)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=par[1]
					beta=par[2]
				}

				
				opt<-.fit_t_DD(phylo=full.phylo,data=data, error= error,model="linear",par=par,subgroup=subgroup,subgroup.map=subgroup.map,method=method,bounds=bounds)
    			
    			sig2<-opt$rates["sigma",subgroup]
				b<-opt$rates["beta",subgroup]
				z0<-opt$anc
				if(!is.null(opt$error)){ 
					mserr = opt$error
					npar=4
					}else{
					mserr = NA
					npar=3
					}
				results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = npar, sig2 = as.numeric(sig2), b = as.numeric(b), z0 = as.numeric(z0), convergence = opt$convergence, nuisance=mserr)
				return(results)
			
			
			}
			
			
		} else { #two-regime models
		
			if(model=="DDexp"){ 
			    if(is.null(par)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=par[1]
					beta=par[2]
				}

				opt<-.fit_t_DD(phylo=full.phylo,data=data, error= error,model="exponential",par=par,subgroup=subgroup,subgroup.map=subgroup.map,regime.map=regime.map,method=method,bounds=bounds)
     			sig2<-opt$rates["sigma",1]
				r1<-opt$rates["beta",1]
				r2<-opt$rates["beta",2]
				z0<-opt$anc
				if(!is.null(opt$error)){ 
					mserr = opt$error
					npar=5
					}else{
					mserr = NA
					npar=4
					}
					
				eval(pare(text=paste0("results<-list(LH = ",opt$LogLik,", aic = ",opt$AIC,", aicc = ",opt$AICc,", free.parameters =",npar,", sig2 = ",as.numeric(sig2),", r1_",colnames(opt$rates)[1]," = ",as.numeric(r1),", r2_",colnames(opt$rates)[2]," = ",as.numeric(r2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,", nuisance = ",mserr,")")))
				return(results)
			} 
			
			if(model=="DDlin"){
			    if(is.null(par)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=par[1]
					beta=par[2]
				}

				opt<-.fit_t_DD(phylo=full.phylo,data=data, error= error,model="linear",par=par,subgroup=subgroup,subgroup.map=subgroup.map,regime.map=regime.map,method=method,bounds=bounds)
    			sig2<-opt$rates["sigma",1]
				b1<-opt$rates["beta",1]
				b2<-opt$rates["beta",2]
				z0<-opt$anc
				if(!is.null(opt$error)){ 
					mserr = opt$error
					npar=5
					}else{
					mserr = NA
					npar=4
					}
					
				eval(pare(text=paste0("results<-list(LH = ",opt$LogLik,", aic = ",opt$AIC,", aicc = ",opt$AICc,", free.parameters =",npar,", sig2 = ",as.numeric(sig2),", b1_",colnames(opt$rates)[1]," = ",as.numeric(b1),", b2_",colnames(opt$rates)[2]," = ",as.numeric(b2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,", nuisance = ",mserr,")")))
				return(results)
			}

		
		
		}	
    	
    }

	GeoByClassObject<-CreateGeobyClassObject(full.phylo,subgroup.map,trim.class=subgroup,ana.events,clado.events,stratified=stratified)

	phylo<-GeoByClassObject$map
	
	root.trimmed.phylo<-max(nodeHeights(phylo))
	root.data<-max(nodeHeights(drop.tip.simmap(phylo,phylo$tip.label[which(!phylo$tip.label%in%names(data))])))
	
	if(round(root.trimmed.phylo,5)!=round(root.data,5)){stop("error where root of trimmed simmap and root of target clade don't match")}
	#if this above error is triggered, need to use scripts written by JPD using JC mvMORPH approach in summer 2018
	
	if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
	geo.object<-GeoByClassObject$geo.object
	#geo.sorted<-.resortGeoObject(phylo,geo.object) 

	if(length(geo.object$geography.object)<phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
	if(length(data)>length(phylo$tip.label)){stop("error: some tips missing from pruned simmap")}
	if(!all(names(data)%in%phylo$tip.label)){stop("error: some tips missing from pruned simmap")}

	if(is.null(par)){par<-c(log(sqrt(var(data)/max(nodeHeights(extract.clade(phylo,getMRCA(phylo,names(data))))))),0)}
    
    # number of parameters estimated by models
    npar = 3
    
    # if NA is provided to "error", then we can estimate nuisance even if we don't have known measurement errors
    if(!is.null(error)){
        
        if(!any(is.na(error))){
            if(is.null(names(error))){
                stop("You should provide a named vector for \"error\" ")
            }else{
                error<-error[phylo$tip.label]
            }
        
        }else{
            error <- rep(0, Ntip(phylo))
            names(error) = phylo$tip.label
        }
        par <- c(par, 0.05*exp(par[1]))
        npar <- npar+1
    }

	
	if(model=="MC"){
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="MC", error=error, method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])^2
		S = -abs(opt$par[2])
        if(!is.null(error)) mserr = exp(opt$par[3]) else mserr = NA
		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="MC",par=opt$par,return.z0=TRUE, error=error)
		results<-list(model = model, LH = -opt$value, aic = (2*npar - 2*(-opt$value)), aicc = (2*npar - 2*(-opt$value))+((2*npar*(npar+1))/(length(phylo$tip.label)-npar-1)), free.parameters = npar, sig2 = sig2, S = S, z0 = as.numeric(z0), convergence = opt$convergence, nuisance=mserr)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="DDexp", error=error, method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])^2
		r = opt$par[2]
        if(!is.null(error)) mserr = exp(opt$par[3]) else mserr = NA

		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDexp",par=opt$par,return.z0=TRUE, error=error)
		results<-list(model = model, LH = -opt$value, aic = (2*npar - 2*(-opt$value)), aicc = (2*npar - 2*(-opt$value))+((2*npar*(npar+1))/(length(phylo$tip.label)-npar-1)), free.parameters = npar, sig2 = sig2, r = r, z0 = as.numeric(z0), convergence = opt$convergence, nuisance=mserr)
		return(results)
		}
	if(model=="DDlin"){
		geography.matrix<-geo.object$geography.object
		maxN<-max(vapply(geography.matrix,function(x)max(rowSums(x)),1))
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,maxN=maxN,data=data,model="DDlin", error=error, method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])^2
		b = opt$par[2]
        if(!is.null(error)) mserr = exp(opt$par[3]) else mserr = NA

		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDlin",par=opt$par,return.z0=TRUE,maxN=maxN, error=error)
		results<-list(model = model, LH = -opt$value, aic = (2*npar - 2*(-opt$value)), aicc = (2*npar - 2*(-opt$value))+((2*npar*(npar+1))/(length(phylo$tip.label)-npar-1)), free.parameters = npar, sig2 = sig2, b = b, z0 = as.numeric(z0), convergence = opt$convergence, nuisance=mserr)
		return(results)
		}

}
