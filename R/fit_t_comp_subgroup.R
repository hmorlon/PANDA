fit_t_comp_subgroup<-function(full.phylo,ana.events,clado.events,stratified=FALSE,map,data,trim.class,model=c("MC","DDexp","DDlin"),error=NULL,par=NULL,method="Nelder-Mead",bounds=NULL){

	if(is.null(names(data))){stop("data missing taxa names")}
	if(!is.null(dim(data))){stop("data needs to be a single trait")}
	is_tip <- full.phylo$edge[,2] <= length(full.phylo$tip.label)
	if(sum(diff(full.phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp_subgroup cannot be used with ladderized full.phylogenies')}
	if(length(unique(ape::branching.times(full.phylo)))<length(ape::branching.times(full.phylo))){stop("fit_t_comp requires phylogenies where no nodes occur at precisely the same time [see ape::branching.times(phylo)]")}

	if(is.null(bounds[["lower"]]) & is.null(bounds[["upper"]])){
        bounds$lower = -Inf
        bounds$upper = Inf
    }
    
	GeoByClassObject<-CreateGeobyClassObject(full.phylo,map,trim.class,ana.events,clado.events,stratified=stratified)

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
