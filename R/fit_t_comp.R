fit_t_comp<-function(phylo,data,model=c("MC","DDexp","DDlin"),par=NULL,geography.object=NULL,method="Nelder-Mead"){

#check to make sure data are univariate, with names matching phylo object
if(length(data)!=length(phylo$tip.label)){stop("length of data does not match length of tree")}
if(!is.null(dim(data))){stop("data needs to be a single trait")}
if(is.null(par)){par<-c(var(data)/max(nodeHeights(phylo)),0)}

if(is.null(geography.object)){
	if(model=="MC"){
		opt<-optim(par,likelihood_t_MC,phylo=phylo,data=data,method=method)
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = abs(opt$par[1]), S = opt$par[2], convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par,likelihood_t_DD,phylo=phylo,data=data,model="DDexp",method=method)
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = abs(opt$par[1]), r = opt$par[2], convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		opt<-optim(par,likelihood_t_DD,phylo=phylo,data=data,model="DDlin",method=method)
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = abs(opt$par[1]), b = opt$par[2], convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object)){
#check to make sure length matches length of nodeDiff
	if(length(geography.object)!=phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
	if(model=="MC"){
		opt<-optim(par,likelihood_t_MC_geog,phylo=phylo,geography.object=geography.object,data=data,method=method)
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = abs(opt$par[1]), S = opt$par[2], convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par,likelihood_t_DD_geog,phylo=phylo,geography.object=geography.object,data=data,model="DDexp",method=method)
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = exp(opt$par[1]), r = opt$par[2], convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		opt<-optim(par,likelihood_t_DD_geog,phylo=phylo,geography.object=geography.object,data=data,model="DDlin",method=method)
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = exp(opt$par[1]), b = opt$par[2], convergence = opt$convergence)
		return(results)
		}

}
}