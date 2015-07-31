fitCompModel<-function(phylo,data,model=c("MC","DDexp","DDlin"),geography.object=NULL,method="Nelder-Mead"){

#check to make sure data are univariate, with names matching phylo object
if(length(data)!=length(phylo$tip.label)){stop("length of data does not match length of tree")}
if(!is.null(dim(data))){stop("data needs to be a single trait")}

if(is.null(geography.object)){
	if(model=="MC"){
		opt<-optim(par=c(0.5,-0.45),likelihood_MC,phylo=phylo,data=data,method=method)
		results<-list(lnL = -opt$value, sig2 = abs(opt$par[1]), S = opt$par[2], aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par=c(var(data)/max(nodeHeights(phylo)),0),likelihood_DD,phylo=phylo,data=data,model="DDexp",method=method)
		results<-list(lnL = -opt$value, sig2 = abs(opt$par[1]), r = opt$par[2], aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		opt<-optim(par=c(var(data)/max(nodeHeights(phylo)),0),likelihood_DD,phylo=phylo,data=data,model="DDlin",method=method)
		results<-list(lnL = -opt$value, sig2 = abs(opt$par[1]), b = opt$par[2], aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object)){
#check to make sure length matches length of nodeDiff
	if(length(geography.object)!=phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
	if(model=="MC"){
		opt<-optim(par=c(0.5,-0.45),likelihood_MC_geog,phylo=phylo,geography.object=geography.object,data=data,method=method)
		results<-list(lnL = -opt$value, sig2 = abs(opt$par[1]), S = opt$par[2], aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par=c(var(data)/max(nodeHeights(phylo)),0),likelihood_DD_geog,phylo=phylo,geography.object=geography.object,data=data,model="DDexp",method=method)
		results<-list(lnL = -opt$value, sig2 = exp(opt$par[1]), r = opt$par[2], aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		opt<-optim(par=c(var(data)/max(nodeHeights(phylo)),0),likelihood_DD_geog,phylo=phylo,geography.object=geography.object,data=data,model="DDlin",method=method)
		results<-list(lnL = -opt$value, sig2 = exp(opt$par[1]), b = opt$par[2], aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, convergence = opt$convergence)
		return(results)
		}

}
}