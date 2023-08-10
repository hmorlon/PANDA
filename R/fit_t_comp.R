fit_t_comp<-function(phylo, data, error=NULL, model=c("MC","DDexp","DDlin"),pars=NULL,geography.object=NULL, regime.map=NULL){

#check to make sure data are univariate, with names matching phylo object
if(length(data)!=length(phylo$tip.label)){stop("length of data does not match length of tree")}
if(is.null(names(data))){stop("data missing taxa names")}
if(!is.null(dim(data))){stop("data needs to be a single trait")}
is_tip <- phylo$edge[,2] <= length(phylo$tip.label)
if(sum(diff(phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp cannot be used with ladderized phylogenies')}
if(length(unique(ape::branching.times(phylo)))<length(ape::branching.times(phylo))){stop("fit_t_comp requires phylogenies where no nodes occur at precisely the same time [see ape::branching.times(phylo)]")}

if(!is.null(regime.map) && (length(regime.map$tip.label)!=length(phylo$tip.label))){stop("regime map does not match phylo")} 

if(is.null(error)){
	
	if(is.null(geography.object) & is.null(regime.map)){ #single-slope version for sympatric clades
		if(model=="MC"){
			if(is.null(pars)){pars<-c(log(sqrt(var(data)/max(nodeHeights(phylo)))),0)}
			params0<-c(0,pars)
			mc.ob<-.createModel_MC(phylo)
			opt<-fitTipData(mc.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			S<-opt$inferredParams[3]
			z0<-opt$inferredParams[1]
			results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), S = as.numeric(S), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
			}
		if(model=="DDexp"){
			if(is.null(pars)){
				sigma=NULL
				beta=NULL
			} else {
				sigma=pars[1]
				beta=pars[2]
			}
			opt<-.fit_t_DD(phylo,data,model="exponential",sigma=sigma,beta=beta)
			sig2<-opt$rates["sigma",]
			r<-opt$rates["beta",]
			z0<-opt$anc
			results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 3, sig2 = as.numeric(sig2), r = as.numeric(r), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
			} 
			
		if(model=="DDlin"){
			if(is.null(pars)){
				sigma=NULL
				beta=NULL
			} else {
				sigma=pars[1]
				beta=pars[2]
			}
			opt<-.fit_t_DD(phylo,data,model="linear",sigma=sigma,beta=beta)
			sig2<-opt$rates["sigma",]
			r<-opt$rates["beta",]
			z0<-opt$anc
			results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 3, sig2 = as.numeric(sig2), r = as.numeric(r), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
			} 
	}
	
	if(!is.null(geography.object) & is.null(regime.map)){ #single-slope version with biogeography
		
		if(class(geography.object)[1]=="list"){
			if(is.null(pars)){pars<-c(log(sqrt(var(data)/max(nodeHeights(phylo)))),0)}
			params0<-c(0,pars)
			if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
			sgeo<-.resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in Marc code
			if(model=="MC"){
				mc.ob<-.createModel_MC_geo(phylo,sgeo)
				opt<-fitTipData(mc.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
				sig2<-(exp(opt$inferredParams[2]))^2
				S<-opt$inferredParams[3]
				z0<-opt$inferredParams[1]
				results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), S = as.numeric(S), z0 = as.numeric(z0), convergence = opt$convergence)
				return(results)
				}
			if(model=="DDexp"){
				ddexp.ob<-.createModel_DDexp_geo(phylo,sgeo)
				opt<-fitTipData(ddexp.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
				sig2<-(exp(opt$inferredParams[2]))^2
				r<-opt$inferredParams[3]
				z0<-opt$inferredParams[1]
				results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), r = as.numeric(r), z0 = as.numeric(z0), convergence = opt$convergence)
				return(results)
				}
			if(model=="DDlin"){
				ddlin.ob<-.createModel_DDlin_geo(phylo,sgeo)
				opt<-fitTipData(ddlin.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
				sig2<-(exp(opt$inferredParams[2]))^2
				b<-opt$inferredParams[3]
				z0<-opt$inferredParams[1]
				results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = as.numeric(sig2), b = as.numeric(b), z0 = as.numeric(z0), convergence = opt$convergence)
				return(results)
				}
		}
		
		if(class(geography.object)[1]=="simmap"){
			if(model=="MC"){stop("MC fits not supported with simmap geo.object (see ?fit_t_comp)")}
			if(model=="DDexp"){
				if(is.null(pars)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=pars[1]
					beta=pars[2]
				}
				opt<-.fit_t_DD(phylo,data,model="exponential",geo.map=geography.object,sigma=sigma,beta=beta)
				sig2<-opt$rates["sigma",]
				r<-opt$rates["beta",]
				z0<-opt$anc
				results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 3, sig2 = as.numeric(sig2)[1], r = as.numeric(r)[1], z0 = as.numeric(z0), convergence = opt$convergence)
				return(results)
				} 
				
			if(model=="DDlin"){
				if(is.null(pars)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=pars[1]
					beta=pars[2]
				}
				opt<-.fit_t_DD(phylo,data,model="linear",geo.map=geography.object,sigma=sigma,beta=beta)
				sig2<-opt$rates["sigma",]
				r<-opt$rates["beta",]
				z0<-opt$anc
				results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 3, sig2 = as.numeric(sig2)[1], r = as.numeric(r)[1], z0 = as.numeric(z0), convergence = opt$convergence)
				return(results)
				} 

		}

	}
	
	if(is.null(geography.object) & !is.null(regime.map)){ #two-slope version for sympatric clades

		if(model=="MC"){
			class.object<-try(CreateClassObject(regime.map))
			if(inherits(class.object, "try-error")){
				class.object<-CreateClassObject(regime.map,rnd=6)
				}
					
			SMatrix<-.CreateSMatrix(class.object)
			smat<-.resortSMatrix(phylo, SMatrix)

			if(is.null(pars)){pars<-c(log(sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1)}
			params0<-c(0,pars)
			mc.ob<-.createModel_MC_twoS(phylo,S.object=smat)
			opt<-fitTipData(mc.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			S1<-opt$inferredParams[3]
			S2<-opt$inferredParams[4]
			z0<-opt$inferredParams[1]
			eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*4 - 2*(-opt$value)),", aicc = ",(2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)),", free.parameters = 4, sig2 = ",as.numeric(sig2),", S1_",SMatrix$S1," = ",as.numeric(S1),", S2_",SMatrix$S2," = ",as.numeric(S2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
			return(results)
			}
		if(model=="DDexp"){
				if(is.null(pars)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=pars[1]
					beta=pars[2]
				}
				opt<-.fit_t_DD(phylo,data,model="exponential",regime.map=regime.map,sigma=sigma,beta=beta)
				sig2<-opt$rates["sigma",1]
				r1<-opt$rates["beta",1]
				r2<-opt$rates["beta",2]
				z0<-opt$anc
				eval(parse(text=paste0("results<-list(LH = ",opt$LogLik,", aic = ",opt$AIC,", aicc = ",opt$AICc,", free.parameters = 4, sig2 = ",as.numeric(sig2),", r1_",colnames(opt$rates)[1]," = ",as.numeric(r1),", r2_",colnames(opt$rates)[2]," = ",as.numeric(r2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
				return(results)
			}
		if(model=="DDlin"){
				if(is.null(pars)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=pars[1]
					beta=pars[2]
				}
				opt<-.fit_t_DD(phylo,data,model="linear",regime.map=regime.map,sigma=sigma,beta=beta)
				sig2<-opt$rates["sigma",1]
				b1<-opt$rates["beta",1]
				b2<-opt$rates["beta",2]
				z0<-opt$anc
				eval(parse(text=paste0("results<-list(LH = ",opt$LogLik,", aic = ",opt$AIC,", aicc = ",opt$AICc,", free.parameters = 4, sig2 = ",as.numeric(sig2),", b1_",colnames(opt$rates)[1]," = ",as.numeric(b1),", b2_",colnames(opt$rates)[2]," = ",as.numeric(b2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
				return(results)
		}
	}
	
	if(!is.null(geography.object) & !is.null(regime.map)){ #two-slope version with biogeography
		if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
		sgeo0<-.resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in code
		
		class.object<-try(CreateClassObject(regime.map))
		if(inherits(class.object, "try-error")){
			class.object<-CreateClassObject(regime.map,rnd=6)
			}
					
		SMatrix<-.CreateSMatrix(class.object)
		smat0<-.resortSMatrix(phylo, SMatrix)
		
		int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo0,S.matrix=smat0))	
		
		#some catches in case there are small rounding issues (happens when events anagenetic in biogeography or regimes happen at a very similar time)				
		if(inherits(int, "try-error")){
			int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=6))
			}	
		if(inherits(int, "try-error")){
			int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=7))
			}	
		if(inherits(int, "try-error")){
			int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=4))
			}	
		
		sgeo<-int$geo.object
		smat<-int$S.matrix

		if(model=="MC"){
			if(is.null(pars)){pars<-c(log(sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1)}
			params0<-c(0,pars)
			mc.ob<-.createModel_MC_twoS_geo(phylo,geo.object=sgeo,S.object=smat)
			opt<-fitTipData(mc.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			S1<-opt$inferredParams[3]
			S2<-opt$inferredParams[4]
			z0<-opt$inferredParams[1]
			eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*4 - 2*(-opt$value)),", aicc = ",(2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)),", free.parameters = 4, sig2 = ",as.numeric(sig2),", S1_",SMatrix$S1," = ",as.numeric(S1),", S2_",SMatrix$S2," = ",as.numeric(S2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
			return(results)
			}
		if(model=="DDexp"){
			if(is.null(pars)){pars<-c(log(sqrt(var(data)/max(nodeHeights(phylo)))),0,0)}
			params0<-c(0,pars)
			ddexp.ob<-.createModel_DDexp_multi_geo(phylo,geo.object=sgeo,r.object=smat)
			opt<-fitTipData(ddexp.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			r1<-opt$inferredParams[3]
			r2<-opt$inferredParams[4]
			z0<-opt$inferredParams[1]
			eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*4 - 2*(-opt$value)),", aicc = ",(2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)),", free.parameters = 4, sig2 = ",as.numeric(sig2),", r1_",SMatrix$S1," = ",as.numeric(r1),", r2_",SMatrix$S2," = ",as.numeric(r2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
			return(results)
			}
		if(model=="DDlin"){
			if(is.null(pars)){pars<-c(log(sqrt(var(data)/max(nodeHeights(phylo)))),0,0)}
			params0<-c(0,pars)
			ddlin.ob<-.createModel_DDlin_multi_geo(phylo,geo.object=sgeo,r.object=smat)
			opt<-fitTipData(ddlin.ob,data,error=NULL,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			b1<-opt$inferredParams[3]
			b2<-opt$inferredParams[4]
			z0<-opt$inferredParams[1]
			eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*4 - 2*(-opt$value)),", aicc = ",(2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)),", free.parameters = 4, sig2 = ",as.numeric(sig2),", b1_",SMatrix$S1," = ",as.numeric(b1),", b2_",SMatrix$S2," = ",as.numeric(b2),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
			return(results)
			}
	
	}
	
	
	
	
	
} else {

  # if NA is provided to "error", then we can estimate nuisance even if we don't have known measurement errors
  if(any(is.na(error))){
    error <- rep(0, Ntip(phylo))
    names(error) = phylo$tip.label
  } 
  
if(is.null(geography.object) & is.null(regime.map)){ #single-slope version for sympatric clades


	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-.createModel_MC_ME(phylo)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*4 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), S = as.numeric(S), nuisance=as.numeric(nuisance),z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
			if(is.null(pars)){
				sigma=NULL
				beta=NULL
			} else {
				sigma=pars[1]
				beta=pars[2]
			}
			opt<-.fit_t_DD(phylo,data,model="exponential",error=error,sigma=sigma,beta=beta)
			sig2<-opt$rates["sigma",]
			r<-opt$rates["beta",]
			z0<-opt$anc
			nuisance<-opt$error
			results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 4, sig2 = as.numeric(sig2), r = as.numeric(r), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
		} 
	if(model=="DDlin"){
			if(is.null(pars)){
				sigma=NULL
				beta=NULL
			} else {
				sigma=pars[1]
				beta=pars[2]
			}
			opt<-.fit_t_DD(phylo,data,model="linear",error=error,sigma=sigma,beta=beta)
			sig2<-opt$rates["sigma",]
			r<-opt$rates["beta",]
			z0<-opt$anc
			nuisance<-opt$error
			results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 4, sig2 = as.numeric(sig2), r = as.numeric(r), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
		} 
}

if(!is.null(geography.object) & is.null(regime.map)){ #single-slope version with biogeography

	if(class(geography.object)[1]=="list"){
		if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
		sgeo<-.resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in code
		if(model=="MC"){
			if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
			params0<-c(0,pars)
			mc.ob<-.createModel_MC_geo_ME(phylo,geo.object=sgeo)
			opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			S<-opt$inferredParams[3]
			nuisance<-exp(opt$inferredParams[4])
			z0<-opt$inferredParams[1]
			results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc =(2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), S = as.numeric(S), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
			}
		if(model=="DDexp"){
			if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
			params0<-c(0,pars)
			ddexp.ob<-.createModel_DDexp_geo_ME(phylo,geo.object=sgeo)
			opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			r<-opt$inferredParams[3]
			nuisance<-exp(opt$inferredParams[4])
			z0<-opt$inferredParams[1]
			results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), r = as.numeric(r), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
			}
		if(model=="DDlin"){
			if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
			params0<-c(0,pars)
			ddlin.ob<-.createModel_DDlin_geo_ME(phylo,geo.object=sgeo)
			opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
			sig2<-(exp(opt$inferredParams[2]))^2
			b<-opt$inferredParams[3]
			nuisance<-exp(opt$inferredParams[4])
			z0<-opt$inferredParams[1]
			results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), b = as.numeric(b), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
			return(results)
			}
		}
	
		if(class(geography.object)[1]=="simmap"){
			if(model=="MC"){stop("MC fits not supported with simmap geo.object (see ?fit_t_comp)")}
			if(model=="DDexp"){
				if(is.null(pars)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=pars[1]
					beta=pars[2]
				}
				opt<-.fit_t_DD(phylo,data,model="exponential",geo.map=geography.object,error=error,sigma=sigma,beta=beta)
				sig2<-opt$rates["sigma",]
				r<-opt$rates["beta",]
				z0<-opt$anc
				nuisance<-opt$error
				results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 4, sig2 = as.numeric(sig2)[1], r = as.numeric(r)[1], nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
				return(results)
				} 
				
			if(model=="DDlin"){
				if(is.null(pars)){
					sigma=NULL
					beta=NULL
				} else {
					sigma=pars[1]
					beta=pars[2]
				}
				opt<-.fit_t_DD(phylo,data,model="linear",geo.map=geography.object,error=error,sigma=sigma,beta=beta)
				sig2<-opt$rates["sigma",]
				r<-opt$rates["beta",]
				z0<-opt$anc
				nuisance<-opt$error
				results<-list(LH = opt$LogLik, aic = opt$AIC, aicc = opt$AICc, free.parameters = 4, sig2 = as.numeric(sig2)[1], r = as.numeric(r)[1], nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
				return(results)
				} 

		}
	
		
		
}

if(is.null(geography.object) & !is.null(regime.map)){ #multi-slope version for sympatric clades (i.e., no biogeography)
	
	
	if(model=="MC"){
		class.object<-try(CreateClassObject(regime.map))
		if(inherits(class.object, "try-error")){
			class.object<-CreateClassObject(regime.map,rnd=6)
			}
					
		SMatrix<-.CreateSMatrix(class.object)
		smat<-.resortSMatrix(phylo, SMatrix)

		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-.createModel_MC_twoS_ME(phylo,S.object=smat)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S1<-opt$inferredParams[3]
		S2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*5 - 2*(-opt$value)),", aicc = ",(2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)),", free.parameters = 5, sig2 = ",as.numeric(sig2),", S1_",SMatrix$S1," = ",as.numeric(S1),", S2_",SMatrix$S2," = ",as.numeric(S2),", nuisance = ",as.numeric(nuisance),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
		return(results)
		}
		
	if(model=="DDexp"){
			if(is.null(pars)){
				sigma=NULL
				beta=NULL
			} else {
				sigma=pars[1]
				beta=pars[2]
			}
			opt<-.fit_t_DD(phylo,data,error=error,model="exponential",regime.map=regime.map,sigma=sigma,beta=beta)
			sig2<-opt$rates["sigma",1]
			r1<-opt$rates["beta",1]
			r2<-opt$rates["beta",2]
			z0<-opt$anc
			nuisance<-opt$error

			eval(parse(text=paste0("results<-list(LH = ",opt$LogLik,", aic = ",opt$AIC,", aicc = ",opt$AICc,", free.parameters = 5, sig2 = ",as.numeric(sig2),", r1_",colnames(opt$rates)[1]," = ",as.numeric(r1),", r2_",colnames(opt$rates)[2]," = ",as.numeric(r2),", nuisance = ", as.numeric(nuisance), ", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
			return(results)
		}
	if(model=="DDlin"){
			if(is.null(pars)){
				sigma=NULL
				beta=NULL
			} else {
				sigma=pars[1]
				beta=pars[2]
			}
			opt<-.fit_t_DD(phylo,data,error=error,model="linear",regime.map=regime.map,sigma=sigma,beta=beta)
			sig2<-opt$rates["sigma",1]
			b1<-opt$rates["beta",1]
			b2<-opt$rates["beta",2]
			z0<-opt$anc
			nuisance<-opt$error
			eval(parse(text=paste0("results<-list(LH = ",opt$LogLik,", aic = ",opt$AIC,", aicc = ",opt$AICc,", free.parameters = 5, sig2 = ",as.numeric(sig2),", b1_",colnames(opt$rates)[1]," = ",as.numeric(b1),", b2_",colnames(opt$rates)[2]," = ",as.numeric(b2),", nuisance = ", as.numeric(nuisance), ", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
			return(results)
	}
		
		
		
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-.createModel_DDexp_multi_ME(phylo,r.object=smat)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r1<-opt$inferredParams[3]
		r2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*5 - 2*(-opt$value)),", aicc = ",(2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)),", free.parameters = 5, sig2 = ",as.numeric(sig2),", r1_",SMatrix$S1," = ",as.numeric(r1),", r2_",SMatrix$S2," = ",as.numeric(r2),", nuisance = ",as.numeric(nuisance),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-.createModel_DDlin_multi_ME(phylo,r.object=smat)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b1<-opt$inferredParams[3]
		b2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*5 - 2*(-opt$value)),", aicc = ",(2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)),", free.parameters = 5, sig2 = ",as.numeric(sig2),", b1_",SMatrix$S1," = ",as.numeric(b1),", b2_",SMatrix$S2," = ",as.numeric(b2),", nuisance = ",as.numeric(nuisance),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
		return(results)
		}
}

if(!is.null(geography.object) & !is.null(regime.map)){ #multi-slope version with biogeography

	if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
	sgeo0<-.resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in code
	
	class.object<-try(CreateClassObject(regime.map))
	if(inherits(class.object, "try-error")){
		class.object<-CreateClassObject(regime.map,rnd=6)
		}
				
	SMatrix<-.CreateSMatrix(class.object)
	smat0<-.resortSMatrix(phylo, SMatrix)
	
	int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo0,S.matrix=smat0))	
	
	#some catches in case there are small rounding issues (happens when events anagenetic in biogeography or regimes happen at a very similar time)				
	if(inherits(int, "try-error")){
		int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=6))
		}	
	if(inherits(int, "try-error")){
		int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=7))
		}	
	if(inherits(int, "try-error")){
		int<-try(.ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=4))
		}	
	
	sgeo<-int$geo.object
	smat<-int$S.matrix
	
	
	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-.createModel_MC_twoS_geo_ME(phylo,geo.object=sgeo,S.object=smat)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S1<-opt$inferredParams[3]
		S2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*5 - 2*(-opt$value)),", aicc = ",(2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)),", free.parameters = 5, sig2 = ",as.numeric(sig2),", S1_",SMatrix$S1," = ",as.numeric(S1),", S2_",SMatrix$S2," = ",as.numeric(S2),", nuisance = ",as.numeric(nuisance),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-.createModel_DDexp_multi_geo_ME(phylo,geo.object=sgeo,r.object=smat)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r1<-opt$inferredParams[3]
		r2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*5 - 2*(-opt$value)),", aicc = ",(2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)),", free.parameters = 5, sig2 = ",as.numeric(sig2),", r1_",SMatrix$S1," = ",as.numeric(r1),", r2_",SMatrix$S2," = ",as.numeric(r2),", nuisance = ",as.numeric(nuisance),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-.createModel_DDlin_multi_geo_ME(phylo,geo.object=sgeo,r.object=smat)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b1<-opt$inferredParams[3]
		b2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		eval(parse(text=paste0("results<-list(LH = ",-opt$value,", aic = ",(2*5 - 2*(-opt$value)),", aicc = ",(2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)),", free.parameters = 5, sig2 = ",as.numeric(sig2),", b1_",SMatrix$S1," = ",as.numeric(b1),", b2_",SMatrix$S2," = ",as.numeric(b2),", nuisance = ",as.numeric(nuisance),", z0 = ",as.numeric(z0),", convergence = ",opt$convergence,")")))
		return(results)
		}

}



}

}
