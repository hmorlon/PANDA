spectR_t<-function(phylo,dat,draw=F){
	if (length(phylo$tip.label)!=length(dat)) 
        stop("dat do not match phylo")
       
    #remove internal nodes from distance matrix    
	phylo$tip.label<-1:length(phylo$tip.label)
	phylomat<-dist.nodes(phylo)
	tipmat<-phylomat[1:length(phylo$tip.label),1:length(phylo$tip.label)]
	difs<-as.matrix(dist(dat))
	left<-phylomat[length(phylo$tip.label):dim(phylomat)[1],length(phylo$tip.label):dim(phylomat)[2]]
	difs*tipmat->difmat
	
	#compute eigenvalues from normalized MGL of tip distance matrix
	x<-eigen(
		graph.laplacian(
			graph.adjacency(difmat,
			weighted=T,mode='undirected'),
		normalized=T),
	only.values=F)
	x<-subset(x$values,x$values>0.1)
	
		l=length(x)
		integr<-function(x,f){
			int=0.5*sum((x[2:l]-x[1:(l-1)])*(f[2:l]+f[1:(l-1)]))			
			return(int)
		}
		
	#estimate summary stats	
		density(x)->d
		d$y/integr(d$x,d$y)->dsc
		splitter=max(x)
		tracer=max(dsc)
		fragmenter=(sum((x-mean(x))^3)/l)/(sum((x-mean(x))^2)/l)^(3/2)
		res<-list(eigenvalues=x,splitter=splitter,tracer=tracer,fragmenter=fragmenter)
		#plot spectral density profile
		if(draw==T){
		par(mar=c(4,5,1,1))
		plot(d$x,dsc,type='l',
			xlab=expression(''[n]*lambda),
			ylab=expression(paste('f(x)/',integral(f(y)*dy)),sep=''))
			polygon(d$x,dsc,col=colors(1)[runif(1,1,500)])	
		return(res)
		}
	else{return(res)}
	}
