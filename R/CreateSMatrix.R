.CreateSMatrix<-function(class.obj,S.cats=NULL){
	if(is.null(S.cats)){
	scl<-as.list(unique(class.obj$class.obj[[length(class.obj$class.obj)]][,2]))
	if(length(scl)!=2){stop("more than two competitive regimes specified in class object")}
	S1=scl[[1]]
	S2=scl[[2]]
	} else {
	if(length(S.cats)!=2){stop("S.cats must contain a list of two vectors specifying which states correspond to names of competitive regimes in simmap")}
	S1=S.cats[[1]]
	S2=S.cats[[2]]
	}
	out.list<-list()
	for(i in 1:length(class.obj$class.object)){
		imat<-class.obj$class.object[[i]]
		out.mat<-matrix(nrow=2,ncol=dim(imat)[1])
		colnames(out.mat)<-imat[,1]
		out.mat[1,]<-imat[,2]%in%S1
		out.mat[2,]<-imat[,2]%in%S2
		out.list[[i]]<-out.mat*1
		}
	return(list(S.matrix=out.list,times=class.obj$times,spans=class.obj$spans, S1=S1, S2=S2))
	}	