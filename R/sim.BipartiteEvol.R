sim.BipartiteEvol=function(nx,ny=nx,NG,dSpace_H=Inf,dSpace_P=Inf,D=3,muP,muH,alphaP=0,alphaH=0,iniP=0,iniH=0,nP=1,nH=1,rP=1,rH=1,effect=1,
                           verbose=100,thin=1,P=NULL,H=NULL){

  N=nx*ny     # number of individuals
  timeStep=floor(sqrt(NG/(nH+nP)))
  oneByOne=FALSE

  parentsP=1:N
  parentsH=1:N
  pmut_act=rep(0,N)
  hmut_act=rep(0,N)
  changeP=c()
  changeH=c()
  
  # kernel=sapply(1:nx,function(x){sapply(1:ny,function(y){exp(-(1/2)*((min(abs(1-x),abs(nx+1-x))^2+min(abs(1-y),abs(ny+1-y))^2)/(dSpace^2)))})})
  
  Ys=rep(1:ny,each=nx,times=1)
  Xs=rep(1:nx,times=ny)
  # choice of the fitness function
  f=function(X,x,r,alpha,pos,dSpace){
    .Call("fitnessFunction", X=as.numeric(t(X)), x=as.numeric(x),r=as.numeric(r),alpha=as.numeric(alpha), Ncol=as.integer(N), D=as.integer(D),
          dSpace=as.numeric(dSpace), Xs=as.numeric(Xs),Ys=as.numeric(Ys),I=as.integer(pos),PACKAGE = "RPANDA")
  }

  
  # trait initial values
  if(is.null(P)) P=matrix(iniP,nrow=D,ncol=N)
  if(is.null(H)) H=matrix(iniH,nrow=D,ncol=N)
  
  # === compressed outputs ===
  init_cap <- 5
  
  out_xP <- lapply(1:D, function(j) list(
    t = grow_init("integer", init_cap), #integer(0),
    ind  = grow_init("integer", init_cap), #integer(0),
    val  = grow_init("numeric", init_cap) #numeric(0)
  ))
  
  out_xH <- lapply(1:D, function(j) list(
    t = grow_init("integer", init_cap), #integer(0),
    ind  = grow_init("integer", init_cap), #integer(0),
    val  = grow_init("numeric", init_cap) #numeric(0)
  ))
  
  out_Pgen <- list(t = grow_init("integer", init_cap), child = grow_init("integer", init_cap), parent = grow_init("integer", init_cap))
  out_Hgen <- list(t = grow_init("integer", init_cap), child = grow_init("integer", init_cap), parent = grow_init("integer", init_cap))
  
  out_Pmut <- list(t = grow_init("integer", init_cap), ind = grow_init("integer", init_cap), count = grow_init("integer", init_cap))
  out_Hmut <- list(t = grow_init("integer", init_cap), ind = grow_init("integer", init_cap), count = grow_init("integer", init_cap))
  
  # creation of lists to help with the execution time
  # Phist=list(a=lapply(1:D,function(i){Matrix::Matrix(0,nrow=1,ncol=N,sparse=TRUE)}))
  # Hhist=list(a=lapply(1:D,function(i){Matrix::Matrix(0,nrow=1,ncol=N,sparse=TRUE)}))
  # for(j in 1:D){
  #   Phist[[1]][[j]][1,]=P[j,]
  #   Hhist[[1]][[j]][1,]=H[j,]
  # }
  # Pmut=list(a=Matrix::Matrix(0,nrow=1,ncol=N,sparse=TRUE))
  # Hmut=Pmut
  # Pgen=Pmut
  # Hgen=Pmut
  
  # maximum vector size needed to record a full inner loop. nP and nH are assumed to be equal. Need to be multiplied by 2 since the dead individuals are replaced by the same number of born individuals
  len_inner <- floor(sqrt(NG/(nH+nP))) * nP * 2
  # number of inner loops contained in a main loop
  len_outer <- max(1,ceiling(NG/(timeStep*thin)))
  
  Phist_start <- matrix(NA, nrow = D, ncol = len_outer)
  Phist_siz <- matrix(NA, nrow = D, ncol = len_outer)
  Phist_t <- matrix(NA, nrow = D, ncol = len_inner * len_outer)
  Phist_ind <- matrix(NA, nrow = D, ncol = len_inner * len_outer)
  Phist_val <- matrix(NA, nrow = D, ncol = len_inner * len_outer)
  
  Hhist_start <- matrix(NA, nrow = D, ncol = len_outer)
  Hhist_siz <- matrix(NA, nrow = D, ncol = len_outer)
  Hhist_t <- matrix(NA, nrow = D, ncol = len_inner * len_outer)
  Hhist_ind <- matrix(NA, nrow = D, ncol = len_inner * len_outer)
  Hhist_val <- matrix(NA, nrow = D, ncol = len_inner * len_outer)
  
  Pmut_start <- rep(NA, len_outer)
  Pmut_siz <- rep(NA, len_outer)
  Pmut_t <- rep(NA, len_inner * len_outer)
  Pmut_ind <- rep(NA, len_inner * len_outer)
  Pmut_count <- rep(NA, len_inner * len_outer)
  
  Hmut_start <- rep(NA, len_outer)
  Hmut_siz <- rep(NA, len_outer)
  Hmut_t <- rep(NA, len_inner * len_outer)
  Hmut_ind <- rep(NA, len_inner * len_outer)
  Hmut_count <- rep(NA, len_inner * len_outer)
  
  Pgen_start <- rep(NA, len_outer)
  Pgen_siz <- rep(NA, len_outer)
  Pgen_t <- rep(NA, len_inner * len_outer)
  Pgen_child <- rep(NA, len_inner * len_outer)
  Pgen_parent <- rep(NA, len_inner * len_outer)
  
  Hgen_start <- rep(NA, len_outer)
  Hgen_siz <- rep(NA, len_outer)
  Hgen_t <- rep(NA, len_inner * len_outer)
  Hgen_child <- rep(NA, len_inner * len_outer)
  Hgen_parent <- rep(NA, len_inner * len_outer)
  
  # preallocated lists to minimize memory use (list-based, avoid)
  # Phist=lapply(as.list(1:D),function(i) vector("list", max(1,ceiling(NG/(timeStep*thin))) ))
  # Hhist=lapply(as.list(1:D),function(i) vector("list", max(1,ceiling(NG/(timeStep*thin))) ))
  # Pmut=vector("list", max(1,ceiling(NG/(timeStep*thin))) )
  # Hmut=vector("list", max(1,ceiling(NG/(timeStep*thin))) )
  # Pgen=vector("list", max(1,ceiling(NG/(timeStep*thin))) )
  # Hgen=vector("list", max(1,ceiling(NG/(timeStep*thin))) )
  
  # fitness=sapply(1:N,function(i){p=P[,i];h=H[,i];c(f(p,h,alphaP,rP),f(p,h,alphaH,rH))}) # vectors of actual fitness
  
  listInd=1
  time=0
  t=0     # number of time steps executed in the loop
  
  start_Phist <- 1
  start_Hhist <- 1
  start_Pmut <- 1
  start_Hmut <- 1
  start_Pgen <- 1
  start_Hgen <- 1
  
  # beginning of the main loop
  for(u in 1:max(1,ceiling(NG/(timeStep*thin)))){
    
    NTime=ceiling((min(u*timeStep*thin,NG)-time)/thin) # number of time steps in one loop run
    
    t=0
    
    # creation of the matrices in which we record what happens 
    # phist=lapply(1:D,function(e){Matrix::Matrix(0,nrow=NTime,ncol=N,sparse=TRUE)})  # trait of P when there is a birth/death
    # hhist=phist         # trait of H when there is a birth/death
    # pmut=phist[[1]]     # number of mutations for P
    # hmut=phist[[1]]     # number of mutations for H
    # pgen=phist[[1]]     # parents for P
    # hgen=phist[[1]]     # parents for H
    
    # begining of the internal loop _ execution between two timeStep
    while(time<min(u*timeStep*thin,NG)){
      
      # initialisation
      if((t) %% thin ==0){
        parentsP=1:N
        parentsH=1:N
        pmut_act=rep(0,N)
        hmut_act=rep(0,N)
        changeP=c()
        changeH=c()
      }
      
      t=t+1
      time=time+1
      
      if(time %% verbose == 0 ) {cat("\r",time,"\r")}
      
      for(k in 1:nH){
        
        deadH=sample(1:N,1)    # which H individual are dying
        newH=0                 # which H individual giving birth
        
        
        # space=t(kernel[c(Y1,Y2),c(X1,X2)])
        
        p=P[,deadH]
        prob=f(H,p,rH,alphaH,deadH,dSpace_H)
        if(any(prob==Inf)|any(is.na(prob))){
          is.inf=(prob==Inf|is.na(prob))
          prob[!is.inf]=0.0000001
          prob[is.inf]=1}
        if(all(prob==0)) {
          y=floor((deadH-1)/nx)+1
          x=deadH-nx*(y-1)
          X1=x:1
          X2=2:(nx+1-x)
          if(x==nx) X2=c()
          Y1=y:1
          Y2=2:(ny+1-y)
          if(y==ny) Y2=c()
          prob=t(kernel[c(Y1,Y2),c(X1,X2)])}
        prob[prob==0]=min(prob[prob>0])/10
        newH=sample(1:N,1,prob = prob)
        
        mutH=(runif(1) < muH)     # mutation ?
        H[,deadH]=H[,newH]        # new H trait values
        changeH=unique(c(changeH,deadH,newH))    # which H changed since the last recording
        parentsH[deadH]=parentsH[newH]           # parents of new Hs
        hmut_act[deadH]=hmut_act[parentsH[newH]]+mutH  # number of mutations for H since last recording
        
        # modification of the trait by mutation
        # for H
        if(sum(mutH)>0){
          if(oneByOne){
            where=sample(D,sum(mutH),replace = TRUE)
            H[D*(deadH[mutH]-1)+where]=H[D*(deadH[mutH]-1)+where]+rnorm(sum(mutH),mean=0,sd=effect)
          }else{
            H[,deadH[mutH]]=H[,deadH[mutH]]+rnorm(D*sum(mutH),mean = 0,sd=effect)
          }
        }
      }
      
      for(k in 1:nP){
        
        deadP=sample(1:N,1)    # which H individual are dying
        newP=0                 # which H individual giving birth
        
        
        # space=t(kernel[c(Y1,Y2),c(X1,X2)])
        h=H[,deadP]
        prob=f(P,h,rP,alphaP,deadP,dSpace_P)
        prob[prob<0]
        if(any(prob==Inf)|any(is.na(prob))){
          is.inf=(prob==Inf|is.na(prob))
          prob[!is.inf]=0.0000001
          prob[is.inf]=1}
        if(all(prob==0)) {
          y=floor((deadP-1)/nx)+1
          x=deadP-nx*(y-1)
          X1=x:1
          X2=2:(nx+1-x)
          if(x==nx) X2=c()
          Y1=y:1
          Y2=2:(ny+1-y)
          if(y==ny) Y2=c()
          prob=t(kernel[c(Y1,Y2),c(X1,X2)])
        }
        prob[prob==0]=min(prob[prob>0])/10
        newP=sample(1:N,1,prob = prob)
        
        mutP=(runif(1) < muP)     # mutation ?
        P[,deadP]=P[,newP]        # new H trait values
        changeP=unique(c(changeP,deadP,newP))    # which H changed since the last recording
        parentsP[deadP]=parentsP[newP]           # parents of new Hs
        pmut_act[deadP]=pmut_act[parentsP[newP]]+mutP  # number of mutations for H since last recording
        
        # modification of the trait by mutation
        # for H
        if(sum(mutP)>0){
          if(oneByOne){
            where=sample(D,sum(mutP),replace = TRUE)
            P[D*(deadP[mutP]-1)+where]=P[D*(deadP[mutP]-1)+where]+rnorm(sum(mutP),mean=0,sd=effect)
          }else{
            P[,deadP[mutP]]=P[,deadP[mutP]]+rnorm(D*sum(mutP),mean = 0,sd=effect)
          }
        }
      }
      
      # recording
      if((t) %% thin ==0){
        # for(i in 1:D){
        #   phist[[i]][t/thin,changeP]=P[i,changeP]
        #   hhist[[i]][t/thin,changeH]=H[i,changeH]
        # }
        
        # compact version of the block above
        for(i in 1:D){
          if (length(changeP)) {
            out_xP[[i]]$t <- grow_append(out_xP[[i]]$t, rep(t, length(changeP)))
            out_xP[[i]]$ind  <- grow_append(out_xP[[i]]$ind,  changeP)
            out_xP[[i]]$val  <- grow_append(out_xP[[i]]$val,  P[i, changeP])
            
            # out_xP[[i]]$t <- out_xP[[i]]$t[out_xP[[i]]$val != 0]
            # out_xP[[i]]$ind <- out_xP[[i]]$ind[out_xP[[i]]$val != 0]
            # out_xP[[i]]$val <- out_xP[[i]]$val[out_xP[[i]]$val != 0]
          }
          out_xP[[i]]$dim <- c(NTime, N)
          
          if (length(changeH)) {
            out_xH[[i]]$t <- grow_append(out_xH[[i]]$t, rep(t, length(changeH)))
            out_xH[[i]]$ind  <- grow_append(out_xH[[i]]$ind,  changeH)
            out_xH[[i]]$val  <- grow_append(out_xH[[i]]$val,  H[i, changeH])
            
            # out_xH[[i]]$t <- out_xH[[i]]$t[out_xH[[i]]$val != 0]
            # out_xH[[i]]$ind <- out_xH[[i]]$ind[out_xH[[i]]$val != 0]
            # out_xH[[i]]$val <- out_xH[[i]]$val[out_xH[[i]]$val != 0]
          }
          out_xH[[i]]$dim <- c(NTime, N)
        }
        
        
        # pgen[t/thin,changeP]=parentsP[changeP]
        # hgen[t/thin,changeH]=parentsH[changeH]
        
        # compact version of the block above
        if (length(changeP)) {
          out_Pgen$t   <- grow_append(out_Pgen$t, rep(t, length(changeP)))
          out_Pgen$child  <- grow_append(out_Pgen$child, changeP)
          out_Pgen$parent <- grow_append(out_Pgen$parent, parentsP[changeP])
        }
        out_Pgen$dim <- c(NTime, N)
        
        if (length(changeH)) {
          out_Hgen$t   <- grow_append(out_Hgen$t, rep(t, length(changeH)))
          out_Hgen$child  <- grow_append(out_Hgen$child, changeH)
          out_Hgen$parent <- grow_append(out_Hgen$parent, parentsH[changeH])
        }
        out_Hgen$dim <- c(NTime, N)
        
        # pmut[t/thin,]=pmut_act
        # hmut[t/thin,]=hmut_act
        
        # compact version of the block above
        mut_indsP <- which(pmut_act > 0)
        if (length(mut_indsP)) {
          out_Pmut$t  <- grow_append(out_Pmut$t,  rep(t, length(mut_indsP)))
          out_Pmut$ind   <- grow_append(out_Pmut$ind,   mut_indsP)
          out_Pmut$count <- grow_append(out_Pmut$count, pmut_act[mut_indsP])
        }
        out_Pmut$dim <- c(NTime, N)
        
        mut_indsH <- which(hmut_act > 0)
        if (length(mut_indsH)) {
          out_Hmut$t  <- grow_append(out_Hmut$t,  rep(t, length(mut_indsH)))
          out_Hmut$ind   <- grow_append(out_Hmut$ind,   mut_indsH)
          out_Hmut$count <- grow_append(out_Hmut$count, hmut_act[mut_indsH])
        }
        out_Hmut$dim <- c(NTime, N)
        
      }
      
    }
    
    # the finalization step
    # Phist[[listInd]]=lapply(out_xP, function(x) list(t = grow_finalize(x$t), ind = grow_finalize(x$ind), val = grow_finalize(x$val), dim = x$dim))
    # Hhist[[listInd]]=lapply(out_xH, function(x) list(t = grow_finalize(x$t), ind = grow_finalize(x$ind), val = grow_finalize(x$val), dim = x$dim))
    # Pmut[[listInd]]=list(t = grow_finalize(out_Pmut$t), ind = grow_finalize(out_Pmut$ind), count = grow_finalize(out_Pmut$count), dim = out_Pmut$dim)
    # Hmut[[listInd]]=list(t = grow_finalize(out_Hmut$t), ind = grow_finalize(out_Hmut$ind), count = grow_finalize(out_Hmut$count), dim = out_Hmut$dim)
    # Pgen[[listInd]]=list(t = grow_finalize(out_Pgen$t), child = grow_finalize(out_Pgen$child), parent = grow_finalize(out_Pgen$parent), dim = out_Pgen$dim)
    # Hgen[[listInd]]=list(t = grow_finalize(out_Hgen$t), child = grow_finalize(out_Hgen$child), parent = grow_finalize(out_Hgen$parent), dim = out_Hgen$dim)
    
    start <- start_Phist
    siz <- length(grow_finalize(out_xP[[1]]$t))
    if(siz == 0){
      for (i in 1:D) {
        Phist_start[i, listInd] <- NA
        Phist_siz[i, listInd] <- NA
      }
    }else{
      for (i in 1:D) {
        Phist_start[i, listInd] <- start
        Phist_siz[i, listInd] <- siz
        Phist_t[i, start: (start+siz-1)] <- grow_finalize(out_xP[[i]]$t)
        Phist_ind[i, start: (start+siz-1)] <- grow_finalize(out_xP[[i]]$ind)
        Phist_val[i, start: (start+siz-1)] <- grow_finalize(out_xP[[i]]$val)
      }
      start_Phist <- start_Phist + siz
    }
    
    start <- start_Hhist
    siz <- length(grow_finalize(out_xP[[1]]$t))
    if(siz == 0){
      for (i in 1:D) {
        Hhist_start[i, listInd] <- NA
        Hhist_siz[i, listInd] <- NA
      }
    }else{
      for (i in 1:D) {
        Hhist_start[i, listInd] <- start
        Hhist_siz[i, listInd] <- siz
        Hhist_t[i, start: (start+siz-1)] <- grow_finalize(out_xP[[i]]$t)
        Hhist_ind[i, start: (start+siz-1)] <- grow_finalize(out_xP[[i]]$ind)
        Hhist_val[i, start: (start+siz-1)] <- grow_finalize(out_xP[[i]]$val)
      }
      start_Hhist <- start_Hhist + siz
    }
    
    start <- start_Pmut
    siz <- length(grow_finalize(out_Pmut$t))
    if(siz == 0){
      Pmut_start[listInd] <- NA
      Pmut_siz[listInd] <- NA
    }else{
      Pmut_start[listInd] <- start
      Pmut_siz[listInd] <- siz
      Pmut_t[start: (start+siz-1)] <- grow_finalize(out_Pmut$t)
      Pmut_ind[start: (start+siz-1)] <- grow_finalize(out_Pmut$ind)
      Pmut_count[start: (start+siz-1)] <- grow_finalize(out_Pmut$count)
      start_Pmut <- start_Pmut + siz
    }
    
    start <- start_Hmut
    siz <- length(grow_finalize(out_Hmut$t))
    if(siz == 0){
      Hmut_start[listInd] <- NA
      Hmut_siz[listInd] <- NA
    }else{
      Hmut_start[listInd] <- start
      Hmut_siz[listInd] <- siz
      Hmut_t[start: (start+siz-1)] <- grow_finalize(out_Hmut$t)
      Hmut_ind[start: (start+siz-1)] <- grow_finalize(out_Hmut$ind)
      Hmut_count[start: (start+siz-1)] <- grow_finalize(out_Hmut$count)
      start_Hmut <- start_Hmut + siz
    }

    
    start <- start_Pgen
    siz <- length(grow_finalize(out_Pgen$t))
    if(siz == 0){
      Pgen_start[listInd] <- NA
      Pgen_siz[listInd] <- NA
    }else{
      Pgen_start[listInd] <- start
      Pgen_siz[listInd] <- siz
      Pgen_t[start: (start+siz-1)] <- grow_finalize(out_Pgen$t)
      Pgen_child[start: (start+siz-1)] <- grow_finalize(out_Pgen$child)
      Pgen_parent[start: (start+siz-1)] <- grow_finalize(out_Pgen$parent)
      start_Pgen <- start_Pgen + siz
    }
    
    start <- start_Hgen
    siz <- length(grow_finalize(out_Hgen$t))
    if(siz == 0){
      Hgen_start[listInd] <- NA
      Hgen_siz[listInd] <- NA
    }else{
      Hgen_start[listInd] <- start
      Hgen_siz[listInd] <- siz
      Hgen_t[start: (start+siz-1)] <- grow_finalize(out_Hgen$t)
      Hgen_child[start: (start+siz-1)] <- grow_finalize(out_Hgen$child)
      Hgen_parent[start: (start+siz-1)] <- grow_finalize(out_Hgen$parent)
      start_Hgen <- start_Hgen + siz
    }
    
    # re-initialize the intermediate variables
    out_Pgen$t <- reset_grow(out_Pgen$t)
    out_Pgen$child <- reset_grow(out_Pgen$child)
    out_Pgen$parent <- reset_grow(out_Pgen$parent)
    out_Pgen$dim <- c(NTime, N)
    out_Hgen$t <- reset_grow(out_Hgen$t)
    out_Hgen$child <- reset_grow(out_Hgen$child)
    out_Hgen$parent <- reset_grow(out_Hgen$parent)
    out_Hgen$dim <- c(NTime, N)
    
    for(i in 1:D){
      out_xP[[i]]$t <- reset_grow(out_xP[[i]]$t)
      out_xP[[i]]$ind <- reset_grow(out_xP[[i]]$ind)
      out_xP[[i]]$val <- reset_grow(out_xP[[i]]$val)
      out_xP[[i]]$dim <- c(NTime, N)
      out_xH[[i]]$t <- reset_grow(out_xH[[i]]$t)
      out_xH[[i]]$ind <- reset_grow(out_xH[[i]]$ind)
      out_xH[[i]]$val <- reset_grow(out_xH[[i]]$val)
      out_xH[[i]]$dim <- c(NTime, N)
    }
    
    out_Pmut$t <- reset_grow(out_Pmut$t)
    out_Pmut$ind <- reset_grow(out_Pmut$ind)
    out_Pmut$count <- reset_grow(out_Pmut$count)
    out_Pmut$dim <- c(NTime, N)
    out_Hmut$t <- reset_grow(out_Hmut$t)
    out_Hmut$ind <- reset_grow(out_Hmut$ind)
    out_Hmut$count <- reset_grow(out_Hmut$count)
    out_Hmut$dim <- c(NTime, N)
    
    listInd=listInd+1
    
  }
  
  Phist=list(
    start = Phist_start,
    siz = Phist_siz,
    t = Phist_t,
    ind = Phist_ind,
    val = Phist_val
  )
  Hhist=list(
    start = Hhist_start,
    siz = Hhist_siz,
    t = Hhist_t,
    ind = Hhist_ind,
    val = Hhist_val
  )
  Pmut=list(
    start = Pmut_start,
    siz = Pmut_siz,
    t = Pmut_t,
    ind = Pmut_ind,
    count = Pmut_count
  )
  Hmut=list(
    start = Hmut_start,
    siz = Hmut_siz,
    t = Hmut_t,
    ind = Hmut_ind,
    count = Hmut_count
  )
  Pgen=list(
    start = Pgen_start,
    siz = Pgen_siz,
    t = Pgen_t,
    child = Pgen_child,
    parent = Pgen_parent
  )
  Hgen=list(
    start = Hgen_start,
    siz = Hgen_siz,
    t = Hgen_t,
    child = Hgen_child,
    parent = Hgen_parent
  )

  return(list(
    Pgenealogy=Pgen,
    Hgenealogy=Hgen,
    xP=Phist,
    xH=Hhist,
    P=P,
    H=H,
    Pmut=Pmut,
    Hmut=Hmut,
    iniP=iniP,
    iniH=iniH,
    thin.factor=thin,
    N = N,
    NG = NG,
    timeStep = timeStep,
    D = D
    ))
  
}


######## Construct the genealogy ###################################

# out an object created by the simulation function
# treeP, treeH former genealogy the tree have to be graphted on

make_gen.BipartiteEvol=function(out, treeP=NULL, treeH=NULL, verbose=TRUE){
  
  N=out$N     # number of individuals
  NG=out$NG    # total time steps number
  D=out$D     # dimension of trait space
  
  # auxiliary function creating a genealogy for one type
  aux=function(P.trait,Gen,X,Mut,ini,thin,tree, verbose){
    
    edgeP=matrix(ncol=2)                       # the edges of the genealogy tree
    xP=lapply(1:D,function(i){c()})            # trait values at the beginning of the edges
    currentP=1:N                               # individuals at current time that have descendants in the present
    node=N+1                                   # current node name (1:N are the tips)
    nodes=1:N                                  # current nodes
    nodes_age=rep(NG+1,N)                      # current nodes age
    edge.length=c()                            # length of the edges of the genealogy
    nMut=rep(0,N)                              # number of mutation between a node and its parent node
    t=NG                                       # actual time
    
    # listInd=length(Mut)+1                    # indice of the current matrix in the list
    listInd=length(Mut$start)+1            # indice of the current matrix in the list
    time=NG                                    # time at the first row of the current matrix
    newInd=NULL                                # who was born at time t
    
    # begining of main loop
    while (t>0 & length(currentP)>1){
      # ie we are not yet at the root and there is more than one individual with descendant in the present
      
      listInd=listInd-1
      # time.step=nrow(Mut[[listInd]])
      time.step=out$timeStep
      time=time-time.step
      if(time<0){
        time.step=time.step+time
        time=0
      }
      
      # x=lapply(1:D,function(i){X[[listInd]][[i]]})  # trait values
      # gen=Gen[[listInd]]                            # parents
      # mut=Mut[[listInd]]                            # mutation number
      
      # matric version
      # trait values (here each row correspond to a dimension in the trait space)
      x_start=X$start[,listInd]
      x_siz=X$siz[,listInd]
      x_t=X$t[,x_start[1]:(x_start[1]+x_siz[1]-1)]
      x_ind=X$ind[,x_start[1]:(x_start[1]+x_siz[1]-1)]
      x_val=X$val[,x_start[1]:(x_start[1]+x_siz[1]-1)]
      # parents
      gen_start=Gen$start[listInd]
      gen_siz=Gen$siz[listInd]
      gen_t=Gen$t[gen_start:(gen_start+gen_siz-1)]
      gen_child=Gen$child[gen_start:(gen_start+gen_siz-1)]
      gen_parent=Gen$parent[gen_start:(gen_start+gen_siz-1)]
      # mutation number
      mut_start=Mut$start[listInd]
      mut_siz=Mut$siz[listInd]
      mut_t <- integer(0)
      mut_ind <- integer(0)
      mut_count <- integer(0)
      if(!is.na(mut_start) && !is.na(mut_siz)){
        mut_t=Mut$t[mut_start:(mut_start+mut_siz-1)]
        mut_ind=Mut$ind[mut_start:(mut_start+mut_siz-1)]
        mut_count=Mut$count[mut_start:(mut_start+mut_siz-1)]
      }
      
      # internal loop
      while(t>time){
        # ie we are not at the first row of the matrix yet
        
        prev.t=t # ???
        
        if(verbose) cat(paste("\r",length(currentP),t,"\r"))
        
        if(is.null(newInd)){
          # newInd=which(gen[t-time,currentP]>0)
          newInd=which(currentP %in% gen_child[which(gen_t == t-time)])
          
        }
        
        # if there was a birth at time t
        if (sum(newInd)>0){
          newCurrentP=currentP
          newNodes=nodes
          # fathers=gen[t-time,currentP[newInd]]    # parents of the new individuals
          if(length(which(gen_t == t-time)) > 0){
            if(length(which(gen_child %in% currentP[newInd])) > 0){
              temp_ind <- intersect(which(gen_t == t-time), which(gen_child %in% currentP[newInd]))
            }
          }
          fathers=gen_parent[temp_ind]    # parents of the new individuals
          
          for(i in 1:length(newInd)){
            father=fathers[i]
            newCurrentP[newInd[i]]=father
            # nMut[nodes[newInd[i]]]=nMut[nodes[newInd[i]]]+(mut[t-time,currentP])[newInd[i]]
            temp <- 0
            if(length(which(mut_t == t-time)) > 0){
              if(length(which(mut_ind == currentP[newInd[i]])) > 0){
                temp_ind <- intersect(which(mut_t == t-time), which(mut_ind == currentP[newInd[i]]))
                if(length(temp_ind) == 1){
                  temp <- mut_count[temp_ind]
                }else{
                  print("wrong")
                }
              }
            }
            nMut[nodes[newInd[i]]]=nMut[nodes[newInd[i]]]+temp
            
            if(sum(fathers[1:i]==father)==1){
              # ie the individual do not have the same parent than another born at this time already treated
              
              if(!(father %in% currentP[newInd]) & father %in% currentP){
                # ie the parent is in currentP (else there is a pb !) and did not die at time t 
                
                # we add a new node 
                newNodes[newInd[i]]=node
                newNodes[currentP==father]=node
                nodes_age=c(nodes_age,t)
                nMut=c(nMut,0)
                
                # and the corresponding offspring edges
                # edge 1
                edgeP=rbind(edgeP,c(node,nodes[newInd[i]]))
                for(j in 1:D){
                  
                  temp <- 0
                  if(length(which(x_t[1,] == t-time)) > 0){
                    if(length(which(x_ind[1,] == currentP[newInd[i]])) > 0){
                      temp_ind <- intersect(which(x_t[1,] == t-time), which(x_ind[1,] == currentP[newInd[i]]))
                      if(length(temp_ind) > 0){
                        temp <- x_val[j,temp_ind]
                      }
                    }
                  }
                  
                  xP[[j]]=c(
                    xP[[j]],
                    # (x[[j]][t-time,currentP])[newInd[i]]
                    temp
                    )
                }
                edge.length=c(edge.length,nodes_age[nodes[newInd[i]]]-t)
                
                # edge 2
                edgeP=rbind(edgeP,c(node,nodes[currentP==father]))
                for(j in 1:D){
                  
                  temp <- 0
                  if(length(which(x_t[1,] == t-time)) > 0){
                    if(length(which(x_ind[1,] == father)) > 0){
                      temp_ind <- intersect(which(x_t[1,] == t-time), which(x_ind[1,] == father))
                      if(length(temp_ind) > 0){
                        temp <- x_val[j,temp_ind]
                      }
                    }
                  }
                  
                  xP[[j]]=c(
                    xP[[j]],
                    # x[[j]][t-time,father]
                    temp
                    )
                }
                edge.length=c(edge.length,nodes_age[nodes[currentP==father]]-t)
                
                node=node+1
                
              }else{
                # ie the parent died at time t
                
                if(sum(fathers==father)>1){
                  # ie there is another current individual born from the same partent at time t
                  
                  # we add a new node
                  newNodes[newInd[i]]=node
                  nodes_age=c(nodes_age,t)
                  nMut=c(nMut,0)
                  
                  # and the first of its offspring edges
                  edgeP=rbind(edgeP,c(node,nodes[newInd[i]]))
                  for(j in 1:D){
                    
                    temp <- 0
                    if(length(which(x_t[1,] == t-time)) > 0){
                      if(length(which(x_ind[1,] == currentP[newInd[i]])) > 0){
                        temp_ind <- intersect(which(x_t[1,] == t-time), which(x_ind[1,] == currentP[newInd[i]]))
                        if(length(temp_ind) > 0){
                          temp <- x_val[j,temp_ind]
                        }
                      }
                    }
                    
                    xP[[j]]=c(
                      xP[[j]],
                      # (x[[j]][t-time,currentP])[newInd[i]])
                      temp
                    )
                      
                  }
                  edge.length=c(edge.length,nodes_age[nodes[newInd[i]]]-t)
                  
                  node=node+1
                }
                # else we don't do anything, the current individual identity is only changed in current P (already done)
              }
            }else{
              # ie there is an individual born at time t (and already seen) with the same parent
              
              if(!(father %in% currentP[newInd]) & father %in% currentP){
                # ie the parent did not die at time t and is in currentP (always, else there is a pb)
                # we select the parent node
                Node=newNodes[currentP==father]
              }else{
                # ie the parent died at time t
                # we select the parent node
                Node=(newNodes[newCurrentP==father])[1]
              }
              
              #and add the corresponding edge
              edgeP=rbind(edgeP,c(Node,nodes[newInd[i]]))
              for(j in 1:D){
                
                temp <- 0
                if(length(which(x_t[1,] == t-time)) > 0){
                  if(length(which(x_ind[1,] == currentP[newInd[i]])) > 0){
                    temp_ind <- intersect(which(x_t[1,] == t-time), which(x_ind[1,] == currentP[newInd[i]]))
                    if(length(temp_ind) > 0){
                      temp <- x_val[j,temp_ind]
                    }
                  }
                }
                
                xP[[j]]=c(
                  xP[[j]],
                  # (x[[j]][t-time,currentP])[newInd[i]])
                  temp
                )
                
              }
              newNodes[newInd[i]]=Node
              edge.length=c(edge.length,nodes_age[nodes[newInd[i]]]-t)
            }
          }
          
          currentP=unique(newCurrentP)
          nodes=unique(newNodes)
        }
        
        # select next time
        if(t==(time+1)){
          # ie the matrix has no line left
          t=t-1
          newInd=NULL
        }else{
          # newInd=which(gen[t-time-1,currentP]>0)
          newInd=which(currentP %in% gen_child[which(gen_t == t-time-1)])
          if(sum(newInd)>0 ){
            # ie there is a change in the next time step
            t=t-1
          }else{
            # we look for the next change
            newInd=NULL
            if(t>1 & length(currentP)>1){
              # mat=Matrix::which(comp_gen[1:(t-time-1),currentP]>0,arr.ind = TRUE)
              ordr <- order(gen_child[gen_t >= 1 & gen_t <= t-time-1])
              row <- gen_t[gen_t >= 1 & gen_t <= t-time-1][ordr]
              col <- gen_child[gen_t >= 1 & gen_t <= t-time-1][ordr]
              mat <- cbind(row, col)
              
              if(length(mat)==0) {
                t=time
              }else{
                t=max(mat[,1])+time
                ind=mat[mat[,1]==(t-time),2]
              }
            }else{
              t=time
            }}}
      } # end internal loop
    } # end main loop
    
    edgeP=edgeP[-1,]   # edgeP was initialized with a row of NA, we suppress it
    
    if(length(currentP)>1){
      # ie all the individual did not coalesce
      
      if(is.null(tree)){
        for(i in 1:length(currentP)){
          edgeP=rbind(edgeP,c(node,nodes[i]))
          for(j in 1:D){xP[[j]]=c(xP[[j]],ini)}
          edge.length=c(edge.length,nodes_age[nodes[i]])
        }
        
      }else{
        # we grapht the new tree on the former tree
        M=max(max(edgeP),N)
        newmut=rep(0,max(tree$edge))
        for(j in 1:nrow(tree$edge)){ newmut[tree$edge[j,2]]=tree$nMut[j]}
        nMut=c(nMut,newmut)
        tree$edge=tree$edge+M
        edgeP=rbind(edgeP,tree$edge)
        edge.length=c(edge.length,tree$edge.length/thin)
        for(j in 1:D){
          xP[[j]]=c(
            xP[[j]],
            tree$x[[j]]
            )
        }
        for (i in 1:length(currentP)){
          ind=which(tree$tip.label==currentP[i])+M
          ind.edge=which(edgeP[,2]==ind)
          nMut[nodes[i]]=nMut[nodes[i]]+nMut[edgeP[ind.edge,2]]
          edgeP[ind.edge,2]=nodes[i]
          edge.length[ind.edge]=edge.length[ind.edge]+nodes_age[nodes[i]]-1
        }
      }
    }
    
    # create a phylo object
    if(is.null(tree) | length(currentP)==1){
      treeP=list(Nnode=max(edgeP)-N,edge=edgeP,tip.label=1:N,edge.length=(edge.length)*thin)
    }else{
      treeP=list(Nnode=length(unique(as.vector(edgeP)))-(2*N-length(currentP)),edge=edgeP,tip.label=c(1:(max(edgeP))),edge.length=(edge.length)*thin)
    }
    class(treeP)="phylo"
    X=rep(0,nrow(edgeP))
    nMut=nMut[treeP$edge[,2]]
    traits=xP
    traits[[D+1]]=nMut
    if(verbose) cat(".")
    
    # make the tree well conformed and suppress extinct lineages
    treeP=rigth.order.BE(treeP,traits)
    if(!is.null(tree) & (2*N-length(currentP))>N){
      extinct=treeP$tree$tip.label[treeP$tree$tip.label>N]
      treeP=prune.with.traits(treeP$tree,extinct,treeP$trait,extensiveTrait = D+1)}
    if(verbose) cat(".")
    
    # add the traits in the phylo object
    xP=treeP$trait[1:D]
    nMut=treeP$trait[[D+1]]
    x.tip=matrix(0,ncol=N,nrow=D)
    for(j in 1:D){
      x.tip[j,]=P.trait[j,treeP$tree$tip.label]
    }
    if(verbose) cat(".")
    treeP=treeP$tree
    treeP$x=xP
    treeP$ini=max(node.depth.edgelength(treeP))
    treeP$nMut=nMut
    treeP$x.tip=x.tip
    
    return(treeP)
  }
  
  if(verbose) print("P")
  P=aux(out$P, out$Pgenealogy, out$xP, out$Pmut, out$iniP, out$thin.factor, treeP, verbose=verbose)
  if(verbose) print("H")
  H=aux(out$H, out$Hgenealogy, out$xH, out$Hmut, out$iniH, out$thin.factor, treeH, verbose=verbose)
  return(list(P=P,H=H))
}


######## Create a species tree from the genealogy #################

# genealogy an object created with make.gen
# threshold number of mutation to be in different species
# distanceP, distanceH distance (ie nb of mutations) matrix between the individual
# if NULL it is computed within the function


define_species.BipartiteEvol=function(genealogy,threshold=1,distanceH=NULL,distanceP=NULL, verbose=TRUE, monophyly=TRUE, seed=NULL){
  
  if (monophyly==TRUE){
    
    D=length(genealogy$P$x)       # dimension of trait space
    
    # auxiliary function to define the species for one type
    aux = function(gen, distance, verbose=TRUE){
      N=length(gen$tip.label)    # number of individuals
      
      # first we define the genetic types of the individuals
      
      if(threshold==1){
        type=rep(1,N)             
        new.type=2
        for(i in 1:nrow(gen$edge)){
          if(gen$nMut[i]>0){
            if(gen$edge[i,2]<=N){
              type[gen$tip.label[gen$edge[i,2]]]=new.type
              new.type=new.type+1
            }else{
              type[as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)]=new.type}
            new.type=new.type+1
          }
        }
      }else{ 
        
        if(is.null(distance)){
          
          # we compute the distance matrix if it is not given as an argument
          distance=matrix(0,nrow=N,ncol=N)
          for(i in 1:nrow(gen$edge)){
            if(gen$nMut[i]>0){
              #add the number of mutation between the individuals separated by this edge
              if (verbose) cat("\r",i,"\r")
              if(gen$edge[i,2]<=N){
                tip=as.integer(gen$tip.label[gen$edge[i,2]])
                distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
                distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
              }else{
                tip=as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)
                distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
                distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
              }
            }
          }
          distance=distance[as.integer(gen$tip.label),as.integer(gen$tip.label)]
        }
        
        # we use the distance matrix to define the genetic type of individuals
        phylo=gen
        phylo$edge.length=phylo$nMut
        type=1:N
        for(i in 1:N){
          ind=((1:i)[distance[i,1:i]<threshold])[1]
          type[phylo$tip.label[i]]=type[phylo$tip.label[ind]]
          if (verbose) cat("\r",i,"\r")
        }
      }
      
      # species are then defined as monophyletic types
      species=rep(1,N)        # the species of each tip
      newSpecies=2
      nextNode=N+1
      
      # test of the tips descendig from a particular node _ main loop
      while(length(nextNode)>0){
        i=nextNode[1]                                       # candidate node
        nextNode=nextNode[-1]                               # update candidate nodes vector
        
        a=extract.clade(gen,i)                              # subtree with root i
        ty=as.integer(gen$tip.label)
        ty=ty[which(!(ty %in% as.integer(a$tip.label)))]
        ty=type[ty]                                         # type of tips not in a
        
        if(length(intersect(type[as.integer(a$tip.label)],ty))==0){
          # seeing how nextNode is built, this should always be the case... Test needed
          n=length(a$tip.label)                             # number of tips in the subtree
          offspring=gen$edge[gen$edge[,1]==i,2]             # offspring nodes of i in the main tree
          offspringSubTree=a$edge[a$edge[,1]==(n+1),2]      # offspring node of i in the subtree
          
          if(length(offspringSubTree)>1){
            
            for(j in 1:(length(offspringSubTree))){
              
              if(offspringSubTree[j]<=n){
                # ie offspringSubTree[j] is a tip
                
                ty=as.integer(a$tip.label)
                ty=ty[which(!(ty %in% as.integer(a$tip.label[offspringSubTree[j]])))]
                ty=type[ty]                                  # type of tips in a that are not offspring[j]
                if(!(type[as.integer(a$tip.label[offspringSubTree[j]])] %in% ty)){
                  # offspringSubTree[j] is a species
                  species[as.integer(a$tip.label[offspringSubTree[j]])]=newSpecies
                  newSpecies=newSpecies+1
                }
                
              }else{
                # ie offspringSubTree[j] is an internal node
                a1=extract.clade(a,offspringSubTree[j])
                ty=as.integer(a$tip.label)                  
                ty=ty[which(!(ty %in% as.integer(a1$tip.label)))]
                ty=type[ty]                                # type of tips in a but not in a1
                if(length(intersect(type[as.integer(a1$tip.label)],ty))==0){
                  # the tips in a1 are in a different species than the other tips in a
                  nextNode=c(nextNode,offspring[j])        # new candidate node
                  species[as.integer(a1$tip.label)]=newSpecies
                  newSpecies=newSpecies+1
                }
              }
            }
          }else{if (verbose) print("error : one node has only one offspring... ???")}
        }else{if (verbose) print(paste("error : this node should not be in next node :",i))}
      }   # end of main loop
      
      # rename the species so that there is no gap in the names
      M=rep(0,max(species))
      name=1
      for(i in 1:length(species)){
        if(M[species[i]]==0){
          M[species[i]]=name
          name=name+1
        }
        species[i]=M[species[i]]
      }
      
      return(species)
    }
    
    # auxiliary function to build the species tree from the genealogy and the species of each tip
    make.phylo=function(gen,spec){
      
      N=length(gen$tip.label)                                                                             # number of individuals
      abundance=sapply(unique(spec),function(i){sum(spec==i)})                                            # abundance of each species
      mean.trait=sapply(unique(spec),function(i){sapply(1:D,function(e){mean(gen$x.tip[e,spec[gen$tip.label]==i])})})    # mean trait of each species
      if(D == 1){mean.trait = matrix(mean.trait,nrow = D)}
      
      if(max(spec)==1){return(list(abundance=abundance,mean.trait=mean.trait))}
      
      keep.tip=c()                                                                                        # which tips will be in the species tree
      keep.tip=sapply(unique(spec),function(i){vect=(1:N)[spec[gen$tip.label]==i];
      if(length(vect)==1){return(vect)};
      sample(vect,1)})
      
      tree=prune.with.traits(gen,gen$tip.label[-keep.tip],gen$x)               # suppress the other tips
      tree$abundance=abundance[spec[as.integer(tree$tree$tip.label)]]
      tree$mean.trait=matrix(mean.trait[,spec[as.integer(tree$tree$tip.label)]], nrow = D, byrow = FALSE)
      tree$tree$tip.label=spec[as.integer(tree$tree$tip.label)]
      return(tree)
    }
    
    P=aux(genealogy$P, distanceP, verbose=verbose)
    H=aux(genealogy$H, distanceH, verbose=verbose)
    Pphylo=make.phylo(genealogy$P, P)
    Hphylo=make.phylo(genealogy$H, H)
    
  } else {   # Monophyly == False
    
    D=length(genealogy$P$x)       # dimention of trait space
    
    # auxiliary function to define the species for one type
    aux = function(gen,distance, verbose=TRUE){
      N=length(gen$tip.label)    # number of individuals
      
      # first we define the genetic types of the individuals
      
      if(threshold!=1){
        if (verbose) print("Threshold set to 1 when monophyly==FALSE")
        threshold <- 1}
      
      
      if(threshold==1){
        type=rep(1,N)             
        new.type=2
        for(i in 1:nrow(gen$edge)){
          if(gen$nMut[i]>0){
            if(gen$edge[i,2]<=N){
              type[gen$tip.label[gen$edge[i,2]]]=new.type
              new.type=new.type+1
            }else{
              type[as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)]=new.type}
            new.type=new.type+1
          }
        }
      }
      
      # rename the species (=type) so that there is no gap in the names
      M=rep(0,max(type))
      name=1
      for(i in 1:length(type)){
        if(M[type[i]]==0){
          M[type[i]]=name
          name=name+1
        }
        type[i]=M[type[i]]
      }
      
      return(type)
    }
    
    # auxiliary function to build the species tree from the genealogy and the species of each tip
    make.phylo=function(gen,spec){
      
      N=length(gen$tip.label)  
      list_species <- unique(spec)
      unique_type <- c()
      pruned_gen <- gen
      pruned_gen$tip.label <- as.character(pruned_gen$tip.label)
      pruned_gen <- multi2di(pruned_gen)
      
      if (!is.null(seed)) {
        set.seed(seed)
        sampled_tips <- sample(1:N)} else {sampled_tips <- 1:N}
      
      for (i in sampled_tips){
        if (spec[i] %in% unique_type){
          pruned_gen <- drop.tip(pruned_gen,tip=as.character(i), trim.internal = TRUE)
        }else{unique_type <- c(unique_type, spec[i])}
      }
      
      pruned_gen$tip.label <- spec[as.numeric(pruned_gen$tip.label)]
      
      abundance=sapply(unique(spec),function(i){sum(spec==i)})   # abundance of each species
      mean.trait=sapply(unique(spec),function(i){sapply(1:D,function(e){mean(gen$x.tip[e,spec[gen$tip.label]==i])})})    # mean trait of each species
      if(D == 1){mean.trait = matrix(mean.trait,nrow = D)}
      if(max(spec)==1){return(list(abundance=abundance,mean.trait=mean.trait))}
      
      
      tree=list(tree=pruned_gen) 
      tree$abundance=abundance[tree$tree$tip.label]
      tree$mean.trait=matrix(mean.trait[,tree$tree$tip.label], nrow = D, byrow = FALSE)
      # tree$mean.trait and tree$abundance order like tree$tree$tip.label
      return(tree)
    }
    
    P=aux(genealogy$P,distanceP, verbose=verbose)
    H=aux(genealogy$H,distanceH, verbose=verbose)
    Pphylo=make.phylo(genealogy$P, P)
    Hphylo=make.phylo(genealogy$H, H)
    
  }
  
  return(list(P=P,H=H,Pphylo=Pphylo,Hphylo=Hphylo))
}


######## Build the network ##################

build_network.BipartiteEvol=function(gen, spec){
  
  
  P=spec$P                                         # what Pspecies is each P individual in ?
  H=spec$H                                         # what Hspecies is each H individual in ?
  N=length(P)                                      # number of individual in each clade
  nSP=max(P)                                       # number of Pspecies
  nSH=max(H) 
  
  link=Matrix::Matrix(0,nrow = nSP,ncol = nSH,sparse = TRUE)  # initialize the link matrix
  
  for(i in 1:N){
      link[P[i],H[i]]=link[P[i],H[i]]+1
  }
  
  return(link)
}




