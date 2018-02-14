################################################################################
##                                                                            ##
##                                RPANDA : Utils                              ##
##                                                                            ##
##   Julien Clavel - 01-02-2018                                               ##
##   S3 methods, simulations, miscellaneous                                   ##
##                                                                            ##
################################################################################

# make a specific S3 method "ancestral" for reconstructing or retrieving ancestral states (see phyl.pca_pl.R)
ancestral <- function(object) UseMethod("ancestral")

# S3 for the fit_t_env class
ancestral.fit_t.env <- function(object){
    
    # extract objects
    if(!inherits(object,"fit_t.env")) stop("only works with \"fit_t.env\" class objects. See ?fit_t_env")
    #tree <- object$scaled_tree
    #n <- object$n
    #Y <- object$Y
    
    # covariance for the nodes
    #V <- .vcvPhyloInternal(tree)
    #indice <- (1:n)
    #AY <- V[-indice,indice]
    #vY <- V[indice,indice]
    
    # Ancestral state at the root
    a <- object$root
    names(a) = "root"
    #one <- matrix(1,ncol=1,nrow=n)
    
    # states at the nodes
    #recons_t <- (AY%*%pseudoinverse(vY)%*%(Y-one%*%a))+(one[-1]%*%a)
    #colnames(recons_t) = colnames(Y)
    #rownames(recons_t) = paste("node_",n+1:Nnode(tree), sep="")
    
    # return the results
    res <- list(root=a, nodes=recons_t)
 return(res)
 warning("only the root state is currently estimated for models of the class \"fit_t.env\"") # To remove later

}

# S3 for the fit_t_comp class?
# TODO
ancestral.fit_t.comp <- function(object){
    
    # extract objects
    if(!inherits(object,"fit_t.comp")) stop("only works with \"fit_t.comp\" class objects. See ?fit_t_comp")
    
    anc <- object$z0
    names(anc) ="root"
    return(root)
    warning("only the root state is currently estimated for models of the class \"fit_t.comp\"")
}

# Build a matrix with tip and internal covariances
.vcvPhyloInternal <- function(tree){
    nbtip <- Ntip(tree)
    dis <- dist.nodes(tree)
    MRCA <- mrca(tree, full = TRUE)
    M <- dis[as.character(nbtip + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    return(M)
}

# ---- Function to simulate random covariance matrices with a specified eigenstructure
# From Uyeda et al. 2015 - Systematic Biology 64(4):677-689.
Posdef <- function (p, ev = rexp(p, 1/100)) {
  Z <- matrix(ncol=p, rnorm(p^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}
