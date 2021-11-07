################################################################################
##                                                                            ##
##                                RPANDA : Utils                              ##
##                                                                            ##
##   Julien Clavel - 01-02-2018                                               ##
##   S3 methods, simulations, miscellaneous                                   ##
##                                                                            ##
################################################################################

# S3 generic method "ancestral" for reconstructing or retrieving ancestral states (see phyl.pca_pl.R) # should I use "predict" instead?
ancestral <- function(object) UseMethod("ancestral")

# S3 for the fit_t_env class
ancestral.fit_t.env <- function(object){
    
    # extract objects
    if(!inherits(object,"fit_t.env")) stop("only works with \"fit_t.env\" class objects. See ?fit_t_env")
    
    # Ancestral state at the root
    a <- object$root
    names(a) = "root"

    res <- a
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
    return(anc)
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

# Build the matrix square root and inverse matrix square root of the phylogenetic covariance matrix
.transformsqrt <- function(tree){
    vcv_tr <- vcv.phylo(tree)
    eig <- eigen(vcv_tr)
    sqrtmValues <- sqrt(eig$values)
    sqrtM1 <- tcrossprod(eig$vectors%*%diag(1/sqrtmValues), eig$vectors)
    sqrtM <- tcrossprod(eig$vectors%*%diag(sqrtmValues), eig$vectors)
    matrices <- list(sqrtM1=sqrtM1, sqrtM=sqrtM)
    return(matrices)
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

# --- Function to simulate multivariate normal distribution
# From the mvtnorm package H.Morlon
rmvnorm_util<-function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, 
    checkSymmetry = TRUE) 
{
    if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) 
        stop("mean and sigma have non-conforming size")
    method <- match.arg(method)
    R <- if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive semidefinite")
        }
        t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 
            0))))
    }
    else if (method == "svd") {
        s. <- svd(sigma)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
            warning("sigma is numerically not positive semidefinite")
        }
        t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
    }
    else if (method == "chol") {
        R <- chol(sigma, pivot = TRUE)
        R[, order(attr(R, "pivot"))]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% 
        R
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}

# --- Function to compute node ages
# From the picante package H.Morlon
node.age_util<-function (phy) 
{
    if (!is.rooted(phy)) {
        stop("function node.age only works with a rooted phylo object")
    }
    rootN = length(phy$tip.label) + 1
    nEdges = nrow(phy$edge)
    ages = rep(NA, nEdges)
    for (n in 1:nEdges) {
        if (phy$edge[n, 1] == rootN) 
            anc.age = 0
        else {
            anc.age = ages[which(phy$edge[, 2] == phy$edge[n, 
                1])]
        }
        ages[n] = anc.age + phy$edge.length[n]
    }
    phy$ages = ages
    return(phy)
}


## Wrapper to compute the inverse fourier transform as in pracma using fft from "stats"
## J. Clavel
ifft_panda <- function(x) return(fft(x, inverse=TRUE)*(1/length(x)))

## Compute the pseudoinverse from svd
## J. Clavel modified from corpcor
## most of the pseudoinverse are applied to square matrices in RPANDA. No need for fat/thin matrices implementations
pseudoinverse = function (m, tol)
{

    # compute the svd
    s = svd(m)
    
    if( missing(tol) )
    tol = max(dim(m))*max(s$d)*.Machine$double.eps
    Positive = s$d > tol
    
    if (length(s$d) == 0)
    {
        return(
        array(0, dim(m)[2:1])
        )
    }
    else
    {
        return(
        s$v[, Positive, drop=FALSE] %*% (1/s$d[Positive] * t(s$u[, Positive, drop=FALSE]))
        )
    }
}

## detect extinct and remove them - from geiger 2.0.7
## J. Clavel
drop.extinct <- function (phy, tol=NULL) {
    if (!"phylo" %in% class(phy)) {
        stop("\"phy\" is not of class \"phylo\".");
    }
    if (is.null(phy$edge.length)) {
        stop("\"phy\" does not have branch lengths.");
    }
    if (is.null(tol)) {
        tol <- min(phy$edge.length)/100;
    }
    aa <- is.extinct(phy=phy, tol=tol);
    if (length(aa) > 0) {
        phy <- drop.tip(phy, aa); # use drop.tip from "ape" => Imports
    }
    return(phy);
}

# return tip.labels, so that tree ordering is not an issue
is.extinct <- function (phy, tol=NULL) {
    if (!"phylo" %in% class(phy)) {
        stop("\"phy\" is not of class \"phylo\".");
    }
    if (is.null(phy$edge.length)) {
        stop("\"phy\" does not have branch lengths.");
    }
    if (is.null(tol)) {
        tol <- min(phy$edge.length)/100;
    }
    phy <- reorder(phy);
    xx <- numeric(Ntip(phy) + phy$Nnode);
    for (i in 1:length(phy$edge[,1])) {
        xx[phy$edge[i,2]] <- xx[phy$edge[i,1]] + phy$edge.length[i];
    }
    aa <- max(xx[1:Ntip(phy)]) - xx[1:Ntip(phy)] > tol;
    if (any(aa)) {
        return(phy$tip.label[which(aa)]);
    } else {
        return(NULL);
    }
}

## EB transform
## part of the code used in fit_t_pl
## J. Clavel

transform_EB <- function(phy, beta, sigmasq){
    
    # parents & descent
    parent <- phy$edge[,1]
    descendent <- phy$edge[,2]

    if (beta!=0){
        distFromRoot <- node.depth.edgelength(phy)
        phy$edge.length = (exp(beta*distFromRoot[descendent])-exp(beta*distFromRoot[parent]))/beta
    }
    
    if(sigmasq!=0) phy$edge.length <- phy$edge.length*sigmasq
    
    return(phy)
}
