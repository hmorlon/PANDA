################################################################################
##                                                                            ##
##                                RPANDA : Utils                              ##
##                                                                            ##
##   Julien Clavel - 01-02-2018                                               ##
##   S3 methods, simulations, miscellaneous                                   ##
##                                                                            ##
################################################################################

# S3 generic method "ancestral" for reconstructing or retrieving ancestral states (see phyl.pca_pl.R) # should I use "predict" instead?
#ancestral <- function(object) UseMethod("ancestral")

# S3 for the fit_t_env class
ancestral.fit_t.env <- function(object, ...){
    
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
ancestral.fit_t.comp <- function(object, ...){
    
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
# From the mvtnorm package 
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

## for the clade-shift model 
## N. Mazet
extract.clade.ln<-function (phy, node, root.edge = 0) 
{
  Ntip <- length(phy$tip.label)
  ROOT <- Ntip + 1
  Nedge <- dim(phy$edge)[1]
  wbl <- !is.null(phy$edge.length)
  if (length(node) > 1) {
    node <- node[1]
    warning("only the first value of 'node' has been considered")
  }
  if (is.character(node)) {
    if (is.null(phy$node.label)) 
      stop("the tree has no node labels")
    node <- which(phy$node.label %in% node) + Ntip
  }
  if (node <= Ntip) 
    stop("node number must be greater than the number of tips")
  if (node == ROOT) 
    return(phy)
  
  
  keep_nodes<-c(node)
  
  while(1)
  {
    keep_nodes1<-unique(c(keep_nodes,phy$edge[which(phy$edge[,1] %in% keep_nodes),2]))
    if (length(keep_nodes1)==length(keep_nodes))
      break
    else
      keep_nodes<-keep_nodes1
    
  }	
  
  #print(keep_nodes)
  
  keep<-which(phy$edge[,1] %in% keep_nodes)
  
  phy$edge <- phy$edge[keep, ]
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep]
  TIPS <- phy$edge[, 2] <= Ntip
  tip <- phy$edge[TIPS, 2]
  phy$tip.label <- phy$tip.label[tip]
  phy$edge[TIPS, 2] <- order(tip)
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[sort(unique(phy$edge[, 
                                                          1])) - Ntip]
  Ntip <- length(phy$tip.label)
  phy$Nnode <- dim(phy$edge)[1] - Ntip + 1L
  newNb <- integer(Ntip + phy$Nnode)
  newNb[node] <- Ntip + 1L
  sndcol <- phy$edge[, 2] > Ntip
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (Ntip + 
                                                          2):(Ntip + phy$Nnode)
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  phy
}

## core analysis of backbone that is parallelized.
## N. Mazet

all_comb_models <- function(to){
  
  # splitting combination into subclades and backbones
  
  comb1 <- strsplit(comb.shift[to], "/")[[1]]
  comb.sub <- strsplit(comb1[[1]], "[.]")[[1]]
  if(length(comb1) == 2){
    comb.bck <- strsplit(comb.shift[to], "/")[[1]][2]
    comb.bck <- strsplit(comb.bck, "[.]")[[1]]
  } else {
    comb.bck <- NULL
  }
  
  message("\n", to, "/", length(comb.shift))
  
  # plot to illustrate
  # plot.phylo.comb(phylo, data, sampling.fractions, comb = comb.shift[to], cex = 0.8, label.offset = 0.2)
  # nodes <- sampling.fractions$nodes[!is.na(sampling.fractions$to_test)]
  # nodelabels(as.character(nodes), nodes)
  
  # create the backbone
  # new way 
  
  int_nodes <- comb.bck
  # order from present to past
  int_nodes <- names(branching.times(phylo)[order(branching.times(phylo))])[names(branching.times(phylo)[order(branching.times(phylo))]) %in% int_nodes]
  branch_times_to_bck <- rep(list(NULL), length(comb.bck)+1)
  phylo_backbone_cut <- rep(list(NULL), length(comb.bck)+1)
  phylo_backbone_core <- drop.tip(phylo, unlist(ALL_clade_names[comb.sub]))
  
  res_bck <- rep(list(NULL), length(comb.bck)+1)
  
  sb.tips <- rep(list(NULL), length(int_nodes))
  sb.desc <- rep(list(NULL), length(int_nodes))
  names(sb.desc) <- int_nodes
  
  for(sb in 1:length(res_bck)){
    
    if(is.null(comb.bck)){ # simple backbone
      
      phylo_backbone_cut <- list(phylo_backbone_core)
      names(phylo_backbone_cut) <- paste0(paste0(comb.sub, collapse = "."),"_bck")
      
      branch_time_sb <- get.branching.nodes(comb.sub, phylo = phylo,
                                            ALL_branch_times_clades = ALL_branch_times_clades,
                                            ALL_clade_names = ALL_clade_names)
      branch_times_to_bck <- list(branch_time_sb)
      names(branch_times_to_bck) <- paste0(comb.sub, collapse = ".")
      
      # check the root? seems ok with parnassiinae
      
    } else { # multibackbone
      
      if(sb < length(phylo_backbone_cut)){ # before deep backbone
        
        sb.desc[[sb]] <- Descendants(phylo, as.numeric(int_nodes[sb]), "all")
        
        if(sb > 1){ # removing descendant in previous int_nodes
          sb.desc[[sb]] <- sb.desc[[sb]][!sb.desc[[sb]] %in% unlist(sb.desc[1:c(sb-1)])]
        }
        
        sb.desc_sb_sp <- phylo$tip.label[sb.desc[[sb]][sb.desc[[sb]] < Ntip(phylo)]]
        sb.desc_sb_sp <- intersect(sb.desc_sb_sp, phylo_backbone_core$tip.label)
        phylo_backbone_cut[[sb]] <- subtree(phylo_backbone_core, sb.desc_sb_sp)
        names(phylo_backbone_cut)[sb] <- paste0(int_nodes[sb],"_sub")
        
        comb.multibackbone <- c(comb.sub[comb.sub %in% sb.desc[[sb]]], int_nodes[int_nodes %in% sb.desc[[sb]]])
        
        branch_time_sb <- get.branching.nodes(comb.multibackbone, phylo = phylo,
                                              ALL_branch_times_clades = ALL_branch_times_clades,
                                              ALL_clade_names = ALL_clade_names)
        
        # check that root of phylo_backbone_cut[[sb]] is int_node
        if(phylo_backbone_cut[[sb]]$node.label[1] != int_nodes[sb] &
           !phylo_backbone_cut[[sb]]$node.label[1] %in% names(branch_time_sb)){
          
          root_sb_to_int_nodes <- c(phylo_backbone_cut[[sb]]$node.label[1], Ancestors(phylo, phylo_backbone_cut[[sb]]$node.label[1]))
          root_sb_to_int_nodes <- root_sb_to_int_nodes[1:c(which(root_sb_to_int_nodes == int_nodes[sb])-1)]
          missed_sb_nodes <- root_sb_to_int_nodes[!root_sb_to_int_nodes %in% as.numeric(names(branch_time_sb))]
          
          for(msb in 1:length(missed_sb_nodes)){
            branch_time_missing_sb <- list(c(missed_sb_nodes[msb], Ancestors(phylo, missed_sb_nodes[msb], "parent")))
            names(branch_time_missing_sb) <- missed_sb_nodes[msb]
            
            branch_time_sb[length(branch_time_sb)+1] <- branch_time_missing_sb
            names(branch_time_sb)[length(branch_time_sb)]<- as.character(missed_sb_nodes[msb])
          }
        }
        
        branch_times_to_bck[sb] <- list(branch_time_sb)
        names(branch_times_to_bck)[sb] <- paste(comb.multibackbone, collapse = ".")
        
      } else {  # deep backbone
        
        tips_up_bck <- unlist(lapply(phylo_backbone_cut, function(x) x$tip.label))
        # remaining comb.sub in the deep backbone
        tips_last_bck <- unlist(ALL_clade_names[comb.sub[!comb.sub %in% unlist(sb.desc, use.names =FALSE)]])
        
        phylo_backbone_cut[[sb]] <- drop.tip(phylo_backbone_core, tips_up_bck)
        names(phylo_backbone_cut)[sb] <- paste(int_nodes[sb-1],"bck", sep = "_")
        
        int_nodes_deep_backbone <- int_nodes[!int_nodes %in% unlist(sapply(branch_times_to_bck, names), use.names =FALSE)]
        
        comb_deep_backbone <- c(comb.sub[!comb.sub %in% unlist(sb.desc, use.names =FALSE)], int_nodes_deep_backbone)
        
        branch_time_sb <- get.branching.nodes(comb_deep_backbone, phylo = phylo,
                                              ALL_branch_times_clades = ALL_branch_times_clades,
                                              ALL_clade_names = ALL_clade_names)
        
        branch_times_to_bck[sb] <- list(branch_time_sb)
        names(branch_times_to_bck)[sb] <- paste(comb_deep_backbone, collapse = ".")
        
      } # deep backbone
    } # multi backbone
  }
  
  branch_nodes_to_bck <- branch_times_to_bck
  for(bck in 1:length(branch_times_to_bck)){
    for(nodeID in 1:length(branch_nodes_to_bck[[bck]])){
      branch_times_to_bck[[bck]][[nodeID]] <- sapply(branch_nodes_to_bck[[bck]][[nodeID]], get.node.ages, phylo = phylo)
    }  
  }
  
  # Sampling fractions ####
  
  lin.node <- data.frame(node = c(comb.sub,comb.bck, Ntip(phylo)+1), n.tips = rep(NA, length(comb.sub) + length(comb.bck)+1))
  lin.node$node <- as.character(lin.node$node)
  lin.node <- merge(lin.node, sampling.fractions[sampling.fractions$nodes %in% lin.node$node, c("nodes", "sp_tt"),],
                    by.x = "node", by.y = "nodes")
  
  node_order <- names(branching.times(phylo)[order(branching.times(phylo))])
  node_order <- node_order[node_order %in% lin.node$node]
  
  lin.node <- lin.node[match(node_order, lin.node$node),]
  
  for(n.lin in 1:nrow(lin.node)){
    desc.n.lin <- length(Descendants(phylo, as.numeric(lin.node$node[n.lin]))[[1]])
    # whether this node is present in an other lineage
    int.n.lin <- Descendants(phylo, as.numeric(lin.node$node[n.lin]), type = "all")
    int.n.lin <- as.character(int.n.lin[int.n.lin > Ntip(phylo)])
    # Ntip
    if(any(comb.sub %in% int.n.lin)){
      lin.node$n.tips[n.lin] <- desc.n.lin - sum(lin.node$n.tips[lin.node$node %in% comb.sub[comb.sub %in% int.n.lin]])
      lin.node$sp_tt[n.lin] <- lin.node$sp_tt[n.lin] - sum(lin.node$sp_tt[lin.node$node %in% comb.sub[comb.sub %in% int.n.lin]])
    } else{
      lin.node$n.tips[n.lin] <- desc.n.lin
    }
  }
  
  lin.node$n.tips_prev <- lin.node$n.tips
  lin.node$sp_tt_prev <- lin.node$sp_tt
  
  lin.node_bck <- lin.node[!lin.node$node %in% comb.sub,]
  
  for(l.n in c(1:nrow(lin.node_bck))){
    int.desc_lin <- unlist(Descendants(phylo, as.numeric(lin.node_bck$node[l.n]), "all"))
    int.desc_lin <- int.desc_lin[int.desc_lin > Ntip(phylo)]
    
    if(any(lin.node_bck$node %in% int.desc_lin)){
      
      bck_up <- lin.node_bck[which(lin.node_bck$node %in% int.desc_lin),]
      
      ntip_bck_up <- sum(bck_up$n.tips_prev)
      ntaxo_bck_up <- sum(bck_up$sp_tt_prev)
      
      lin.node_bck$n.tips_prev[l.n] <- lin.node_bck$n.tips[l.n] - ntip_bck_up
      lin.node_bck$sp_tt_prev[l.n] <- lin.node_bck$sp_tt[l.n] - ntaxo_bck_up  
      
    }
  }
  lin.node[lin.node$node %in% lin.node_bck$node,] <- lin.node_bck
  
  lin.node <- lin.node[-(1:length(comb.sub)),]
  
  f <- as.list(lin.node$n.tips_prev/lin.node$sp_tt_prev)
  names(f) <- names(phylo_backbone_cut)
  
  for(btb in 1:length(phylo_backbone_cut)){
    
    # by default backbone.option = "crown.shift"
    backbone <- backbone.option
    spec_times <- NULL
    cond <- "crown"
    
    # CHECKED!
    tot_time3 <- max(c(node.age(phylo_backbone_cut[[btb]])$ages, unlist(branch_times_to_bck[[btb]])))
    
    # for converting in stem.shift
    if(backbone.option == "stem.shift"){
      
      spec_times <- sapply(branch_times_to_bck[[btb]], "[[", 2)
      cond <- "stem"
      
      if(!is.null(phylo_backbone_cut[[btb]]$root.edge)){
        tot_time3 <- max(node.age(phylo_backbone_cut[[btb]])$ages) + phylo_backbone_cut[[btb]]$root.edge
      }
      
      # if deep backbone, conditioning backbone at crown 
      if(length(grep("_bck", names(phylo_backbone_cut[btb]))) == 1){
        cond <- "crown"
      }
      branch_times_to_bck[[btb]] <- rep(list(NULL),1)
    }
    
    ##################################### models
    
    results <- div.models(phylo = phylo_backbone_cut[[btb]], tot_time = tot_time3, f = f[[btb]],
                          backbone = backbone, spec_times = spec_times, branch_times = branch_times_to_bck[[btb]],
                          cond = cond, models = models, n.max = n.max, rate.max = rate.max, verbose = TRUE)
    if(btb < length(phylo_backbone_cut)){
      # cond has to be changed to properly estimate likelihood of each part if they are not the last part
      results1 <- div.models(phylo = phylo_backbone_cut[[btb]], tot_time = tot_time3, f = f[[btb]],
                             backbone = backbone, spec_times = spec_times, branch_times = branch_times_to_bck[[btb]],
                             cond =FALSE, models = models, n.max = n.max, rate.max = rate.max, verbose =FALSE)
      
      results2 <- merge(results1[,c(1:4)], results[,c(1,5:8)], by="Models")
      results <- results2[match(results$Models, results2$Models),]
      
      # adding a parameter for the location of the shift (to modify for the printing)
      results$AICc <- 2 * -results$logL + 2 * (results$Parameters+1) + (2 * (results$Parameters+1) * ((results$Parameters+1) + 1))/(Ntip(phylo_backbone_cut[[btb]]) - (results$Parameters+1) - 1)
      results$Parameters <- results$Parameters+1
    }
    
    results[,-1] <- apply(results[,-1], 2, as.numeric)
    res_bck[btb] <- list(results)
  }
  
  desc_comb.sub <-  Descendants(phylo, as.numeric(comb.sub), "all")
  desc_comb.sub <- lapply(desc_comb.sub, function(x) x[x > Ntip(phylo)])
  
  nodes_backbone_th <- setdiff(phylo$node.label, unlist(desc_comb.sub))
  
  nodes_backbone_obs <- unlist(lapply(phylo_backbone_cut, function(x) x$node.label), use.names =FALSE)
  all_branching_nodes_to <- unlist(lapply(branch_nodes_to_bck, function(x) unique(sapply(x, "[[", 2))), use.names =FALSE)
  branch_nodes_to_bck <- unlist(lapply(branch_nodes_to_bck, names), use.names =FALSE)
  
  nodes_backbone_obs <- as.numeric(c(nodes_backbone_obs,
                                     branch_nodes_to_bck,
                                     all_branching_nodes_to))
  
  if(all(nodes_backbone_th %in% unique(nodes_backbone_obs))){
    check <- T
  }
  
  names(res_bck) <- names(phylo_backbone_cut)
  
  if(!check){
    stop("\n#### Some branches are missing... ####\n")
  }
  return(res_bck)
  # Multi merge
}

## Get node ages for the backbone in the clade-shift model.
## N. Mazet
get.node.ages <-  function(nodes, ...){
  
  dots <- list(...)
  if(!hasArg(phylo)) stop()
  phylo <- dots$phylo
  
  ALL_nodes_ages <- as.data.frame(apply(data.frame(nodesID=names(branching.times(phylo)),ages=branching.times(phylo)), 2, as.numeric))
  nodes_ages_selected <- sort(ALL_nodes_ages$ages[ALL_nodes_ages$nodesID %in% nodes])
  return(nodes_ages_selected)
}

## Get branching nodes from a combination for the clade-shift model.
## N. Mazet

get.branching.nodes <- function(comb, ...){
  
  dots <- list(...)
  if(!hasArg(phylo)) stop()
  phylo <- dots$phylo
  if(!hasArg(ALL_branch_times_clades)) stop()
  ALL_branch_times_clades <- dots$ALL_branch_times_clades
  if(!hasArg(ALL_clade_names)) stop()
  ALL_clade_names <- dots$ALL_clade_names
  
  root_ID = phylo$node.label[1]
  
  root_clade <- 0
  root_node <-  NULL
  
  # account for poor backbone resulting in a subclade
  phylo_backbone_sb <- drop.tip(phylo, unlist(ALL_clade_names[comb]))
  sibling_shift_nodes <- unlist(Siblings(phylo, as.numeric(comb)))
  
  shift <- ALL_branch_times_clades[comb]
  
  if(phylo_backbone_sb$node.label[1] != phylo$node.label[1]){
    
    root_clade_sb <- list(list(c(phylo_backbone_sb$node.label[1],
                                 Ancestors(phylo, phylo_backbone_sb$node.label[1], type = "parent"))))
    names(root_clade_sb) <- phylo_backbone_sb$node.label[1]
    shift <- c(shift, root_clade_sb)
  }
  
  # coalescence (core of the function)
  df_ALL <- as.data.frame(sapply(unlist(shift,recursive =FALSE), function(m) m[2]))
  colnames(df_ALL) <- "node"
  row.names(df_ALL) <- 1:nrow(df_ALL)
  
  # detect the root in the clades TO REMOVE BECAUSE ONLY ON PARENTAL NODES
  #if(any(df_ALL$node == Ntip(phylo) + 1)){
  # root_clade <- 0
  #root_node <-  NULL # because already in df_all
  #} 
  
  df_ALL <- data.frame(node = df_ALL[which(!df_ALL$node %in% c(root_ID)),])
  
  if(nrow(df_ALL) > 1){
    
    all_ancestors <- unlist(list(rep(list(NULL), nrow(df_ALL))),recursive =FALSE)
    
    for(df_l in 1:nrow(df_ALL)){
      
      all_ancestors[df_l] <- list(c(df_ALL$node[df_l],Ancestors(phylo, df_ALL$node[df_l], type = "all")))
      
    }
    
    # removing root node
    all_ancestors <- lapply(all_ancestors, function(x) x[1:c(which(x == root_ID)-1)])
    
    # counting parental nodes
    ALL_par_nodes <- NULL
    
    coal <- as.data.frame(table(sapply(all_ancestors, function(m) m[1])))
    
    while(any(coal$Freq == 2) & is.null(all_ancestors) ==FALSE){
      
      if(any(coal$Freq == 2)){
        ALL_par_nodes <- c(ALL_par_nodes,as.numeric(as.character(coal$Var1[coal$Freq == 2])))
        all_ancestors <- unique(all_ancestors)
        
        all_ancestors <- lapply(all_ancestors, function(x) x[x %in% ALL_par_nodes == FALSE])
        
        coal <- as.data.frame(table(sapply(all_ancestors, function(m) m[1])))
        
      } else {
        ALL_par_nodes <- NULL
        all_ancestors <- NULL
      }
    }
    
  } else {ALL_par_nodes <- NULL}
  
  
  if(length(ALL_par_nodes) != 0){
    
    parental_nodes <- unlist(list(rep(list(NULL),length(ALL_par_nodes))),recursive = FALSE)
    
    # ALL OTHER NODES
    if(length(parental_nodes) != 0){
      
      for(df_l in 1:length(parental_nodes)){
        
        if(ALL_par_nodes[df_l] != root_ID){
          
          parental_nodes[[df_l]] <- c(ALL_par_nodes[df_l], Ancestors(phylo, ALL_par_nodes[df_l], type = "parent"))
        }
      }
      # WHETHER PARENTAL NODES ARE THE ROOT
      for(p in 1:length(parental_nodes)){
        
        if(parental_nodes[[p]][2] == root_ID){
          root_node <- parental_nodes[[p]][1]
          root_clade <- 1
        }
        
      }
    }
    
  } else {
    parental_nodes <- NULL
  }
  
  df_ALL <- t(as.data.frame(unlist(shift,recursive = FALSE)))
  
  branches_df_all <- apply(df_ALL, 1, paste, collapse = ".")
  if(!is.null(parental_nodes)){
    branches_parental <- apply(do.call(rbind, parental_nodes), 1, paste, collapse = ".")
    parental_nodes <- parental_nodes[!branches_parental %in% branches_df_all]
  }
  
  #
  branch_times_to <- unlist(list(rep(list(NULL),nrow(df_ALL) + length(parental_nodes) + root_clade)),recursive = FALSE)
  
  bt_1 <-  unlist(shift,recursive = FALSE)
  for(bt in 1:length((bt_1))){
    branch_times_to[bt] <- bt_1[bt]
  }
  
  p = 0
  if(length(parental_nodes) != 0){
    for(p in 1:length(parental_nodes)){
      branch_times_to[bt + p] <- parental_nodes[p]
    }
  }
  
  branch_root <- c(Siblings(phylo, root_node),Ancestors(phylo, root_node, type = "parent"))
  
  if(root_clade == 1 & paste(branch_root, collapse = ".") %in% sapply(branch_times_to, paste0, collapse= ".") == FALSE){
    branch_times_to[bt + p + root_clade] <- list(branch_root)
  }
  
  #names(branch_times_to) <- c(names(bt_1),rep("parental_node",length(parental_nodes)),rep("root",length(root_node)))
  names(branch_times_to) <- c(names(bt_1), sapply(parental_nodes, function(x) ifelse(!is.null(x), x[1], NULL)), Siblings(phylo, root_node))
  branch_times_to <- branch_times_to[!sapply(branch_times_to, is.null)]
  
  return(branch_times_to)
}

## Used in clade.shift model to isolate subclade
## N. Mazet

subtree<-function(tree,species_list)
  
{
  # find the MRCA of the species in the list
  node<-unique(tree$edge[which(tree$edge[,2] %in% which(tree$tip.label %in% species_list)),1])
  subtree<-extract.clade.ln(tree,min(node))
  
  while(sum(!(species_list %in% subtree$tip.label))>0)
  {
    node<-unique(tree$edge[which(tree$edge[,2] %in% node),1])
    subtree<-extract.clade.ln(tree,min(node))
    #print(node)
  }
  
  root_length<-tree$edge.length[which(tree$edge[,2] == min(node))]
  subtree$root.edge<-root_length
  
  return(subtree)}


# function from ParallelLogger

#clusterApply <- function(cluster, x, fun, ..., stopOnError = FALSE, progressBar = TRUE) 
