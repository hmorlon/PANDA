mantel <-
function (formula = formula(data), data = sys.parent(), nperm = 1000, correlation = c("Pearson", "Spearman", "Kendall")) {
  
  m <- match.call(expand.dots = FALSE)
  m2 <- match(c("formula", "data"), names(m), nomatch = 0)
  m <- m[c(1, m2)]
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  m <- as.matrix(m)
  n <- (1 + sqrt(1 + 8 * nrow(m)))/2
  if (abs(n - round(n)) > 1e-12) 
    stop("Matrix not square.\n")
  n <- round(n)
  if (ncol(m) < 2) 
    stop("Not enough data. \n")
  
  ymat <- as.vector(m[, 1])
  xmat <- as.vector(m[, 2])
     
  if (correlation %in% c("Spearman", "Kendall")){
      ymat <- rank(ymat)
      xmat <- rank(xmat)
  }

  ycor <- ymat
  xcor <- xmat
  
  if (correlation=="Pearson"){mantelr <- cor(xcor, ycor, use = "everything", method = "pearson")}
  if (correlation=="Spearman"){mantelr <- cor(xcor, ycor, use = "everything", method = "spearman")}
  if (correlation=="Kendall"){mantelr <- cor(xcor, ycor, use = "everything", method = "kendall")}
  
  xmat <- full(xmat)
  ymat <- full(ymat)
  xmat <- xmat[col(xmat) > row(xmat)]
  ymat <- ymat[col(ymat) > row(ymat)]
  if (nperm > 0) {
    zstats <- numeric(nperm)
    tmat <- matrix(0, n, n)
    rarray <- rep(0, n)
    ncor <- length(xmat)
    w1 <- sum(xmat)/ncor
    w2 <- sum(xmat^2)
    w2 <- sqrt(w2/ncor - w1^2)
    xmat <- (xmat - w1)/w2
    w1 <- sum(ymat)/ncor
    w2 <- sum(ymat^2)
    w2 <- sqrt(w2/ncor - w1^2)
    ymat <- (ymat - w1)/w2
    
    if (correlation %in% c("Pearson", "Spearman")){  # sum of the cross products
      
      print("delete 2")
      system.time(1)
      
      cresults <- .C("permute", as.double(xmat), as.double(ymat), 
                     as.integer(n), as.integer(length(xmat)), as.integer(nperm), 
                     zstats = as.double(zstats), as.double(as.vector(tmat)), 
                     as.integer(rarray),
                     PACKAGE = "RPANDA")
      
      print("delete 3")
      system.time(1)
      
    }
      
    if (correlation=="Kendall"){
      cresults <- .C("permuteKendall", as.double(xmat), as.double(ymat),
                     as.integer(n), as.integer(length(xmat)), as.integer(nperm), 
                     zstats = as.double(zstats), as.double(as.vector(tmat)), 
                     as.integer(rarray),
                     PACKAGE = "RPANDA")
    }
      
    zstats <- cresults$zstats
    pval1 <- length(which(zstats >= zstats[1]))/nperm
    pval2 <- length(which(zstats <= zstats[1]))/nperm
    pval3 <- length(which(abs(zstats) >= abs(zstats[1])))/nperm
  } else {
    pval1 <- 0
    pval2 <- 0
    pval3 <- 0
  }

  c(mantelr = mantelr, pval1 = pval1, pval2 = pval2, pval3 = pval3)
}
