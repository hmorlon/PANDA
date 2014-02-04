fit_coal_cst_mod <- function (phylo, tau0=1, meth = "Nelder-Mead", N0=0, Vtimes=FALSE)
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  if (Vtimes==TRUE)
  {
    Vtimes<-sort(phylo)
    ntips<-length(phylo)+1
  }
  else
  {
    Vtimes <- sort(branching.times(phylo))
    ntips<-Ntip(phylo)
  }
  if (N0==0) {N0<-ntips}

  init<-c(tau0)

  nbpar<-length(init)
  nbobs<-length(Vtimes)-1

  optimLH.MoranCST <- function(init)
  {
    tau0 <- init[1]
    LH <- likelihood_coal_cst_mod(Vtimes,ntips,tau0,N0)$res
    return(-LH)
  }

  temp<-suppressWarnings(optim(init, optimLH.MoranCST, method = meth))

  res <- list(model = "MoranCST", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), tau0 = temp$par[1])

  return(res)
}
