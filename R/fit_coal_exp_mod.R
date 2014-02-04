fit_coal_exp_mod <- function (phylo, tau0=1.e-2, gamma=1, meth = "Nelder-Mead", N0=0, Vtimes=FALSE)
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  if (Vtimes==TRUE)
  {
    Vtimes<-sort(phylo)
    ntips<-length(phylo)+1
    if (N0==0)
      {
        N0<-ntips
      }
  }

  else
  {
    Vtimes <- sort(branching.times(phylo))
    ntips<-Ntip(phylo)
    if (N0==0)
    {
      N0<-ntips
    }
  }

  init<-c(tau0,gamma)
  nbpar<-length(init)
  nbobs<-length(Vtimes)-1

  optimLH.MoranEXP <- function(init)
  {
    tau0 <- init[1]
    gamma <- init[2]
    LH <- likelihood_coal_exp_mod(Vtimes,ntips,tau0,gamma,N0)$res
    return(-LH)
  }

  temp<-suppressWarnings(optim(init, optimLH.MoranEXP, method = meth))

  res <- list(model = "MoranEXP", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), tau0 = temp$par[1], gamma=temp$par[2])
  return(res)
}
