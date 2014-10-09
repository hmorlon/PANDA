print.fit.bd <- function(x,...)
{
  if (!inherits(x, "fit.bd"))
      stop("object \"x\" is not of class \"fit.bd\"")
  print("\n\tFit Birth Death Model\n\t---------------------\n\n")
  print("\tModel validation\n\t---------------- \n\n")
  cat("\tLog-likelihood:", x$LH, "\n")
  cat("\tAICC:", x$aicc, "\n\n")
  print("\tParameters estimate \n\t-------------------\n")
  cat("    Birth:", x$lamb_par, "\n\n")
  if(!is.null(x$mu_par))
    cat("    Death:", x$mu_par, "\n\n")
}
