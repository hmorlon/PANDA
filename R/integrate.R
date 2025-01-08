.Integrate <- function(...)
{
  op <- getOption("show.error.messages")
  on.exit( options(show.error.messages=op) )
  options(show.error.messages=FALSE)
  
  res <- try((integrate(...)));

  if (inherits(res, "try-error"))
  {
    return(-Inf)
  }
  else
  {
    return(res$value)
  }
}
