.Integrate <- function(...)
{
  op <- getOption("show.error.messages")
  options(show.error.messages=FALSE)
  res <- try((integrate(...)));
  options(show.error.messages=op)
  if (inherits(res, "try-error"))
  {
    return(-Inf)
  }
  else
  {
    return(res$value)
  }
}
