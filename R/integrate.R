.Integrate <- function(...)
{
  res <- try(integrate(...));
  if (class(res) == 'try-error')
  {
    return(-Inf)
  }
  else
  {
    return(res$value)
  }
}
