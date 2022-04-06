\name{remove.model}
\alias{remove.model}

\title{
  Removing a model from shift.estimates output
}
\description{
  Allows to remove a model from the model comparisons of shift.estimates output.
}
\usage{
  remove.model(shift.res, model)
}

\arguments{
  \item{shift.res}{
    the output resulting from shift.estimates.
  }
  \item{model}{
    character. Specifies the model to remove from the set of model of diversification applied in shift.res.
  }
}

\value{
  the same output resulting from shift.estimates but without the chosen model in model comparisons.
  
}
\details{
  This function allow to remove model one at a time. The idea is to remove a model without having to reanalyse the phylogeny and all the combinations of shifts if a model (e.g. BVAR_DVAR) behaves strangely on the studied phylogeny.
  
}
\references{
	Mazet, N., Morlon, H., Fabre, P-H. Condamine, F. (in prep) 
	
}
\author{Nathan Mazet}
\seealso{
\code{\link{shift.estimates}}
}
\examples{
\dontrun{

# To remove "BVAR_DVAR" models if it was included before
# with shifts_cetacea the output from shift.estimates
shifts_cetacea_noBVAR_DVAR <- remove.models(shift.res = shifts_cetacea,
                                            model = "BVAR_DVAR")

}
}