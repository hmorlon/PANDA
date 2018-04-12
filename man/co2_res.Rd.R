\name{co2_res}
\alias{co2_res}
\docType{data}
\title{co2 data since the beginning of the Cenozoic}
\description{Atmospheric co2 data since the beginning of the Cenozoic}
\usage{data(co2_res)}
\details{
  Implied co2 data since the beginning of the Cenozoic taken from Hansen et al., (2013). The data are the amount of co2 in ppm reuquired to yield observed global temperature throughout the Cenozoic:
  \describe{
    \item{\code{age}}{a numeric vector corresponding to the geological age, in Myrs before the present}
    \item{\code{co2}}{a numeric vector corresponding to the estimated co2 at that age}
  }
}
\references{
Hansen, J., Sato, M., Russell, G., Kharecha, P. (2013) Climate sensitivity, sea level and atmospheric carbon dioxide \emph{Philosophoical Transactions of the Royal Society A \emph{371:2001}
}
\examples{
data(co2_res)
plot(co2_res)
}
\keyword{datasets}