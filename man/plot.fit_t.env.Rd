\name{plot.fit_t.env}
\alias{plot}
\alias{plot.fit_t.env}

\title{
Plot the phenotypic evolutionary rate through time estimated by the \emph{fit_t_env} function
}
\description{
Fits model of trait evolution for which evolutionary rates depends on an environmental function.
}
\usage{
## S3 method for class 'fit_t.env'
plot(x, steps = 100, ...)

}
\arguments{
  \item{x}{
an object of class 'fit_t.env' obtained from a \emph{fit_t_env} fit.
}
  \item{steps}{
the number of steps from the root to the present used to compute the evolutionary rate \eqn{\sigma2} through time.
}
  \item{...}{
further arguments to be passed to \code{plot}. See ?\code{plot}.
}

}


\value{
	\code{plot.fit_t.env} returns invisibly a list with the following components used in the current plot:
  
  \item{time_steps}{the times steps where the climatic function was evaluated to compute the rate. The number of steps is controlled through the argument \code{steps}.}
  \item{rates}{the estimated evolutionary rate through time estimated at each \code{time_steps}}
  
}


\note{
	All the graphical parameters (see \code{\link{par}}) can be passed through (e.g. line type: \code{lty}, line width: \code{lwd}, color: \code{col} ...)
}
\references{
Clavel, J. & Morlon, H. (in prep) 
}

\author{J. Clavel}

\seealso{
\code{\link{lines.fit_t.env}}, \code{\link{simulate_t_env}}, \code{\link{likelihood_t_env}}
}

\examples{
data(Cetacea)
data(InfTemp)

# Simulate trait on the Cetacean tree
set.seed(123)
trait<-rTraitCont(Cetacea)

## Fit the Environmental-exponential model

result1=fit_t_env(Cetacea, trait, env_data=InfTemp, par=list(scale=TRUE))
plot(result1)

# further options
plot(result1, lty=2, lwd=2, col="red")


}
