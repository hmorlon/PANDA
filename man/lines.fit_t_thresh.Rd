\name{lines.threshML}
\alias{lines.threshML}

\title{
Add to a plot line segments joining the evolutionary rate through time estimated by the "ClimExp" model in the \emph{fit_t_thresh} function.
}
\description{
Plot the estimated evolutionary rate of the latent variable in the "environmental" threshold model as a function of the environmental data and time.
}
\usage{

\method{lines}{threshML}(x, steps = 100, ...)

}
\arguments{
  \item{x}{
an object of class 'threshML' obtained from a \emph{fit_t_thresh} fit.
}
  \item{steps}{
the number of steps from the root to the present used to compute the relative rate \eqn{\sigma(t)} through time.
}
  \item{...}{
further arguments to be passed to \code{lines}. See ?\code{lines}.
}

}


\value{
	\code{lines.threshML} returns invisibly a list with the following components used to add the line segments to the current plot:
  
  \item{time_steps}{the times steps where the climatic function was evaluated to compute the rate. The number of steps is controlled through the argument \code{steps}.}
  \item{values}{the estimated (relative) rate values through time estimated at each \code{time_steps}}
  
}


\note{
	All the graphical parameters (see \code{\link{par}}) can be passed through (e.g. line type: \code{lty}, line width: \code{lwd}, color: \code{col} ...)
}
\references{
Clavel, J. & Morlon, H., 2017. Accelerated body size evolution during cold climatic periods in the Cenozoic. Proceedings of the National Academy of Sciences,  114(16): 4183-4188.

Melendez-Vazquez, F., Lucaci, A. G., Selberg A., Clavel, J., Rincon-Sandoval, M., Santaquiteria, A., White, W. T., Drabeck, D., Carnevale, G., Duarte-Ribeiro, E., Miya, M., Westneat, M. W., Baldwin, C. C., Hughes, L. C., Ortí, G., Kosakovsky Pond, S. L., Betancur-R, R., Arcila, D., 2025. Ecological interactions and genomic innovation fueled the evolution of ray-finned fish endothermy. Sci. Adv., 11 (eads8488):1-16.

}

\author{J. Clavel}

\seealso{
%% add later  \code{\link{simulate_t_env}}  
\code{\link{lines.threshML}}, \code{\link{fit_t_thresh}}
}

\examples{
\dontshow{test = FALSE}
if(test){
data(InfTemp)

\donttest{

## Comparison between models
## require(phytools)
data(InfTemp)
set.seed(1)
tree <- pbtree(n=150, scale=67)
tree <- reorder(tree, "postorder")

beta = -0.2 # strengh and direction of the climatic effect
# The (latent) trait is simulated using the function "sim_t_env" (see also ?fit_t_env)
x <- sim_t_env(tree, param=c(0.1, beta), env_data=InfTemp, model="EnvExp", root.value=0,
               step=0.001, plot=TRUE)

# Use the "threshState" function to simulate the discrete traits from the continuous one 'x'.
th <- sapply(x, threshState, thresholds = setNames(c(0,Inf), letters[1:2]))


# Fit the climatic model
clim_fit1 <- fit_t_thresh(tree, th, model="Clim", env_data=InfTemp, df=50)
plot(clim_fit1, lty=2)

clim_fit2 <- fit_t_thresh(tree, th, model="Clim", env_data=InfTemp, df=10)
lines(clim_fit2, col="red")


}
}

}
