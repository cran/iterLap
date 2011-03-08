\name{GRApprox}
\alias{GRApprox}
\title{
  Gelman-Rubin mode approximation
}
\description{
  Performs the multiple mode approximation of Gelman-Rubin (applies a
  Laplace approximation to each mode). The weights are determined
  corresponding to the height of each mode.
}
\usage{
GRApprox(post, start, grad, method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
         control = list(), ...)
}
\arguments{
  \item{post}{
    log-posterior density.
}
  \item{start}{
    vector of starting values if dimension=1
    otherwise matrix of starting values with the starting values in the rows
  }
  \item{grad}{
    gradient of log-posterior
  }
  \item{method}{
    Which optimizer to use
}
  \item{control}{
    Control list for the chosen optimizer
}
  \item{\dots}{
    Additional arguments for log-posterior density specified in \code{post}
}
}
\value{
  Produces an object of class mixDist. That a list mit entries\cr
  \code{weights} Vector of weights for individual components\cr
  \code{means} Matrix of component medians of components\cr
  \code{sigmas} List containing scaling matrices\cr
  \code{eigenHess} List containing eigen decompositions of scaling matrices\cr
  \code{dets} Vector of determinants of scaling matrix\cr
  \code{sigmainv} List containing inverse scaling matrices\cr
}
\references{
  Gelman, A., Carlin, J. B., Stern, H. S. & Rubin, D. B. (2003) Bayesian Data
  Analysis, 2nd edition, Chapman and Hall. (Chapter 12)

  Bornkamp, B. (2011). Approximating Probability Densities by Iterated
  Laplace Approximations, \emph{Journal of Computational and Graphical
  Statistics}, \emph{to appear}
}
\author{
  Bjoern Bornkamp
}

\seealso{
 \code{\link{iterLap}}
}
\examples{

  ## log-density for banana example
  banana <- function(pars, b, sigma12){
    dim <- 10
    y <- c(pars[1], pars[2]+b*(pars[1]^2-sigma12), pars[3:dim])
    cc <- c(1/sqrt(sigma12), rep(1, dim-1))
    return(-0.5*sum((y*cc)^2))
  }

  start <- rbind(rep(0,10),rep(-1.5,10),rep(1.5,10))
  ## multiple mode Laplace approximation
  aa <- GRApprox(banana, start, b = 0.03, sigma12 = 100)
  ## print mixDist object
  aa
  ## summary method
  summary(aa)
  ## importance sampling using the obtained mixDist object 
  ## using a mixture of t distributions with 10 degrees of freedom
  dd <- IS(aa, nSim=1000, df = 10, post=banana, b = 0.03,
           sigma12 = 100)
  ## effective sample size
  dd$ESS
}
\keyword{ misc }
