\name{iterLap-package}
\alias{iterLap-package}
\docType{package}
\title{
  iterLap package information
}
\description{
  Implementation of iterLap
}
\details{
This package implements the multiple mode Laplace approximation
by Gelman and Rubin (via function GRApprox) and the iterated
Laplace approximation (via the function iterLap). Both functions
return objects of class mixDist, which contain the fitted mode
vectors and covariance matrices. Print and summary methods exist
to display the contents of a mixDist object in human-readable form.
Function IS performs importance sampling, using a mixDist object
as input parameter.

}
\author{
Bjoern Bornkamp

Maintainer: Bjoern Bornkamp <bbnkmp@gmail.com>
}
\references{
  Bornkamp, B. (2011). Approximating Probability Densities by Iterated
  Laplace Approximations, \emph{Journal of Computational and Graphical
  Statistics}, \bold{20}(3), 656--669.
}

\keyword{ package }

\examples{
  ## banana example
  banana <- function(pars, b, sigma12){
    dim <- 10
    y <- c(pars[1], pars[2]+b*(pars[1]^2-sigma12), pars[3:dim])
    cc <- c(1/sqrt(sigma12), rep(1, dim-1))
    return(-0.5*sum((y*cc)^2))
  }

  start <- rbind(rep(0,10),rep(-1.5,10),rep(1.5,10))
  ## multiple mode Laplace approximation
  gr <- GRApprox(banana, start, b = 0.03, sigma12 = 100)
  ## print mixDist object
  gr
  ## summary method
  summary(gr)
  ## importance sampling using the obtained mixDist object 
  ## using a mixture of t distributions with 10 degrees of freedom
  issamp <- IS(gr, nSim=1000, df = 10, post=banana, b = 0.03,
               sigma12 = 100)
  ## effective sample size
  issamp$ESS

  ## now use iterated Laplace approximation (using gr mixDist object
  ## from above as starting approximation)
  iL <- iterLap(banana, GRobj = gr, b = 0.03, sigma12 = 100)
  ISObj <- IS(iL, nSim=10000, df = 100, post=banana, b = 0.03,
              sigma12 = 100)
  ## residual resampling to obtain unweighted sample
  sims <- resample(1000, ISObj)
  plot(sims[,1], sims[,2], xlim=c(-40,40), ylim = c(-40,20))
}
