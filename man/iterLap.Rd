\name{iterLap}
\alias{iterLap}
\title{
  Iterated Laplace Approximation
}
\description{
  Iterated Laplace Approximation
}
\usage{
iterLap(post, ..., GRobj = NULL, vectorized = FALSE, startVals = NULL,
        method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"), control = NULL,
        nlcontrol = list())
}
\arguments{
  \item{post}{
    log-posterior density
}
  \item{\dots}{
    additional arguments to log-posterior density
}
  \item{GRobj}{
    object of class mixDist, for example resulting from a call to GRApprox
}
  \item{vectorized}{
    Logical determining, whether \code{post} is vectorized    
}
  \item{startVals}{
    Starting values for GRApprox, when GRobj is not specified.
    Vector of starting values if dimension=1
    otherwise matrix of starting values with the starting values in the rows
  }
  \item{method}{
    Type of optimizer to be used.
  }
  \item{control}{
    List with entries:\cr
    \code{gridSize} Determines the size of the grid for each component\cr
    \code{delta} Stopping criterion based on the maximum error on the grid\cr
    \code{maxDim} Maximum number of components allowed (default 20)\cr
    \code{eps} Stopping criterion based on normalization
    constant of approximation\cr
    \code{info} How much information should be displayed during
    iterations: 0 - none, 1 - minimum information, 2 - maximum information\cr
}
  \item{nlcontrol}{
    Control list for the used optimizer.
}
}
\value{
  Produces an object of class mixDist: A list with entries\cr
  \code{weights} Vector of weights for individual components\cr
  \code{means} Matrix of component medians of components\cr
  \code{sigmas} List containing scaling matrices\cr
  \code{eigenHess} List containing eigen decompositions of scaling matrices\cr
  \code{dets} Vector of determinants of scaling matrix\cr
  \code{sigmainv} List containing inverse scaling matrices\cr
}
\author{
  Bjoern Bornkamp
}
\references{
  Bornkamp, B. (2011). Approximating Probability Densities by Iterated
  Laplace Approximations, \emph{Journal of Computational and Graphical
  Statistics}, \emph{to appear}
}

\examples{
  #### banana example
  banana <- function(pars, b, sigma12){
    dim <- 10
    y <- c(pars[1], pars[2]+b*(pars[1]^2-sigma12), pars[3:dim])
    cc <- c(1/sqrt(sigma12), rep(1, dim-1))
    return(-0.5*sum((y*cc)^2))
  }

  ###############################################################
  ## first perform multi mode Laplace approximation
  start <- rbind(rep(0,10),rep(-1.5,10),rep(1.5,10))
  grObj <- GRApprox(banana, start, b = 0.03, sigma12 = 100)
  ## print mixDist object
  grObj
  ## summary method
  summary(grObj)
  ## importance sampling using the obtained mixDist object 
  ## using a mixture of t distributions with 10 degrees of freedom
  isObj <- IS(grObj, nSim=1000, df = 10, post=banana, b = 0.03,
              sigma12 = 100)
  ## effective sample size
  isObj$ESS
  ## independence Metropolis Hastings algorithm
  imObj <- IMH(grObj, nSim=1000, df = 10, post=banana, b = 0.03,
               sigma12 = 100)
  ## acceptance rate
  imObj$accept

  ###############################################################
  ## now use iterated Laplace approximation
  ## and use Laplace approximation above as starting point
  iL <- iterLap(banana, GRobj = grObj, b = 0.03, sigma12 = 100)
  isObj2 <- IS(iL, nSim=10000, df = 100, post=banana, b = 0.03,
               sigma12 = 100)
  ## residual resampling to obtain unweighted sample
  samples <- resample(1000, isObj2)
  ## plot samples in the first two dimensions
  plot(samples[,1], samples[,2], xlim=c(-40,40), ylim = c(-40,20))
  ## independence Metropolis algorithm
  imObj2 <- IMH(iL, nSim=1000, df = 10, post=banana, b = 0.03,
                sigma12 = 100)
  imObj2$accept
  plot(imObj2$samp[,1], imObj2$samp[,2], xlim=c(-40,40), ylim = c(-40,20))

  ## IMH and IS can exploit multiple cores, example for two cores
  \dontrun{
  isObj3 <- IS(iL, nSim=10000, df = 100, post=banana, b = 0.03,
               sigma12 = 100, cores = 2)
  }

}
\keyword{ misc }

