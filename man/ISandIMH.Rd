\name{Importance Sampling and independence Metropolis Hastings sampling}
\alias{IS}
\alias{IMH}
\title{
  Monte Carlo sampling using the iterated Laplace approximation.
}
\description{
  Use iterated Laplace approximation as a proposal for importance
  sampling or the independence Metropolis Hastings algorithm.
}
\usage{
IS(obj, nSim, df = 4, post, vectorized = FALSE, cores = 1, ...)

IMH(obj, nSim, df = 4, post, vectorized = FALSE, cores = 1, ...)
}
\arguments{
  \item{obj}{
    an object of class "mixDist"
}
  \item{nSim}{
    number of simulations
}
  \item{df}{
    degrees of freedom of the mixture of t distributions proposal
}
  \item{post}{
    log-posterior density
  }
  
  \item{vectorized}{
        Logical determining, whether \code{post} is vectorized
      }
      \item{cores}{
	number of cores you want to use to evaluate the target density
        (uses the mclapply function from the parallel package). For
        Windows machines, a value > 1 will have no effect, see mclapply
        help for details.
}
  \item{\dots}{
    additional arguments passed to \code{post}.
}
}
\value{
  A list with entries:\cr
  \code{samp}: Matrix containing sampled values\cr
  \code{w}: Vector of weights for values in samp\cr
  \code{normconst}: normalization constant estimated based on importance
  sampling\cr
  \code{ESS}: Effective sample size (for IS)\cr
  \code{accept}: Acceptance rate (for IMH)\cr
}
\author{
  Bjoern Bornkamp
}

\examples{
  ## see function iterLap for an example on how to use IS and IMH
}
\keyword{ misc }

