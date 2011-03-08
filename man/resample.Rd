\name{resample}
\alias{resample}
\title{
  Residual resampling
}
\description{
  Perform residual resampling to the result of importance sampling
}
\usage{
resample(n, obj)
}
\arguments{
  \item{n}{
    Number of resamples to draw
}
  \item{obj}{
    An object of class \code{IS}, as produced by the \code{IS} function
}
}
\value{
  Matrix with resampled values
}
\author{
  Bjoern Bornkamp
}
\examples{
  ## see function iterLap for an example on how to use resample
}
\keyword{ misc }