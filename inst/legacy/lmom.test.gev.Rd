\name{lmom.test.gev}
\alias{lmom.test.gev}
\title{ Test L-moment and Parameter Algorithms of the Generalized Extreme Value Distribution}
\description{
This function computes the L-moments of the data and the parameters of the
Generalized Extreme Value distribution and in turn computes the L-moments from the fitted 
parameters.
}
\usage{
lmom.test.gev(data, digits=4)
}
\arguments{
  \item{data}{A vector of data.}
  \item{digits}{Number of digits to pass internally to the \code{signif} function for rounding of results.}
}
\value{
  Comparison of the median of the distribution and reverse computation
  of the median from the 0.5 nonexceedance probability.

  Output from \code{\link{lmom.diff}}.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lmomgev}}, \code{\link{pargev}} }
\examples{
lmom.test.gev(c(123,34,4,654,37,78))

lmom.test.gev(rnorm(50))
}
\keyword{distribution}
