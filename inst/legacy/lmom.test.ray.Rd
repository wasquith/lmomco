\name{lmom.test.ray}
\alias{lmom.test.ray}
\title{ Test L-moment and Parameter Algorithms of the Rayleigh Distribution }
\description{
This function computes the L-moments of the data and the parameters of the
Rayleigh distribution and in turn computes the L-moments from the fitted 
parameters.
}
\usage{
lmom.test.ray(data, digits=4)
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
\seealso{ \code{\link{lmomray}}, \code{\link{parray}} }
\examples{
lmom.test.ray(c(123,34,4,654,37,78))

lmom.test.ray(rnorm(50))
}
\keyword{distribution}
