\name{lmom.test.gum}
\alias{lmom.test.gum}
\title{ Test L-moment and Parameter Algorithms of the Gumbel Distribution }
\description{
This function computes the L-moments of the data and the parameters of the
Gumbel distribution and in turn computes the L-moments from the fitted 
parameters.
}
\usage{
lmom.test.gum(data, digits=4)
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
\seealso{ \code{\link{lmomgum}}, \code{\link{pargum}} }
\examples{
lmom.test.gum(c(123,34,4,654,37,78))

lmom.test.gum(rnorm(50))
}
\keyword{distribution}
