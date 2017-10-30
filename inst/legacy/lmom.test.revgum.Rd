\name{lmom.test.revgum}
\alias{lmom.test.revgum}
\title{ Test L-moment and Parameter Algorithms of the Reverse Gumbel Distribution }
\description{
This function computes the L-moments of the data and the parameters of the
Reverse Gumbel distribution and in turn computes the L-moments from the fitted 
parameters.
}
\usage{
lmom.test.revgum(data, zeta=1, digits=4)
}
\arguments{
  \item{data}{A vector of data.}
  \item{zeta}{The right-tail censoring fraction (see \code{\link{parrevgum}}).}
  \item{digits}{Number of digits to pass internally to the \code{signif} function for rounding of results.}
}
\value{
  Comparison of the median of the distribution and reverse computation
  of the median from the 0.5 nonexceedance probability.

  Output from \code{\link{lmom.diff}}.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lmomrevgum}}, \code{\link{parrevgum}} }
\examples{
lmom.test.revgum(c(123,34,4,654,37,78))

lmom.test.revgum(rnorm(50))
}
\keyword{distribution}
