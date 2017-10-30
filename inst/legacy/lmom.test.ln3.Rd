\name{lmom.test.ln3}
\alias{lmom.test.ln3}
\title{ Test L-moment and Parameter Algorithms of the 3-Parameter Log-Normal Distribution }
\description{
This function computes the L-moments of the data and the parameters of the
Log-Normal3 distribution and in turn computes the L-moments from the fitted 
parameters.
}
\usage{
lmom.test.ln3(data, digits=4)
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
\seealso{ \code{\link{lmomln3}}, \code{\link{parln3}} }
\examples{
lmom.test.ln3(c(123,34,4,654,37,78))

lmom.test.ln3(exp(rnorm(50)))
}
\keyword{distribution}
