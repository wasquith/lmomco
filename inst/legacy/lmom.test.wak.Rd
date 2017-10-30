\name{lmom.test.wak}
\alias{lmom.test.wak}
\title{ Test L-moment and Parameter Algorithms of the Wakeby Distribution }
\description{
This function computes the L-moments of the data and the parameters of the
Wakeby distribution and in turn computes the L-moments from the fitted 
parameters.
}
\usage{
lmom.test.wak(data, digits=4)
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
\seealso{ \code{\link{lmomwak}}, \code{\link{parwak}} }
\examples{
lmom.test.wak(c(123,34,4,654,37,78))

lmom.test.wak(rnorm(50))
}
\keyword{distribution}
