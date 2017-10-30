\name{lmom.test.all}
\alias{lmom.test.all}
\title{Test All lmom.CCC.test Functions}
\description{
This function is a dispatcher on top of \code{lmom.CCC.test} functions, where
\code{CCC} represents distribution: \code{exp}, \code{gam}, \code{gev}, \code{glo}, \code{gno}, \code{gpa},
\code{gum}, \code{kap}, \code{nor}, \code{pe3}, \code{revgum}, \code{wak}, and \code{wei}. The motivation for the function is provide an example and builtin tool
to assess the performance of the algorithms implements the L-moments for the
supported univariate distributions functions. This function is to broadly call each supported 
distribution in the library through (1) computation of the L-moments of the data, 
(2) computation of the corresponding parameters of the distribution, and
(3) return computation of the of the theoretical L-moments of the distribution.
The differences between the sample and theoretical L-moments are produced by
\code{lmom.diff}. Further the median quantile of the distribution is computed
through the quantile function, and in turn, the median nonexceedance probability
is computed through the cumulative distribution function.
}
\usage{
lmom.test.all(data)
}
\arguments{
  \item{data}{A vector of data.}
}
\value{
  This is a high level function and is not intended to return anything other
than output to the user.
}
\author{ W.H. Asquith}
\seealso{
\code{\link{lmom.diff}}, 
\code{\link{lmom.test.exp}}, 
\code{\link{lmom.test.gam}}, 
\code{\link{lmom.test.gev}}, 
\code{\link{lmom.test.glo}}, \cr
\code{\link{lmom.test.gno}}, 
\code{\link{lmom.test.gpa}}, 
\code{\link{lmom.test.gum}}, 
\code{\link{lmom.test.kap}}, \cr
\code{\link{lmom.test.nor}}, 
\code{\link{lmom.test.pe3}},
\code{\link{lmom.test.revgum}}, \cr
\code{\link{lmom.test.wak}}, and
\code{\link{lmom.test.wei}}. 
}
\examples{
lmom.test.all(c(123,34,4,654,37,78))
}
\keyword{distribution}
