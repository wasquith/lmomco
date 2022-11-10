\name{are.pargep.valid}
\alias{are.pargep.valid}
\title{Are the Distribution Parameters Consistent with the Generalized Exponential Poisson Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfgep}}, \code{\link{pdfgep}}, \code{\link{quagep}}, and
\code{\link{lmomgep}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.pargep.valid}}
function. The parameters must be \eqn{\beta > 0}, \eqn{\kappa > 0}, and \eqn{h > 0}.
}
\usage{
are.pargep.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{pargep}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{gep} consistent.}
  \item{FALSE}{If the parameters are not \code{gep} consistent.}
}
\references{
Barreto-Souza, W., and Cribari-Neto, F., 2009, A generalization of the exponential-Poisson distribution: Statistics and Probability, 79, pp. 2493--2500.
}
\note{
This function calls \code{\link{is.gep}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\seealso{\code{\link{is.gep}}, \code{\link{pargep}} }
\examples{
#para <- pargep(lmoms(c(123,34,4,654,37,78)))
#if(are.pargep.valid(para)) Q <- quagep(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Generalized Exponential Poisson}
\keyword{utility (distribution/parameter validation)}

