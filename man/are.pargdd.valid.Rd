\name{are.pargdd.valid}
\alias{are.pargdd.valid}
\title{Are the Distribution Parameters Consistent with the Gamma Difference Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfgdd}}, \code{\link{pdfgdd}}, \code{\link{quagdd}}, and \code{\link{lmomgdd}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively. These functions internally use the \code{\link{are.pargdd.valid}} function.
}
\usage{
are.pargdd.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{pargdd}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{gdd} consistent.}
  \item{FALSE}{If the parameters are not \code{gdd} consistent.}
}
\references{
Klar, B., 2015, A note on gamma difference distributions: Journal of Statistical Computation and Simulation v. 85, no. 18, pp. 1--8, \doi{10.1080/00949655.2014.996566}.
}
\note{
This function calls \code{\link{is.gdd}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\seealso{\code{\link{is.gdd}}, \code{\link{pargdd}} }
\examples{
#
}
\keyword{utility (distribution)}
\keyword{Distribution: Gamma Difference}
\keyword{utility (distribution/parameter validation)}

