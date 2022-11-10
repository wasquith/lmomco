\name{are.parwei.valid}
\alias{are.parwei.valid}
\title{Are the Distribution Parameters Consistent with the Weibull Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions (\code{\link{cdfwei}}, \code{\link{pdfwei}}, \code{\link{quawei}}, and \code{\link{lmomwei}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively.  These functions internally use the \code{\link{are.parwei.valid}} function.
}
\usage{
are.parwei.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parwei}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{wei} consistent.}
  \item{FALSE}{If the parameters are not \code{wei} consistent.}
}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\note{
This function calls \code{\link{is.wei}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\seealso{\code{\link{is.wei}}, \code{\link{parwei}} }
\examples{
para <- parwei(lmoms(c(123,34,4,654,37,78)))
if(are.parwei.valid(para)) Q <- quawei(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Weibull}
\keyword{utility (distribution/parameter validation)}

