\name{are.parrice.valid}
\alias{are.parrice.valid}
\title{Are the Distribution Parameters Consistent with the Rice Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfrice}}, \code{\link{pdfrice}}, \code{\link{quarice}}, and
\code{\link{lmomrice}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parrice.valid}}
function.
}
\usage{
are.parrice.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parrice}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{rice} consistent.}
  \item{FALSE}{If the parameters are not \code{rice} consistent.}
}
\note{
This function calls \code{\link{is.rice}} to verify consistency between the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{ \code{\link{is.rice}}, \code{\link{parrice}} }
\examples{
#para <- parrice(lmoms(c(123,34,4,654,37,78)))
#if(are.parrice.valid(para)) Q <- quarice(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Rice}
\keyword{utility (distribution/parameter validation)}

