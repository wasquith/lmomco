\name{are.parst3.valid}
\alias{are.parst3.valid}
\title{Are the Distribution Parameters Consistent with the 3-Parameter Student t Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfst3}}, \code{\link{pdfst3}}, \code{\link{quast3}}, and \code{\link{lmomst3}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively. These functions internally use the \code{\link{are.parst3.valid}} function.
}
\usage{
are.parst3.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parst3}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{st3} consistent.}
  \item{FALSE}{If the parameters are not \code{st3} consistent.}
}
\note{
This function calls \code{\link{is.st3}} to verify consistency between the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{ \code{\link{is.st3}}, \code{\link{parst3}} }
\examples{
para <- parst3(lmoms(c(90,134,100,114,177,378)))
if(are.parst3.valid(para)) Q <- quast3(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Student t (3-parameter)}
\keyword{utility (distribution/parameter validation)}

