\name{are.parln3.valid}
\alias{are.parln3.valid}
\title{Are the Distribution Parameters Consistent with the 3-Parameter Log-Normal Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfln3}}, \code{\link{pdfln3}}, \code{\link{qualn3}}, and \code{\link{lmomln3}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively.  These functions internally use the \code{\link{are.parln3.valid}} function.
}
\usage{
are.parln3.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parln3}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{ln3} consistent.}
  \item{FALSE}{If the parameters are not \code{ln3} consistent.}
}
\note{
This function calls \code{\link{is.ln3}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{is.ln3}}, \code{\link{parln3}} }
\examples{
para <- parln3(lmoms(c(123,34,4,654,37,78)))
if(are.parln3.valid(para)) Q <- qualn3(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Log-Normal (3-parameter)}
\keyword{utility (distribution/parameter validation)}

