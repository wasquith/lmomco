\name{are.parray.valid}
\alias{are.parray.valid}
\title{Are the Distribution Parameters Consistent with the Rayleigh Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfray}}, \code{\link{pdfray}}, \code{\link{quaray}}, and
\code{\link{lmomray}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parray.valid}}
function.
}
\usage{
are.parray.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parray}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{ray} consistent.}
  \item{FALSE}{If the parameters are not \code{ray} consistent.}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments:
Research Report RC12210, IBM Research Division, Yorkton Heights, N.Y.
}
\note{
This function calls \code{\link{is.ray}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.ray}}, \code{\link{parray}} }
\examples{
para <- parray(lmoms(c(123,34,4,654,37,78)))
if(are.parray.valid(para)) Q <- quaray(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Rayleigh}
\keyword{utility (distribution/parameter validation)}

