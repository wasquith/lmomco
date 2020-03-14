\name{are.partexp.valid}
\alias{are.partexp.valid}
\title{Are the Distribution Parameters Consistent with the Truncated Exponential Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdftexp}}, \code{\link{pdftexp}}, \code{\link{quatexp}}, and \code{\link{lmomtexp}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively. These functions internally use the \code{\link{are.partexp.valid}} function.
}
\usage{
are.partexp.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parexp}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{texp} consistent.}
  \item{FALSE}{If the parameters are not \code{texp} consistent.}
}
\note{
This function calls \code{\link{is.texp}} to verify consistency between the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\references{
Vogel, R.M., Hosking, J.R.M., Elphick, C.S., Roberts, D.L., and Reed, J.M., 2008, Goodness of fit of probability distributions for sightings as species approach extinction: Bulletin of Mathematical Biology, DOI 10.1007/s11538-008-9377-3, 19 p.
}
\seealso{ \code{\link{is.texp}}, \code{\link{partexp}} }
\examples{
para <- partexp(lmoms(c(90,134,100,114,177,378)))
if(are.partexp.valid(para)) Q <- quatexp(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Exponential (trimmed)}
\keyword{utility (distribution/parameter validation)}

