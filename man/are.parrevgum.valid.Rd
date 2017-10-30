\name{are.parrevgum.valid}
\alias{are.parrevgum.valid}
\title{Are the Distribution Parameters Consistent with the Reverse Gumbel Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfrevgum}}, \code{\link{pdfrevgum}}, \code{\link{quarevgum}}, and
\code{\link{lmomrevgum}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parrevgum.valid}}
function.
}
\usage{
are.parrevgum.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parrevgum}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{revgum} consistent.}
  \item{FALSE}{If the parameters are not \code{revgum} consistent.}
}
\references{
Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data,
in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan,
chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\note{
This function calls \code{\link{is.revgum}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.revgum}}, \code{\link{parrevgum}} }
\examples{
para <- vec2par(c(.9252, .1636, .7),type='revgum')
if(are.parrevgum.valid(para)) Q <- quarevgum(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Gumbel (reversed)}
\keyword{utility (distribution, parameter validation)}

