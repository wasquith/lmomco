\name{are.parkur.valid}
\alias{are.parkur.valid}
\title{Are the Distribution Parameters Consistent with the Kumaraswamy Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfkur}}, \code{\link{pdfkur}}, \code{quakur}, and
\code{\link{lmomkur}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parkur.valid}}
function.
}
\usage{
are.parkur.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parkur}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{kur} consistent.}
  \item{FALSE}{If the parameters are not \code{kur} consistent.}
}
\references{
Jones, M.C., 2009, Kumaraswamy's distribution---A beta-type distribution with
some tractability advantages: Statistical Methodology, v. 6, pp. 70--81.
}
\note{
This function calls \code{\link{is.kur}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.kur}}, \code{\link{parkur}} }
\examples{
para <- parkur(lmoms(c(0.25, 0.4, 0.6, 0.65, 0.67, 0.9)))
if(are.parkur.valid(para)) Q <- quakur(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Kumaraswamy}
\keyword{utility (distribution/parameter validation)}

