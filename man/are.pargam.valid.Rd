\name{are.pargam.valid}
\alias{are.pargam.valid}
\title{Are the Distribution Parameters Consistent with the Gamma Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfgam}}, \code{\link{pdfgam}}, \code{\link{quagam}}, and
\code{\link{lmomgam}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.pargam.valid}}
function.  The parameters are restricted to the following conditions.
\deqn{\alpha > 0 \mbox{ and } \beta > 0\mbox{.}}

Alternatively, a three-parameter version is available following the parameterization of the Generalized Gamma distribution used in the \pkg{gamlss.dist} package and  and for \pkg{lmomco} is documented under \code{\link{pdfgam}}. The parameters for this version are 
\deqn{\mu > 0;\;\; \sigma > 0;\;\; -\infty < \nu < \infty}
for parameters number 1, 2, 3, respectively.
}
\usage{
are.pargam.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{pargam}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{gam} consistent.}
  \item{FALSE}{If the parameters are not \code{gam} consistent.}
}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\note{
This function calls \code{\link{is.gam}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.gam}}, \code{\link{pargam}} }
\examples{
para <- pargam(lmoms(c(123,34,4,654,37,78)))
if(are.pargam.valid(para)) Q <- quagam(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Gamma}
\keyword{utility (distribution/parameter validation)}
\keyword{Distribution: Generalized Gamma}

