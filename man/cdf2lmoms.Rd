\name{cdf2lmoms}
\alias{cdf2lmoms}
\title{ Compute L-moments from Cumulative Distribution Function }
\description{
Compute the L-moments from a cumulative distribution function. For \eqn{r \ge 1}, the L-moments can be computed by sequential calling of \code{\link{cdf2lmom}}. Consult the documentation of that function for mathematical definitions.
}
\usage{
cdf2lmoms(para, nmom=6, fdepth=0, silent=TRUE, lambegr=1, ...)
}
\arguments{
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
  \item{nmom}{The number of moments to compute. Default is 6.}
  \item{fdepth}{The depth of the nonexceedance/exceedance probabilities to determine the lower and upper integration limits through a call to the \code{\link{par2qua}} function. The default of 0 implies the quantile for \eqn{F=0} and quantile for \eqn{F=1} as the respective lower and upper limits.}
  \item{silent}{A logical to be passed into \code{\link{cdf2lmom}} and then onto the  \code{try} functions encompassing the \code{integrate} function calls.}
  \item{lambegr}{The \eqn{r}th order to begin the sequence for L-moment computation. Can be used as a means to bypass a mean computation if the user has an alternative method for the mean or other central tendency characterization in which case \code{lambegr = 2}.}
  \item{...}{Additional arguments to pass to \code{\link{cdf2lmom}}.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\hat{\lambda}^{(0,0)}_1}, second element is \eqn{\hat{\lambda}^{(0,0)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\hat{\tau}^{(0,0)}}, third element is \eqn{\hat{\tau}^{(0,0)}_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal \code{NULL} is not support for trimming is provided by this function.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which will equal \code{NULL} is not support for trimming is provided by this function.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which will equal \code{NULL} is not support for trimming is provided by this function.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{cdf2lmoms}.}
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdf2lmom}}, \code{\link{lmoms}}}
\examples{
cdf2lmoms(vec2par(c(10,40), type="ray"))
\dontrun{
# relatively slow computation
vec2par(c(.9,.4), type="emu"); cdf2lmoms(para, nmom=4)
vec2par(c(.9,.4), type="emu"); cdf2lmoms(para, nmom=4, fdepth=0)}
}
\keyword{L-moment (theoretical)}
