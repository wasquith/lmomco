\name{pargdd}
\alias{pargdd}
\title{Estimate the Parameters of the Gamma Difference Distribution}
\description{
This function estimates the parameters of the Gamma Difference distribution given the L-moments of the data in an ordinary L-moment object (\code{\link{lmoms}}). The relations between distribution parameters and L-moments are complex (see \code{\link{lmomgdd}}). The distribution has four parameters. The vector \code{para} in the parameter object with a fifth parameter uses that as a trigger between a symmetrical distribution with \code{para[3:4]} equals \code{para[1:2]} if \code{para[5] = 1}. If \code{para[5]} is not present, then the distribution can be asymmetrical, or if \code{para[5]} is present and set to any value that is not 1, then the distribution can be asymmetrical.
}
\usage{
pargdd(lmom, checklmom=TRUE, symgdd=FALSE, init.para=NULL,
             silent=FALSE,  trace=FALSE, control=list(abstol=0.001, maxit=1000), ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function.}
  \item{symgdd}{A logical to trigger a symmetrical distribution by \eqn{\alpha_2 = \alpha_1} and \eqn{\beta_1 = \beta_1} and the fifth element of \code{para} on the return will be set to 1.}
  \item{init.para}{Optional initial values for the parameters used for starting values for the \code{optim} function. If this argument is not set, then an unrigorous attempt is made to guess at the initial parameters using some poor admittedly heuristics. The fifth element, if present, and set to 1, then the \code{symdd} is internally set to true.}
  \item{silent}{The argument \code{silent} for \code{try()}.}
  \item{control}{The argument \code{control} for \code{optim()}.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gdd}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{pargdd}.}
  \item{optim}{The results of the parameter optimization call.}
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomgdd}}, \code{\link{cdfgdd}}, \code{\link{pdfgdd}}, \code{\link{quagdd}}
}
\examples{
\dontrun{
# Example of the symmetrical case, see lmomgdd-Note section.
x <- seq(-20, 20, by=0.1); para <- list(para=c(3, 0.4, NA, NA, 1), type="gdd")
slmr  <- lmomgdd(  para);  nara <- pargdd(slmr, symgdd=TRUE)
given <- pdfgdd(x, para);  fit  <- pdfgdd(x, nara)
plot( x, given, type="l", col=8, lwd=4, ylim=range(c(given, fit)))
lines(x, fit,   col="red") # }

\dontrun{
# Example of the asymmetrical case, and as of Summer 2024 experiments, it seems
# the author does not quite have limits of GDD implementation known. Though this
# example works, we do not always L-moment recreation from fitted parameters.
x <- seq(-10, 30, by=0.1); para <- list(para=c(3, 0.95, 1, 1), type="gdd")
slmr  <- lmomgdd(  para);  nara <- pargdd(slmr)
given <- pdfgdd(x, para);  fit  <- pdfgdd(x, nara)
plot( x, given, type="l", col=8, lwd=4, ylim=range(c(a, fit)))
lines(x, fit,   col="red") # }
}
\keyword{distribution (parameters)}
\keyword{Distribution: Gamma Difference}
