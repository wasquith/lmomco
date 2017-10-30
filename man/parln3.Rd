\name{parln3}
\alias{parln3}
\title{Estimate the Parameters of the 3-Parameter Log-Normal Distribution}
\description{
This function estimates the parameters (\eqn{\zeta}, lower bounds; \eqn{\mu_{\mathrm{log}}}, location; and \eqn{\sigma_{\mathrm{log}}}, scale) of the Log-Normal3 distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomln3}}. The function uses algorithms of the Generalized Normal for core computations. Also, if \eqn{\tau_3 \le 0}, then the Log-Normal3 distribution can not be fit, however reversing the data alleviates this problem.
}
\usage{
parln3(lmom, zeta=NULL, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{zeta}{Lower bounds, if \code{NULL} then solved for.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\details{
Let the L-moments by in variable \code{lmr}, if the \eqn{\zeta} (lower bounds) is unknown, then the algorithms return the same fit as the Generalized Normal will attain. However, \code{\link{pargno}} does not have intrinsic control on the lower bounds and \code{\link{parln3}} does. The \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\tau_3} are used in the fitting for \code{\link{pargno}} and \code{\link{parln3}} but only \eqn{\lambda_1} and \eqn{\lambda_2} are used when the \eqn{\zeta} is provided as in \code{parln3(lmr, zeta=0)}. In otherwords, if \eqn{\zeta} is known, then \eqn{\tau_3} is not used and shaping comes from the choice of \eqn{\zeta}. 
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{ln3}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parln3}.}
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{lmomln3}},
\code{\link{cdfln3}}, \code{\link{pdfln3}}, \code{\link{qualn3}}, \code{\link{pargno}}
}
\examples{
lmr <- lmoms(rnorm(20))
parln3(lmr)

\dontrun{
# Handling condition of negative L-skew
# Data reversal looks like: Y <- -X, but let us use an example
# on the L-moments themselves.
lmr.pos <- vec2lmom(c(100, 45, -0.1)) # parln3(lmr.pos) fails
lmr.neg <- lmr.pos
lmr.neg$lambdas[1] <- -lmr.neg$lambdas[1]
lmr.neg$ratios[3]  <- -lmr.neg$ratios[3]
F <- nonexceeds()
plot(F, -qualn3(1-F, parln3(lmr.neg)), type="l", lwd=3, col=2) # red line
lines(F, quagno(F, pargno(lmr.pos))) # black line
}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Log-Normal (3-parameter)}
