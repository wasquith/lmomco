\name{parrice}
\alias{parrice}
\title{Estimate the Parameters of the Rice Distribution}
\description{
This function estimates the parameters (\eqn{\nu} and \eqn{\alpha}) of the Rice distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are complex and tabular lookup is made using a relation between \eqn{\tau} and a form of signal-to-noise ratio \eqn{\mathrm{SNR}} defined as \eqn{\nu/\alpha} and a relation between \eqn{\tau} and precomputed Laguerre polynomial (\code{\link{LaguerreHalf}}).

The \eqn{\lambda_1} (mean) is most straightforward
\deqn{
\lambda_1 = \alpha \times \sqrt{\pi/2} \times L_{1/2}(-\nu^2/[2\alpha^2])\mbox{,}
}
for which the terms to the right of the multiplication symbol are uniquely a function of \eqn{\tau} and precomputed for tabular lookup and interpolation from \file{sysdata.rdb} (\env{.lmomcohash$RiceTable}). Parameter estimation also relies directly on tabular lookup and interpolation to convert \eqn{\tau} to \eqn{\mathrm{SNR}}. The file \file{SysDataBuilder.R} provides additional technical details.
}
\usage{
parrice(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check. However, the end point of the Rice distribution for high \eqn{\nu/\alpha} is not determined here, so it is recommended to leave \code{checklmom} turned on.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{rice}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parrice}.}
  \item{ifail}{A numeric failure mode.}
  \item{ifailtext}{A helpful message on the failure.}
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{lmomrice}}, \code{\link{cdfrice}}, \code{\link{pdfrice}}, \code{\link{quarice}}
}
\examples{
\dontrun{ 
  parrice(lmomrice(vec2par(c(10,50),   type="rice"))) # Within Rician limits
  parrice(lmomrice(vec2par(c(100,0.1), type="rice"))) # Beyond Rician limits

plotlmrdia(lmrdia(), xlim=c(0,0.2), ylim=c(-0.1,0.22),
           autolegend=TRUE, xleg=0.05, yleg=0.05)
lines(.lmomcohash$RiceTable$TAU3, .lmomcohash$RiceTable$TAU4, lwd=5, col=8)
legend(0.1,0, "RICE DISTRIBUTION", lwd=5, col=8, bty="n")
text(0.14,-0.04, "Normal distribution limit on left end point")
text(0.14,-0.055, "Rayleigh distribution limit on right end point")

# check parrice against a Maximum Likelihood method in VGAM
set.seed(1)
library(VGAM) # now example from riceff() of VGAM
vee <- exp(2); sigma <- exp(1); y <- rrice(n <- 1000, vee, sigma)
fit <- vglm(y ~ 1, riceff, trace=TRUE, crit="c")
Coef(fit)
# NOW THE MOMENT OF TRUTH, USING L-MOMENTS
parrice(lmoms(y))
# VGAM package 0.8-1 reports
#     vee    sigma
# 7.344560 2.805877
# lmomco 1.2.2 reports
#      nu    alpha
# 7.348784 2.797651}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Rice}
