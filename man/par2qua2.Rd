\name{par2qua2}
\alias{par2qua2}
\title{ Equivalent Quantile Function of Two Distributions}
\description{
This function computes the nonexceedance probability of a given quantile from a linear weighted combination of two quantile functions---a mixed distribution:
\deqn{Q_\mathrm{mixed}(F; \Theta_1, \Theta_2, \omega) = (1-\omega)Q_1(F, \Theta_1) + \omega Q_2(F, \Theta2)\mbox{,}}
where \eqn{Q} is a quantile function for nonexceedance probability \eqn{F}, the distributions have parameters \eqn{\Theta_1} and \eqn{\Theta_2}, and \eqn{\omega} is a weight factor.


The distributions are specified by the two parameter object arguments in usual \pkg{lmomco} style. When proration by the nonexceedance probability is desired (\code{weight=}\code{NULL}, default), the left-tail parameter object (\code{para1}) is the distribution obviously governing the left tail; the right-tail parameter object (\code{para2}) is of course governs the right tail. The quantile function algebra is
\deqn{ Q(F) = (1-F^\star) \times {\triangleleft}Q(F) + F^\star \times Q(F){\triangleright}\mbox{,}}
where \eqn{Q(F)} is the mixed quantile for nonexceedance probability \eqn{F}. \eqn{{\triangleleft}Q(F)} is the first or left-tail quantile function; \eqn{Q(F){\triangleright}} is the second or right-tail quantile function. In otherwords, if \code{weight = NULL}, then \eqn{F^\star = F =} \code{f} and the weight between the two quantile functions thus continuously varies from left to right. This is a probability proration from one to the other. A word of caution in this regard. The resulting weighted- or mixed-quantile function is not rigorously checked for monotonic increase with \eqn{F}, which is a required property of quantile functions. However, a first-order difference on the mixed quantiles with the probabilities is computed and a warning issued if not monotonic increasing.

If the optional \code{weight} argument is provided with length 1, then \eqn{\omega} equals that weight. If \code{weight = 0}, then only the quantiles for \eqn{Q_1(F)} are returned, and if \code{weight = 1}, then only the quantiles for the left tail \eqn{Q_2(F)} are returned.

If the optional \code{weight} argument is provided with length 2, then \eqn{(1 - \omega)} is replaced by the first weight and \eqn{\omega} is replaced by the second weight. These are internally rescaled to sum to unity before use and a warning is issued that this was done. Finally, the \code{\link{par2cdf2}} function inverses the above equation for \eqn{F}.
}
\usage{
par2qua2(f, para1, para2, wfunc=NULL, weight=NULL, as.list=FALSE, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para1}{The first or left-tail parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{para2}{The second or right-tail parameters from \code{\link{lmom2par}} or similar.}
  \item{wfunc}{A function taking the argument \code{f} and computing a weight for the \code{para2} curve for which the complement of the computed weight is used for the weight on \code{para1}.}
  \item{weight}{An optional weighting argument to use in lieu of \code{F}. If \code{NULL} then prorated by the \code{f}, if \code{weight} has length 1, then weight on left distribution is the complement of the weight and weight on right distribution is \code{weight[1]}, and if \code{weight} had length 2, then \code{weight[1]} is the weight on the left distribution and \code{weight[2]} is the weight on the right distribution.}
  \item{as.list}{A logical to control whether an \R \code{data.frame} is returned having a column for \code{f} and for the mixed quantiles. This feature is provided for some design consistency with \code{\link{par2qua2lo}}, which mandates a \code{data.frame} return.}
  \item{...}{The additional arguments are passed to the quantile function.}
}
\value{
  The weighted quantile value for \eqn{F} from the two distributions.
}
\author{ W.H. Asquith}
\seealso{\code{\link{par2qua}}, \code{\link{par2cdf2}}, \code{\link{par2qua2lo}}
}
\examples{
lmr <- lmoms(rnorm(20)); left <- parnor(lmr); right <- pargev(lmr)
mixed.median <- par2qua2(0.5, left, right)

# Bigger example--using Kappa fit to whole sample for the right tail and
# Normal fit to whole sample for the left tail
D   <- c(123, 523, 345, 356, 2134, 345, 2365, 235, 12, 235, 61, 432, 843)
lmr <- lmoms(D); KAP <- parkap(lmr); NOR <- parnor(lmr); PP <- pp(D)
plot( PP, sort(D), ylim=c(-500, 2300))
lines(PP, par2qua( PP, KAP),      col=2)
lines(PP, par2qua( PP, NOR),      col=3)
lines(PP, par2qua2(PP, NOR, KAP), col=4)
}
\keyword{quantile mixture function}
\keyword{distribution (mixture)}
\keyword{mixed distribution}
