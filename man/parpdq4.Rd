\encoding{utf8}
\name{parpdq4}
\alias{parpdq4}
\title{Estimate the Parameters of the Polynomial-Density Quantile4 Distribution}
\description{
This function estimates the parameters of the Polynomial-Density Quantile4 distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between the distribution parameters and L-moments are seen under \code{\link{lmompdq4}}.
}
\usage{
parpdq4(lmom, checklmom=TRUE)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is unlikely that the L-moments will not be viable. However, for some circumstances or large simulation exercises then one might want to bypass this check.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{pdq4}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parpdq4}.}
}
\note{
The following is a study of the performance of \code{parpdq4} as the upper limit of the shape parameter \eqn{\kappa} is approached. The algorithms have the ability to estimate the \eqn{\kappa} reliabily, it is the scale parameter \eqn{\alpha} that breaks down and hence there is a hard-wired setting of \eqn{\kappa > 0.98} in which a warning is issue in \code{parpdq4} about \eqn{\alpha} reliability:
\preformatted{
  A <- 10
  K <- seq(0.8, 1, by=0.0001)
  As <- Ks <- rep(NA, length(K))
  for(i in 1:length(K)) {
    para <- list(para=c(0, A, K[i]), type="pdq4")
    As[i] <- parpdq4(lmompdq4(para))$para[2]
    Ks[i] <- parpdq4(lmompdq4(para))$para[3]
  }
  plot( K, (As-A)/A, type="l", col="red")
  abline(v=0.98) # heuristically determined threshold
}

The lower limit of \eqn{\kappa} does not really exist but as \eqn{\kappa \rightarrow -\infty}, the qualty of the \eqn{\alpha} operation will degrade. The approach in the code involves an \R function \code{uniroot()} operation and the lower limit is not set to \code{-Inf} but is set to \code{-.Machine$double.xmax^(1/64)}, which is not too small of a number, but the \eqn{\tau_4} associated with this limit is -0.2499878576145593, which is extremely close to \eqn{\tau_4 > -1/4} lower limit. The implementation here will snap incoming \eqn{\tau_4} to a threshold towards zero as
\preformatted{
  smallTAU4 <- -0.2499878576145593
  if(TAU4 < smallTAU4) TAU4 <- smallTAU4 + sqrt(.Machine$double.eps)
  print(TAU4, 16) # -0.2499878427133981
}
and this snapping produces an operational lower bounds of \eqn{\kappa} of
-65455.6715146775. This topic can be explored by operations such as
\preformatted{
  # Have tau4 but with internals to protect quality of the
  # alpha estimation and speed root-solving the kappa, there
  # is an operational lower bounds of tau4. Here lower limit
  # tau4 = -0.25 and the operations below return -0.2499878.
  lmompdq4(parpdq4(vec2lmom(c(0, 100, 0, -1/4))))$ratios[4]
}
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2,870--2,891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmompdq4}},
\code{\link{cdfpdq4}}, \code{\link{pdfpdq4}}, \code{\link{quapdq4}}
}
\examples{
# Normal, Hosking (2007, p.2883)
para <- list(para=c(0, 0.4332,-0.7029), type="pdq4")
parpdq4(lmompdq4(para))$para
# parameter reversion shown

para <- list(para=c(0, 0.4332, 0.7029), type="pdq4")
parpdq4(lmompdq4(para))$para
# parameter reversion shown with sign change kappa

\dontrun{
# other looks disabled for check --timings
para <- list(para=c(0, 0.4332, 0.97), type="pdq4")
parpdq4(lmompdq4(para))$para
# see now that alpha changing in fourth decimal as kappa
# approaches the 0.98 threshold (see Note)

# make two quick checks near zero and then zero
para <- list(para=c(0, 0.4332, +0.0001), type="pdq4")
parpdq4(lmompdq4(para))$para
para <- list(para=c(0, 0.4332, -0.0001), type="pdq4")
parpdq4(lmompdq4(para))$para
para <- list(para=c(0, 0.4332, 0), type="pdq4")
parpdq4(lmompdq4(para))$para }
}
\keyword{distribution (parameters)}
\keyword{Distribution: Polynomial-Density Quantile4}
