\encoding{utf8}
\name{parpdq4}
\alias{parpdq4}
\title{Estimate the Parameters of the Polynomial Density-Quantile4 Distribution}
\description{
This function estimates the parameters of the Polynomial Density-Quantile4 distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between the distribution parameters and L-moments are seen under \code{\link{lmompdq4}}.
}
\usage{
parpdq4(lmom, checklmom=TRUE, snapt4uplimit=TRUE)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is unlikely that the L-moments will not be viable. However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{snapt4uplimit}{A logical controlling the behavior of the function for \eqn{\tau_4} exceeding an operational upper margin and whether the incoming \eqn{\tau_4} can be snapped down to this margin (see \bold{Note}).}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{pdq4}.}
  \item{para}{The parameters of the distribution.}
  \item{ifail}{A numeric field connected to the \code{ifailtext}; a value of 0 indicates fully successful operation of the function.}
  \item{ifailtext}{A message, instead of a warning, about the internal operations or operational limits of the function.}
  \item{source}{The source of the parameters: \dQuote{parpdq4}.}
}
\note{
\bold{Upper Limit of the Shape Parameter}---The following is a study of the performance of \code{parpdq4} as the upper limit of the shape parameter \eqn{\kappa} is approached. The algorithms have the ability to estimate the \eqn{\kappa} reliabily, it is the scale parameter \eqn{\alpha} that breaks down and hence there is a hard-wired setting of \eqn{\kappa > 0.99} in which a message is issued to \code{ifail} about \eqn{\alpha} reliability:
\preformatted{
  A <- 100
  K <- seq(0.8, 1, by=0.0001)
  As <- Ks <- rep(NA, length(K))
  for(i in 1:length(K)) {
    para  <- list(para=c(0, A, K[i]), type="pdq4")
    pdq4  <- parpdq4(lmompdq4(para), snapt4uplimit=FALSE)
    As[i] <- pdq4$para[2]
    Ks[i] <- pdq4$para[3]
  }
  plot( K, (As-A)/A, type="l", col="red")
  abline(v=0.99) # heuristically determined threshold
}

\bold{Lower Limit of the Shape Parameter}---The lower limit of \eqn{\kappa} does not really exist but as \eqn{\kappa \rightarrow -\infty}, the qualty of the \eqn{\alpha} operation will degrade. The approach in the code involves an \R function \code{uniroot()} operation and the lower limit is not set to \code{-Inf} but is set within sources as the value \cr \code{-.Machine$double.xmax^(1/64)},\cr which is not too small of a number, but the \eqn{\tau_4} associated with this limit is -0.2499878576145593, which is extremely close to \eqn{\tau_4 > -1/4} lower limit. The implementation here will snap incoming \eqn{\tau_4} to a threshold towards zero as
\preformatted{
  TAU4 <- "users tau4"
  smallTAU4 <- -0.2499878576145593
  if(TAU4 < smallTAU4) TAU4 <- smallTAU4 + sqrt(.Machine$double.eps)
  print(TAU4, 16) # -0.2499878427133981
}
and this snapping produces an operational lower bounds of \eqn{\kappa} of -65455.6715146775. This topic can be explored by operations such as
\preformatted{
  # Have tau4 but with internals to protect quality of the
  # alpha estimation and speed root-solving the kappa, there
  # is an operational lower bounds of tau4. Here lower limit
  # tau4 = -0.25 and the operations below return -0.2499878.
  lmompdq4(parpdq4(vec2lmom(c(0, 100, 0, -1/4))))$ratios[4]
}

\bold{Upper Operational Limit of L-kurtosis}---The script below explores the operational limit of \eqn{\tau_4} within the algorithms themselves. It is seen in the computations that breakdown in the reverse computation of the \eqn{\tau_4} from the parameters begins at \eqn{\tau_4 >= 0.867}. As a result, the argument \code{snapt4upmargin} by default and convenience could trigger snapping the solution to this upper limit (see section \bold{Even Lower Maximum Operational Limit of L-kurtosis}).
\preformatted{
  T4s <- seq(0.8, 0.9, by=0.001) # sweeping through very high Tau4
  unit_std <- 1/sqrt(pi)
  FF <- pnorm(seq(-6, 6, by=0.01))
  plot(0,0, type="n", xlim=range(qnorm(FF)), ylim=c(-6, 6),
            xlab="Standard Normal Variate", ylab="Quantile")
  for(i in 1:length(T4s)) {
    lmr  <- vec2lmom(c(0, unit_std, 0, T4s[i]))
    pdq4 <- parpdq4(lmr, snapt4uplimit=FALSE)
    lmr4 <- lmompdq4(pdq4)
    lines(qnorm(FF), quapdq4(FF, pdq4))
    err1 <- theoLmoms(pdq4)$lambdas[2] - unit_std
    err2 <-            lmr4$lambdas[2] - unit_std
    vals <- c(T4s[i], pdq4$para[3], err1, err2)
    names(vals) <- c("Tau4", "Kappa", "Err1(theoLmoms)", "Err2(lmompdq4)")
    print(vals) # both methods of Lambda2 estimation
  } # working and degenerates at Tau4 >= 0.867, so use 0.866 as a margin
}
The problem geometrically is, as the \eqn{\tau_4} becomes very \dQuote{large}, that the distribution is become so peaked that its variation will be degenerating to zero, which is not compatible with the infinite limits of the distribution. Presumably beyond \eqn{\tau_4 >= 0.867}, the TL-moments could be used with further algorithmic development. There are other difficulties though in the next example as \eqn{\tau_4} gets large.

\bold{Even Lower Maximum Operational Limit of L-kurtosis}---Further study of the limits of maximum operational limit of \eqn{\tau_4} can be made for reliable use of the basic internal functions of \R. Consider the following code:
\preformatted{
  T4s <- seq(0.4, 0.9, by=0.002)
  errs <- vector(mode="numeric", length(T4s))
  for(i in 1:length(T4s)) {
    lmra <- vec2lmom(c(0, 1, 0, T4s[i]))
    para <- parpdq4(lmra, snapt4uplimit=FALSE)
    lmrb <- lmompdq4(para)
    errs[i] <- abs(lmra$lambdas[4] - lmrb$lambdas[4])/lmra$lambdas[4]
    print(c(T4s[i], errs[i], para$para[3]))
  }
  plot(T4s, errs, ylab="abs(Lambda4 - EstLambda4)/Lambda4", col="red")
  abline(v=0.845) # so use 0.845 as a lower margin
}
The \eqn{\tau_4 >= 0.845} is therefore a more defensive upper limit for operational purposes of the \pkg{lmomco} package.

\bold{Lower Limit Performance of L-kurtosis}---The lower limit of \eqn{\tau_4 = -1/4} for the distribution is a statement of pure bimodality (two sides of a coin, as a matter of speaking). Visualization of the quantile function at the lower limit of \eqn{\tau_4} in the recipe that follows shows this fact with two flat-line segments of solid red lines with the change over at right angles at standard normal variate of zero. Then the \eqn{\tau_4} is nudged away from the lower limit just a little and replotted as the dashed line. Two other lines, but still for \eqn{\tau_4 < 0}, are shown in red and dark green. Finally, the demonstration ends with a magenta line for \eqn{\tau_4 = 0}.
\preformatted{
  FF <- pnorm(seq(-6, 6, by=0.01))
  plot(0,0, type="n", xlim=range(qnorm(FF)), ylim=c(-6, 6),
            xlab="Standard Normal Variate", ylab="Quantile")
  pdq4 <- parpdq4(vec2lmom(c(0, 1/sqrt(pi), 0, -1/4     )))
  lines(qnorm(FF), quapdq4(FF, pdq4), col="red"   )
  pdq4 <- parpdq4(vec2lmom(c(0, 1/sqrt(pi), 0, -1/4+0.03)))
  lines(qnorm(FF), quapdq4(FF, pdq4), col="red", lty=2) # dashed
  pdq4 <- parpdq4(vec2lmom(c(0, 1/sqrt(pi), 0, -1/8     )))
  lines(qnorm(FF), quapdq4(FF, pdq4), col="darkgreen")
  pdq4 <- parpdq4(vec2lmom(c(0, 1/sqrt(pi), 0, -1/16    )))
  lines(qnorm(FF), quapdq4(FF, pdq4), col="blue"   )
  pdq4 <- parpdq4(vec2lmom(c(0, 1/sqrt(pi), 0, 0        )))
  lines(qnorm(FF), quapdq4(FF, pdq4), col="magenta")
}
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmompdq4}}, \code{\link{cdfpdq4}}, \code{\link{pdfpdq4}}, \code{\link{quapdq4}}
}
\examples{
# Normal, Hosking (2007, p.2883)
para <- list(para=c(0, 0.4332, -0.7029), type="pdq4")
parpdq4(lmompdq4(para))$para
# parameter reversion shown

para <- list(para=c(0, 0.4332,  0.7029), type="pdq4")
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
  parpdq4(lmompdq4(para))$para # }
}
\keyword{distribution (parameters)}
\keyword{Distribution: Polynomial Density-Quantile4}
