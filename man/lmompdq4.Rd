\encoding{utf8}
\name{lmompdq4}
\alias{lmompdq4}
\title{L-moments of the Polynomial Density-Quantile4 Distribution}
\description{
This function estimates the L-moments of the Polynomial Density-Quantile4 distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) from \code{\link{parpdq4}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi\mbox{,}}
\deqn{\lambda_2 = \frac{\alpha}{\kappa} \bigl(1-\kappa^2\bigr)\, \mathrm{atanh}(\kappa)\mathrm{\ for\ } \kappa > 0\mbox{,}}
\deqn{\lambda_2 = \frac{\alpha}{\kappa} \bigl(1+\kappa^2\bigr)\, \mathrm{atan}(\kappa)\mathrm{\ for\ } \kappa < 0\mbox{,}}
\deqn{\tau_3 = 0 \mbox{, and}}
\deqn{\tau_4 = -\frac{1}{4} + \frac{5}{4\kappa}\biggl(\frac{1}{\kappa} - \frac{1}{\mathrm{atanh}(\kappa)} \biggr) \mathrm{\ for\ } \kappa > 0\mbox{,}}
\deqn{\tau_4 = -\frac{1}{4} - \frac{5}{4\kappa}\biggl(\frac{1}{\kappa} - \frac{1}{\mathrm{atan}(\kappa)} \biggr) \mathrm{\ for\ } \kappa < 0\mbox{,}}
}
\usage{
lmompdq4(para, paracheck=TRUE)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{ifail}{A numeric field connected to the \code{ifailtext}; a value of 0 indicates fully successful operation of the function.}
  \item{ifailtext}{A message, instead of a warning, about the internal operations or operational limits of the function.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmompdq4}.}
}
\note{
\bold{What L-kurtosis produces the widest 95th-percentile bounds?}---Study of the shapes of the PDQ4 will show that with support for \eqn{\tau_4} much less and even negative and much more than the \eqn{\tau_4 = 0.122602} defined into the Normal distribution considerable variation. The widths or spreads between quantiles moderately deep into the tails might be interesting to study. Consider the code that follows that seeks the \eqn{\tau_4} that will produce the widest 95th-percentile bounds:
\preformatted{
  ofunc <- function(t4,  lscale=NA) {
    lmr <- vec2lmom(c(0, lscale, 0, t4))
    if(! are.lmom.valid(lmr)) return(-Inf)
    pdq4  <- lmomco::parpdq4(lmr, snapt4uplimit=FALSE)
    return(-diff(lmomco::quapdq4(c(0.025, 0.975), pdq4)))
  }
  optim(0.2, ofunc, lscale=1)$par # [1] 0.4079688
}

The code maximizes at about \eqn{\tau_4 = 0.4079688}. It is informative to visualizing the nature of the objective function. In the code below, we standardize the width by division of the \eqn{\lambda_2 = 1} for generality and because of symmetry only the 97.5th percentile requires study:
\preformatted{
  lscale <- 1
  tau4s  <- seq(-1/4, 0.9, by=0.01)
  qua975s <- rep(NA, length(tau4s))
  for(i in 1:length(tau4s)) {
    lmr <- vec2lmom(c(0, lscale, 0, tau4s[i]))
    if(! are.lmom.valid(lmr)) next
    pdq4 <- lmomco::parpdq4(lmr, snapt4uplimit=FALSE)
    quas <- lmomco::quapdq4(c(0.025, 0.975), pdq4)
    qua975s[i] <- quas[2] / lscale
  }
  plot(tau4s, qua975s, ylim=c(-0.1, 5), col="blue")
  abline(v=0.845, lty=2) # supporting the "snaptau4uplimit" in parpdq4().
  abline(v=0.4079688, col=2, lwd=2)
  abline(h=qnorm(0.975, sd=sqrt(pi)), col="green", lty=3, lwd=3)
}

The figure so produces shows that the maximum at the red vertical line for \eqn{\tau_4} is at the crest of the blue points. The figure shows that for \eqn{\tau_4 >= 0.845} that numerical problems manifest and contribute to an snapping limit of \eqn{\tau_4} in \code{\link{parpdq4}}. The figure also shows with a dotted green line that the equivalent percentile of the Normal distribution with a standard deviation equivalent to the \eqn{\lambda_2 = 1} has two intersections on the widths of the PDQ4.

Now some further experiments on the apparent computational limits to \eqn{\tau_4} can be made using the code that follows. This support the threshold of \eqn{\tau_4 \le 0.845} embedded into \code{\link{parpdq4}} through the use of the \code{\link{theoTLmoms}} function.
\preformatted{
  t4s <- seq(-1/4, 1, by=0.02)
  t4s <- t4s[t4s > -1/4 & t4s < 1]
  l2s_theo <- t4s_theo <- t6s_theo <- rep(NA, length(t4s))
  for(i in 1:length(t4s)) {
    lmr  <- vec2lmom(c(0, 1, 0, t4s[i]))
    suppressWarnings(par <- parpdq4(lmr, snapt4uplimit=FALSE))
    tlmr <- theoTLmoms(par, nmom=6, trim=0)
    l2s_theo[i] <- tlmr$lambdas[2]
    t4s_theo[i] <- tlmr$ratios[ 4]
    t6s_theo[i] <- tlmr$ratios[ 6]
  }
  plot(  t4s_theo, l2s_theo, type="l")
  points(t4s_theo, l2s_theo)
    abline(v=0.864, lty=2) # see "snaptau4uplimit" in parpdq4()
    abline(v=0.845, lty=2) # see "snaptau4uplimit" in parpdq4()
  plot(  t4s_theo, t4s,      type="l")
  points(t4s_theo, t4s)
    abline(v=0.864, lty=2) # see "snaptau4uplimit" in parpdq4()
    abline(v=0.845, lty=2) # see "snaptau4uplimit" in parpdq4()
  plot(  t4s_theo, t6s_theo, type="l")
  points(t4s_theo, t6s_theo)
    abline(v=0.864, lty=2) # see "snaptau4uplimit" in parpdq4()
    abline(v=0.845, lty=2) # see "snaptau4uplimit" in parpdq4()
}
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{parpdq4}}, \code{\link{cdfpdq4}}, \code{\link{pdfpdq4}}, \code{\link{quapdq4}} }
\examples{
para <- vec2par(c(0, 1, -100), type="pdq4")
lmompdq4(  para)$ratios[4]                 # -0.2421163
theoTLmoms(para, nmom=6, trim=0)$ratios[4] # -0.2421163
theoTLmoms(para, nmom=6, trim=1)$ratios[4] # -0.2022106
theoTLmoms(para, nmom=6, trim=2)$ratios[4] # -0.1697186

\dontrun{
  para <- list(para=c(20, 1, -0.5), type="pdq4")
  lmoms(quapdq4(runif(100000), para))$lambdas
  lmompdq4(para)$lambdas #}

\dontrun{
  para <- list(para=c(20, 1, +0.5), type="pdq4")
  lmoms(quapdq4(runif(100000), para))$lambdas
  lmompdq4(para)$lambdas #}

\dontrun{
  K1 <- seq(-5, 0, by=0.001)
  K2 <- seq( 0, 1, by=0.001)
  suppressWarnings(mono_decrease_part1 <- -(1/4) + (5/(4*K1)) * (1/K1 - 1/atanh(K1)))
                   mono_increase_part2 <- -(1/4) - (5/(4*K1)) * (1/K1 - 1/atan( K1))
                   mono_increase_part1 <- -(1/4) + (5/(4*K2)) * (1/K2 - 1/atanh(K2))
                   mono_decrease_part2 <- -(1/4) - (5/(4*K2)) * (1/K2 - 1/atan( K2))

  plot( 0, 0, type="n", xlim=range(c(K1, K2)), ylim=c(-0.25, 1),
       xlab="Kappa shape parameter PDQ4 distribution", ylab="L-kurtosis (Tau4)")
  lines(K1, mono_decrease_part1, col=4, lwd=0.3)
  lines(K2, mono_increase_part1, col=4, lwd=3)
  lines(K2, mono_decrease_part2, col=2, lwd=0.3)
  lines(K1, mono_increase_part2, col=2, lwd=3)

  abline(h= 1/6, lty=2, lwd=0.6)
  abline(h=-1/4, lty=2, lwd=0.6)
  text(-5, -1/4, "Tau4 lower bounds", pos=4, cex=0.8)
  abline(v=0,    lty=2, lwd=0.6)
  abline(v=1,    lty=1, lwd=0.9)
  points(-0.7029, 0.1226, pch=15, col="darkgreen")

  # bigTAU4 <- 0.845 # see parpdq4.R and parpdq4.Rd
  pdq4 <- parpdq4(vec2lmom(c(0, 1, 0, 0.845)), snapt4uplimit=FALSE)
  points(pdq4$para[3], 0.845, cex=1.5, pch=17, col="blue")

  legend("topleft", c("Monotonic increasing for kappa < 0 (used for PDQ4)",
                      "Monotonic increasing for kappa > 0 (used for PDQ4)",
                      "Monotonic decreasing for kappa > 0 (not used for PDQ4)",
                      "Monotonic decreasing for kappa < 0 (not used for PDQ4)",
                      "Normal distribution (Tau4=0.122602 by definition)",
                      "Operational upper limit of Tau4 before numerical problems"), cex=0.8,
     pch=c(NA, NA, NA, NA, 15, 17), lwd=c(3,3, 0.3, 0.3, NA, NA),
     pt.cex=c(NA, NA, NA, NA, 1, 1.5), col=c(2, 4, 2, 4, "darkgreen", "blue")) # }
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Polynomial Density-Quantile4}
