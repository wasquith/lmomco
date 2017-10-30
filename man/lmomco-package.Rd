\name{lmomco-package}
\docType{package}
\alias{lmomco-package}
\title{L-moments, Censored L-moments, Trimmed L-moments, L-comoments, and Many Distributions}
\description{
   The \pkg{lmomco} package is a comparatively comprehensive implementation of L-moments in addition to probability-weighted moments, and parameter estimation for numerous familiar and not-so-familiar distributions. L-moments and their cousins are based on certain linear combinations of order statistic expectations. Being based on linear mathematics and thus especially robust compared to conventional moments, they are particular suitable for analysis of rare events of non-Normal data. L-moments are consistent and often have smaller sampling variances than maximum likelihood in small to moderate sample sizes. L-moments are especially useful in the context of quantile functions. The method of L-moments (\code{\link{lmr2par}}) is augmented here with access to the methods of maximum likelihood (\code{\link{mle2par}}) and maximum product of spacings (\code{\link{mps2par}}) as alternatives bound against the distributions of the \pkg{lmomco} package.

About 350 user-level functions are implemented in \pkg{lmomco} that range from low-level utilities forming an application programming interface (API) to high-level sophisticated data analysis and visualization operators. The \dQuote{See Also} section lists recommended function entry points for new users. The nomenclature (\code{d}, \code{p}, \code{r}, \code{q})-\code{lmomco} is directly analogous to that for distributions built-in to \R. To conclude, the \R packages \pkg{lmom} (Hosking), \pkg{lmomRFA} (Hosking), \pkg{Lmoments} (Karvanen) might also be of great interest.

How does \pkg{lmomco} basically work? The design of \pkg{lmomco} is to fit distributions to the L-moments of sample data. Distributions are specified by a \code{type} argument for very many functions. The package stores both L-moments (see \code{\link{vec2lmom}}) and parameters (see \code{\link{vec2par}}) in simple \R \code{list} structures---very elementary. The following code shows a comparison of parameter estimation for a random sample (\code{\link{rlmomco}}) of a GEV distribution using L-moments (\code{\link{lmoms}} coupled with \code{\link{lmom2par}} or simply \code{\link{lmr2par}}), maximum likelihood (MLE, \code{\link{mle2par}}), and maximum product of spacings (MPS, \code{\link{mps2par}}). (A note of warning, the MLE and MPS algorithms might not converge with the initial parameters---for purposes of \dQuote{learning} about this package just rerun the code below again for another random sample.)
\preformatted{
  parent.lmoments <- vec2lmom(c(3.08, 0.568, -0.163)); ty <- "gev"
  Q <- rlmomco(63, lmom2par(parent.lmoments, type=ty)) # random sample
  init <- lmoms(Q); init$ratios[3] <- 0 # failure rates for mps and mle are
  # substantially lowered if starting from the middle of the distribution's
  # shape to form the initial parameters for para.int
  lmr  <- lmr2par(Q, type=ty)                 # method of L-moments
  mle  <- mle2par(Q, type=ty, para.init=init) # method of MLE
  mps  <- mps2par(Q, type=ty, para.init=init) # method of MPS
  lmr1 <- lmr$para; mle1 <- mle$para; mps1 <- mps$para
}
The \code{lmr1}, \code{mle1}, and \code{mps1} variables each contain distribution parameter estimates, but before they are inspected, how about quick comparison to another \R package (\pkg{eva})?
\preformatted{
  lmr2 <- eva::gevrFit(Q, method="pwm")$par.ests # PWMs == L-moments
  mle2 <- eva::gevrFit(Q, method="mle")$par.ests # method of MLE
  mps2 <- eva::gevrFit(Q, method="mps")$par.ests # method of MPS
  # Package eva uses a different sign convention on the GEV shape parameter
  mle2[3] <- -mle2[3]; mps2[3] <- -mps2[3]; lmr2[3] <- -lmr2[3];
}
Now let us inspect the contents of the six estimates of the three GEV parameters by three different methods:
\preformatted{
  message("LMR(lmomco): ", paste(round(lmr1, digits=5), collapse="  "))
  message("LMR(   eva): ", paste(round(lmr2, digits=5), collapse="  "))
  message("MLE(lmomco): ", paste(round(mle1, digits=5), collapse="  "))
  message("MLE(   eva): ", paste(round(mle2, digits=5), collapse="  "))
  message("MPS(lmomco): ", paste(round(mps1, digits=5), collapse="  "))
  message("MPS(   eva): ", paste(round(mps2, digits=5), collapse="  "))
}
The results show compatible estimates between the two packages. Lastly, let us plot what these distributions look like using the \pkg{lmomco} functions: \code{\link{add.lmomco.axis}}, \code{\link{nonexceeds}}, \code{\link{pp}}, and \code{\link{qlmomco}}.
\preformatted{
  par(las=2, mgp=c(3,0.5,0)); FF <- nonexceeds(); qFF <- qnorm(FF)
  PP <- pp(Q); qPP <- qnorm(PP); Q <- sort(Q)
  plot(  qFF, qlmomco(FF, lmr), xaxt="n", xlab="", tcl=0.5,
                                ylab="QUANTILE", type="l")
  lines( qFF, qlmomco(FF, mle), col="blue")
  lines( qFF, qlmomco(FF, mps), col="red")
  points(qPP, Q, lwd=0.6, cex=0.8, col=grey(0.3)); par(las=1)
  add.lmomco.axis(las=2, tcl=0.5, side.type="NPP")
}
}
\references{
Asquith, W.H., 2007, L-moments and TL-moments of the generalized lambda distribution: Computational Statistics and Data Analysis, v. 51, no. 9, pp. 4484--4496, \url{http://dx.doi.org/10.1016/j.csda.2006.07.016}

Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8, \url{http://www.amazon.com/dp/1463508417/}

Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970, \url{http://dx.doi.org/10.1016/j.csda.2012.12.013}

Dey, D.K., Roy, Dooti, Yan, Jun, 2016, Univariate extreme value analysis, chapter 1, \emph{in} Dey, D.K., and Yan, Jun, eds., Extreme value modeling and risk analysis---Methods and applications: Boca Raton, FL, CRC Press, pp. 1--22.

Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational statistics and data analysis, vol. 43, pp. 299-314,  \url{http://dx.doi.org/10.1016/S0167-9473(02)00250-5}

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124, \url{http://www.jstor.org/stable/2345653}

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press, \url{http://www.amazon.com/dp/0521019400/}

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York, \url{http://www.amazon.com/dp/0817683607/}

Serfling, R., and Xiao, P., 2007, A contribution to multivariate L-moments---L-comoment matrices: Journal of Multivariate Analysis, v. 98, pp. 1765--1781, \url{http://dx.doi.org/10.1016/j.jmva.2007.01.008}
}
\author{William Asquith \email{william.asquith@ttu.edu} }
\seealso{
   \code{\link{lmoms}}, \code{\link{dlmomco}}, \code{\link{plmomco}}, \code{\link{rlmomco}}, \code{\link{qlmomco}}, \code{\link{lmom2par}},
\code{\link{plotlmrdia}}, \code{\link{lcomoms2}}
}
\keyword{distribution}
\keyword{mixture}
\keyword{L-moments}
\concept{Lmoments}
\concept{LMoments}
\concept{linear moments}
\concept{L moments}
\concept{L Moments}
\concept{method of L moments}
\concept{method of L-moments}
\concept{method of L Moments}
\concept{method of maximum likelihood}
\concept{method of maximum of product spacings}
\concept{maximum product estimation}
