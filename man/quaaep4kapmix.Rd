\name{quaaep4kapmix}
\alias{quaaep4kapmix}
\title{Quantile Function Mixture Between the 4-Parameter Asymmetric Exponential Power and Kappa Distributions}
\description{
This function computes the quantiles of a mixture as needed between the 4-parameter Asymmetric Exponential Power (AEP4) and Kappa distributions given L-moments (\code{\link{lmoms}}). The quantile function of a two-distribution mixture is supported by \code{\link{par2qua2}} and is
\deqn{x(F) = (1-w) \times A(F) + w \times K(F)\mbox{,} }
where \eqn{x(F)} is the mixture for nonexceedance probability \eqn{F}, \eqn{A(F)} is the AEP4 quantile function (\code{\link{quaaep4}}), \eqn{K(F)} is the Kappa quantile function (\code{\link{quakap}}), and \eqn{w} is a weight factor.

Now, the above mixture is only applied if the \eqn{\tau_4} for the given \eqn{\tau_3} is within the overlapping region of the AEP4 and Kappa distributions. For this condition, the \eqn{w} is computed by proration between the upper Kappa distribution bound (same as the \eqn{\tau_3} and \eqn{\tau_4} of the Generalized Logistic distribution, see \code{\link{lmrdia}}) and the lower bounds of the AEP4. For \eqn{\tau_4} above the Kappa, then the AEP4 is exclusive and conversely, for \eqn{\tau_4} below the AEP4, then the Kappa is exclusive.

The \eqn{w} therefore is the proration
\deqn{w = [\tau^{K}_4(\hat\tau_3) - \hat\tau_4] / [\tau^{K}_4(\hat\tau_3) - \tau^{A}_4(\hat\tau_3)]\mbox{,}}
where \eqn{\hat\tau_4} is the sample L-kurtosis, \eqn{\tau^{K}_4} is the upper bounds of the Kappa and \eqn{\tau^{A}_4} is the lower bounds of the AEP4 for the sample L-skew (\eqn{\hat\tau_3}).

The parameter estimation for the AEP4 by \code{\link{paraep4}} can fall back to pure Kappa if argument \code{kapapproved=TRUE} is set. Such a fall back is unrelated to the mixture described here.
}
\usage{
quaaep4kapmix(f, lmom, checklmom=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{lmom}{A L-moment object created by \code{\link{lmoms}} or similar.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.
}
\seealso{\code{\link{par2qua2}}, \code{\link{quaaep4}}, \code{\link{quakap}}, \code{\link{paraep4}}, \code{\link{parkap}}}
\examples{
\dontrun{
FF <- c(0.0001, 0.0005, 0.001, seq(0.01,0.99, by=0.01), 0.999,
       0.9995, 0.9999); Z <- qnorm(FF)
t3s <- seq(0, 0.5, by=0.1); T4step <- 0.02
pdf("mixture_test.pdf")
for(t3 in t3s) {
   T4low <- (5*t3^2 - 1)/4; T4kapup <- (5*t3^2 + 1)/6
   t4s <- seq(T4low+T4step, T4kapup+2*T4step, by=T4step)
   for(t4 in t4s) {
      lmr <- vec2lmom(c(0,1,t3,t4)) # make L-moments for lmomco
      if(! are.lmom.valid(lmr)) next # for general protection
      kap  <- parkap(lmr)
      if(kap$ifail == 5) next # avoid further work if numeric problems
      aep4 <- paraep4(lmr, method="A")
      X <- quaaep4kapmix(FF, lmr)
      if(is.null(X)) next # one last protection
      plot(Z, X, type="l", lwd=5, col=1, ylim=c(-15,15),
           xlab="STANDARD NORMAL VARIATE",
           ylab="VARIABLE VALUE")
      mtext(paste("L-skew =",lmr$ratios[3],
                  "  L-kurtosis = ",lmr$ratios[4]))
      # Now add two more quantile functions for reference and review
      # of the mixture. These of course would not be done in practice
      # only quaaep4kapmix() would suffice.
      if(! as.logical(aep4$ifail)) {
         lines(Z, qlmomco(F,aep4), lwd=2, col=2)
      }
      if(! as.logical(kap$ifail)) {
         lines(Z, qlmomco(F,kap),  lwd=2, col=3)
      }
      message("t3=",t3,"  t4=",t4) # stout for a log file
  }
}
dev.off()
}
}
\keyword{quantile function}
\keyword{quantile mixture function}
\keyword{distribution (mixture)}
\keyword{Distribution: Kappa}
\keyword{Distribution: Asymmetric Exponential Power}




