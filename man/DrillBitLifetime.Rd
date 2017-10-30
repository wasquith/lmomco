\name{DrillBitLifetime}
\docType{data}
\alias{DrillBitLifetime}
\title{Lifetime of Drill Bits}
\description{
Hamada (1995, table 9.3) provides a table of lifetime to breakage measured in cycles for drill bits used for producing small holes in printed circuit boards. The data were collected under various control and noise factors to perform reliability assessment to maximize bit reliability with minimization of hole diameter. Smaller holes permit higher density of placed circuitry, and are thus economically attractive. The testing was completed at 3,000 cycles---the right censoring threhold.
}
\usage{
data(DrillBitLifetime)
}
\format{
A data frame with
\describe{
  \item{LIFETIME}{Measured in cycles.}
}
}
\references{
Hamada, M., 1995, Analysis of experiments for reliability improvement and robust reliability: in Balakrishnan, N. (ed.) Recent Advances in Life-Testing and Reliability: Boca Raton, Fla., CRC Press, ISBN 0--8493--8972--0, pp. 155--172.
}
\examples{
data(DrillBitLifetime)
summary(DrillBitLifetime)
\dontrun{
data(DrillBitLifetime)
X     <- DrillBitLifetime$LIFETIME
lmr   <- lmoms(X); par <- lmom2par(lmr,  type="gpa")
pwm   <- pwmRC(X, threshold=3000); zeta <- pwm$zeta
lmrrc <- pwm2lmom(pwm$Bbetas)
rcpar <- pargpaRC(lmrrc, zeta=zeta)
XBAR  <- lmomgpa(rcpar)$lambdas[1]
F <- nonexceeds(); P <- 100*F; x <- seq(min(X), max(X))
plot(sort(X), 100*pp(X), xlab="LIFETIME", ylab="PERCENT", xlim=c(1,10000))
rug(X, col=rgb(0,0,0,0.5))
lines(c(XBAR, XBAR), range(P), lty=2) # mean (expectation of life)
lines(cmlmomco(F, rcpar),  P,  lty=2) # conditional mean
points(XBAR, 0, pch=16)
lines(x, 100*plmomco(x, par),   lwd=2, col=8) # fitted dist.
lines(x, 100*plmomco(x, rcpar), lwd=3, col=1) # fitted dist.

lines( rmlmomco(F, rcpar), P,   col=4) # residual mean life
lines(rrmlmomco(F, rcpar), P,   col=4, lty=2) # rev. residual mean life
lines(x, 1E4*hlmomco(x, rcpar), col=2) # hazard function
lines(x, 1E2*lrzlmomco(plmomco(x, rcpar), rcpar), col=3) # Lorenz func.
legend(4000, 40,
       c("Mean (vertical) or conditional mean (dot at intersect.)",
         "Fitted GPA naively to all data",
         "Fitted GPA to right-censoring PWMs",
         "Residual mean life", "Reversed residual mean life",
         "Hazard function x 1E4", "Lorenz curve x 100"
        ), cex=0.75,
       lwd=c(1, 2, 3, 1, 1, 1, 1), col=c(1, 8, 1, 4, 4, 2, 3),
       lty=c(2, 1, 1, 1, 2, 1, 1), pch=rep(NA, 8))
}
}
\keyword{datasets}
