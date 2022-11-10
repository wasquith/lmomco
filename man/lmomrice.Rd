\name{lmomrice}
\alias{lmomrice}
\title{L-moments of the Rice Distribution}
\description{
This function estimates the L-moments of the Rice distribution given the parameters (\eqn{\nu} and \eqn{\alpha}) from \code{\link{parrice}}. The L-moments in terms of the parameters are complex. They are computed here by the system of maximum order statistic expectations from \code{\link{theoLmoms.max.ostat}}, which uses \code{\link{expect.max.ostat}}. The connection between \eqn{\tau_2} and \eqn{\nu/\alpha} and a special function (the Laguerre polynomial, \code{\link{LaguerreHalf}}) of \eqn{\nu^2/\alpha^2} and additional algebraic terms is tabulated in the \R \code{data.frame} located within \env{.lmomcohash$RiceTable}. The file \file{SysDataBuilder.R} provides additional details.
}
\usage{
lmomrice(para, ...)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{...}{Additional arguments passed to \code{\link{theoLmoms.max.ostat}}.}
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
\item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomrice}, but the exact contents of the remainder of the string might vary as limiting distributions of Normal and Rayleigh can be involved for \eqn{\nu/\alpha > 52} (super high SNR, Normal) or \eqn{24 < \nu/\alpha \le 52} (high SNR, Normal) or \eqn{\nu/\alpha < 0.08} (very low SNR, Rayleigh).}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\author{W.H. Asquith}
\seealso{\code{\link{parrice}}, \code{\link{cdfrice}}, \code{\link{cdfrice}}, \code{\link{quarice}} }
\examples{
\dontrun{
lmomrice(vec2par(c(65,34), type="rice"))

# Use the additional arguments to show how to avoid unnecessary overhead
# when using the Rice, which only has two parameters.
  rice <- vec2par(c(15,14), type="rice")
  system.time(lmomrice(rice, nmom=2)); system.time(lmomrice(rice, nmom=6))

  lcvs <- vector(mode="numeric"); i <- 0
  SNR  <- c(seq(7,0.25, by=-0.25), 0.1)
  for(snr in SNR) {
    i <- i + 1
    rice    <- vec2par(c(10,10/snr), type="rice")
    lcvs[i] <- lmomrice(rice, nmom=2)$ratios[2]
  }
  plot(lcvs, SNR,
       xlab="COEFFICIENT OF L-VARIATION",
       ylab="LOCAL SIGNAL TO NOISE RATIO (NU/ALPHA)")
  lines(.lmomcohash$RiceTable$LCV,
        .lmomcohash$RiceTable$SNR)
  abline(1,0, lty=2)
  mtext("Rice Distribution")
  text(0.15,0.5, "More noise than signal")
  text(0.15,1.5, "More signal than noise")
}
\dontrun{
# A polynomial expression for the relation between L-skew and 
# L-kurtosis for the Rice distribution can be readily constructed.
T3 <- .lmomcohash$RiceTable$TAU3
T4 <- .lmomcohash$RiceTable$TAU4
LM <- lm(T4~T3+I(T3^2)+I(T3^3)+I(T3^4)+
               I(T3^5)+I(T3^6)+I(T3^7)+I(T3^8))
summary(LM) # note shown}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Rice}
