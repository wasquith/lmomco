\name{lmomkur}
\alias{lmomkur}
\title{L-moments of the Kumaraswamy Distribution}
\description{
This function estimates the L-moments of the Kumaraswamy distribution given the parameters (\eqn{\alpha} and \eqn{\beta}) from \code{\link{parkur}}. The L-moments in terms of the parameters with \eqn{\eta = 1 + 1/\alpha} are
\deqn{\lambda_1 = \beta B(\eta, \beta) \mbox{,}}
\deqn{\lambda_2 = \beta [B(\eta, \beta) - 2B(\eta, 2\beta)] \mbox{,}}
\deqn{\tau_3 = \frac{B(\eta,\beta) - 6B(\eta,2\beta) + 6B(\eta,3\beta)}{B(\eta,\beta) - 2B(\eta,2\beta)} \mbox{,}}
\deqn{\tau_4 = \frac{B(\eta,\beta) - 12B(\eta,2\beta) + 30B(\eta,3\beta) - 40B(\eta,4\beta)}{B(\eta,\beta) - 2B(\eta,2\beta)} \mbox{, and}}
\deqn{\tau_5 = \frac{B(\eta,\beta) - 20B(\eta,2\beta) + 90B(\eta,3\beta) - 140B(\eta,4\beta) + 70B(\eta,5\beta)}{B(\eta,\beta) - 2B(\eta,2\beta)} \mbox{.}}
where \eqn{B(a,b)} is the complete beta function or \code{beta()}.
}
\usage{
lmomkur(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
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
source of the L-moments: \dQuote{lmomkur}.}
}
\references{
Jones, M.C., 2009, Kumaraswamy's distribution---A beta-type distribution with some tractability advantages: Statistical Methodology, v. 6, pp. 70--81.
}
\author{W.H. Asquith}
\seealso{\code{\link{parkur}},  \code{\link{cdfkur}}, \code{\link{pdfkur}}, \code{\link{quakur}} }
\examples{
lmr <- lmoms(c(0.25, 0.4, 0.6, 0.65, 0.67, 0.9))
lmomkur(parkur(lmr))
\dontrun{
A <- B <- exp(seq(-3,5, by=.05))
logA <- logB <- T3 <- T4 <- c();
i <- 0
for(a in A) {
  for(b in B) {
    i <- i + 1
    parkur <- list(para=c(a,b), type="kur");
    lmr <- lmomkur(parkur)
    logA[i] <- log(a); logB[i] <- log(b)
    T3[i] <- lmr$ratios[3]; T4[i] <- lmr$ratios[4]
  }
}
library(lattice)
contourplot(T3~logA+logB, cuts=20, lwd=0.5, label.style="align",
            xlab="LOG OF ALPHA", ylab="LOG OF BETA",
            xlim=c(-3,5), ylim=c(-3,5),
            main="L-SKEW FOR KUMARASWAMY DISTRIBUTION")
contourplot(T4~logA+logB, cuts=10, lwd=0.5, label.style="align",
            xlab="LOG OF ALPHA", ylab="LOG OF BETA",
            xlim=c(-3,5), ylim=c(-3,5),
            main="L-KURTOSIS FOR KUMARASWAMY DISTRIBUTION")
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Kumaraswamy}
