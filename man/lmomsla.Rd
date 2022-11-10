\name{lmomsla}
\alias{lmomsla}
\title{Trimmed L-moments of the Slash Distribution}
\description{
This function estimates the trimmed L-moments of the Slash distribution given the parameters (\eqn{\xi} and \eqn{\alpha}) from \code{\link{parsla}}. The relation between the TL-moments (\code{trim=1}) and the parameters have been numerically determined and are 
\eqn{\lambda^{(1)}_1 = \xi}, 
\eqn{\lambda^{(1)}_2 = 0.9368627\alpha},  
\eqn{\tau^{(1)}_3 = 0}, 
\eqn{\tau^{(1)}_4 = 0.3042045}, 
\eqn{\tau^{(1)}_5 = 0}, and 
\eqn{\tau^{(1)}_6 = 0.1890072}.
These TL-moments (trim=1) are symmetrical for the first L-moments defined because \eqn{\mathrm{E}[X_{1:n}]} and \eqn{\mathrm{E}[X_{n:n}]} are undefined expectations for the Slash.
}
\usage{
lmomsla(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the trimmed L-moments. First element is \eqn{\lambda^{(1)}_1}, second element is \eqn{\lambda^{(1)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is \eqn{\tau^{(1)}}, third element is \eqn{\tau^{(1)}_3} and so on. }
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{1}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{1}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{1}.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomsla}}
  \item{trim}{Level of symmetrical trimming used.}
}
\references{
Rogers, W.H., and Tukey, J.W., 1972, Understanding some long-tailed symmetrical distributions: Statistica Neerlandica, v. 26, no. 3, pp. 211--226.
}
\author{W.H. Asquith}
\seealso{\code{\link{parsla}}, \code{\link{cdfsla}}, \code{\link{pdfsla}}, \code{\link{quasla}} }
\examples{
\dontrun{
# This example was used to numerically back into the TL-moments and the 
# relation between \alpha and \lambda_2.
"lmomtrim1" <- function(para) {
    bigF <- 0.999
    minX <- para$para[1] - para$para[2]*qnorm(1 - bigF) / qunif(1 - bigF)
    maxX <- para$para[1] + para$para[2]*qnorm(    bigF) / qunif(1 - bigF)
    minF <- cdfsla(minX, para); maxF <- cdfsla(maxX, para)
    lmr <- theoTLmoms(para, nmom = 6, leftrim = 1, rightrim = 1)
}

U <- -10; i <- 0
As <- seq(.1,abs(10),by=.2)
L1s <- L2s <- T3s <- T4s <- T5s <- T6s <- vector(mode="numeric", length=length(As))
for(A in As) {
   i <- i + 1
   lmr <- lmomtrim1(vec2par(c(U, A), type="sla"))
   L1s[i] <- lmr$lambdas[1]; L2s[i] <- lmr$lambdas[2]
   T3s[i] <- lmr$ratios[3];  T4s[i] <- lmr$ratios[4]
   T5s[i] <- lmr$ratios[5];  T6s[i] <- lmr$ratios[6]
}
print(summary(lm(L2s~As-1))$coe)
print(mean(T4s))
print(mean(T6s))
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Slash}
