\name{reslife.lmoms}
\alias{reslife.lmoms}
\title{L-moments of Residual Life}
\description{
This function computes the L-moments of residual life for a quantile function \eqn{x(F)} for an exceedance threshold in probabiliy of \eqn{u}. The L-moments of residual life are thoroughly described by Nair et al. (2013, p. 202).  These L-moments are define as
\deqn{\lambda(u)_r = \sum_{k=0}^{r-1} (-1)^k {r-1 \choose k}^2 \int_u^1 \left(\frac{p-u}{1-u}\right)^{r-k-1} \left(\frac{1-p}{1-u}\right)^k \frac{x(p)}{1-u}\,\mathrm{d}p \mbox{,}}
where \eqn{\lambda(u)_r} is the \eqn{r}th L-moment at residual life probability \eqn{u}.  The L-moment ratios \eqn{\tau(u)_r} have the usual definitions. The implementation here exclusively uses the quantile function of the distribution. If \eqn{u=0}, then the usual L-moments of the quantile function are returned because the integration domain is the entire potential lifetime range.  If \eqn{u=1}, then \eqn{\lambda(1)_1 = x(1)} is returned, which is the maximum lifetime of the distribution (the value for the upper support of the distribution), and the remaining \eqn{\lambda(1)_r} for \eqn{r \ge 2} are set to \code{NA}. Lastly, the notation \eqn{(u)} is neither super or subscripted to avoid confusion with L-moment order \eqn{r} or the TL-moments that indicate trimming level as a superscript (see \code{\link{TLmoms}}).
}
\usage{
reslife.lmoms(f, para, nmom=5)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{life.exceeds}{The value for \eqn{x(F)} for \eqn{F=} \code{f}.}
  \item{life.percentile}{The value \eqn{100\times}\code{f}.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{NULL} because no trimming theory for L-moments of residual life have been developed or researched.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL} because no trimming theory for L-moments of residual life have been developed or researched.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL} because no trimming theory for L-moments of residual life have been developed or researched.}
  \item{source}{An attribute identifying the computational source of the L-moments: \cr \dQuote{reslife.lmoms}.}
}
\author{W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{rmlmomco}}, \code{\link{rreslife.lmoms}}
}
\examples{
A <- vec2par(c(230, 2649, 3), type="gov") # Set lower bounds = 230 hours
F <- nonexceeds(f01=TRUE)
plot(F, rmlmomco(F,A), type="l", ylim=c(0,3000), # mean residual life [black]
     xlab="NONEXCEEDANCE PROBABILITY",
     ylab="LIFE, RESIDUAL LIFE (RL), RL_L-SCALE, RL_L-skew (rescaled)")
L1 <- L2 <- T3 <- vector(mode="numeric", length=length(F))
for(i in 1:length(F)) {
  lmr <- reslife.lmoms(F[i], A, nmom=3)
  L1[i] <- lmr$lambdas[1]; L2[i] <- lmr$lambdas[2]; T3[i] <- lmr$ratios[3]
}
lines(c(0,1), c(1500,1500),  lty=2) # Origin line (to highlight T3 crossing "zero")
lines(F, L1,          col=2, lwd=3) # Mean life (not residual, that is M(u)) [red]
lines(F, L2,          col=3, lwd=3) # L-scale of residual life [green]
lines(F, 5E3*T3+1500, col=4, lwd=3) # L-skew of residual life (re-scaled) [blue]
\dontrun{
# Nair et al. (2013, p. 203), test shows L2(u=0.37) = 771.2815
A <- vec2par(c(230, 2649, 0.3), type="gpa"); F <- 0.37
"afunc" <- function(p) { return((1-p)*rmlmomco(p,A)) }
L2u1 <- (1-F)^(-2)*integrate(afunc,F,1)$value
L2u2 <- reslife.lmoms(F,A)$lambdas[2]
}
}
\keyword{lifetime/reliability analysis}
\keyword{L-moment (residual life)}
\keyword{L-moment (theoretical)}
\concept{residual life}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
