\name{lmomemu}
\alias{lmomemu}
\title{L-moments of the Eta-Mu Distribution}
\description{
This function estimates the L-moments of the Eta-Mu (\eqn{\eta:\mu}) distribution given the parameters (\eqn{\eta} and \eqn{\mu}) from \code{\link{paremu}}. The L-moments in terms of the parameters are complex.  They are computed here by the \eqn{\alpha_r} probability-weighted moments in terms of the Yacoub integral (see \code{\link{cdfemu}}). The linear combination relating the L-moments to the conventional \eqn{\beta_r} probability-weighted moments is
\deqn{
\lambda_{r+1} = \sum_{k=0}^{r} (-1)^{r-k} {r \choose k} { r + k \choose k } \beta_k\mbox{,}
}
for \eqn{r \ge 0} and the linear combination relating the less common \eqn{\alpha_r} to \eqn{\beta_r} is
\deqn{
\alpha_r = \sum_{k=0}^r (-1)^k { r \choose k } \beta_k\mbox{,}
}
and by definition the \eqn{\alpha_r} are the expectations
\deqn{
\alpha_r \equiv E\{ X\,[1-F(X)]^r\}\mbox{,}
}
and thus
\deqn{
\alpha_r =  \int_{-\infty}^{\infty} x\, [1 - F(x)]^r f(x)\; \mathrm{d}x\mbox{,}
}
in terms of \eqn{x}, the PDF \eqn{f(x)}, and the CDF \eqn{F(x)}. Lastly, the \eqn{\alpha_r} for the Eta-Mu distribution with substitution of the Yacoub integral are
\deqn{
\alpha_r = \int_{-\infty}^{\infty} Y_\mu\biggl( \eta,\; x\sqrt{2h\mu} \biggr)^r\,x\, f(x)\; \mathrm{d}x\mbox{.}
}

Yacoub (2007, eq. 21) provides an expectation for the \eqn{j}th moment of the distribution as given by
\deqn{
\mathrm{E}(x^j) = \frac{\Gamma(2\mu+j/2)}{h^{\mu+j/2}(2\mu)^{j/2}\Gamma(2\mu)}\times {}_2F_1(\mu+j/4+1/2, \mu+j/4; \mu+1/2; (H/h)^2)\mbox{,}
}
where \eqn{{}_2F_1(a,b;c;z)} is the Gauss hypergeometric function of Abramowitz and Stegun (1972, eq. 15.1.1) and \eqn{h = 1/(1-\eta^2)} (format 2 of Yacoub's paper and the format exclusively used by \pkg{lmomco}). The \code{lmomemu} function optionally solves for the mean (\eqn{j=1}) using the above equation in conjunction with the mean as computed by the order statistic minimums. The \eqn{{}_2F_1(a,b;c;z)} is defined as
\deqn{
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)}{\Gamma(a)\Gamma{(b)}} \sum_{i=0}^\infty \frac{\Gamma(a+i)\Gamma{(b+i)}}{\Gamma{(c+i)}}\frac{z^i}{n!}\mbox{.}
}
Yacoub (2007, eq. 21) is used to compute the mean.
}
\usage{
lmomemu(para, nmom=5, paracheck=TRUE, tol=1E-6, maxn=100)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{nmom}{The number of L-moments to compute.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity.}
  \item{tol}{An absolute tolerance term for series convergence of the Gauss hypergeometric function when the Yacoub (2007) mean is to be computed.}
  \item{maxn}{The maximum number of interations in the series of the Gauss hypergeometric function when the Yacoub (2007) mean is to be computed.}
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
source of the L-moments: \dQuote{lmomemu}.}
   \item{yacoubsmean}{A list containing the mean, convergence error, and number of iterations in the series until convergence.}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81
}
\author{W.H. Asquith}
\seealso{\code{\link{paremu}}, \code{\link{cdfemu}}, \code{\link{pdfemu}}, \code{\link{quaemu}}  }
\examples{
\dontrun{
emu <- vec2par(c(.19,2.3), type="emu")
lmomemu(emu)

par <- vec2par(c(.67, .5), type="emu")
lmomemu(par)$lambdas
cdf2lmoms(par, nmom=4)$lambdas
system.time(lmomemu(par))
system.time(cdf2lmoms(par, nmom=4))

# This extensive sequence of operations provides very important
# perspective on the L-moment ratio diagram of L-skew and L-kurtosis.
# But more importantly this example demonstrates the L-moment
# domain of the Kappa-Mu and Eta-Mu distributions and their boundaries.
#
t3 <- seq(-1,1,by=.0001)
plotlmrdia(lmrdia(), xlim=c(-0.05,0.5), ylim=c(-0.05,.2))
# The following polynomials are used to define the boundaries of
# both distributions. The applicable inequalities for these
# are not provided for these polynomials as would be in deeper
# implementation---so don't worry about wild looking trajectories.
"KMUup" <- function(t3) {
             return(0.1227 - 0.004433*t3 - 2.845*t3^2 +
                    + 18.41*t3^3 - 50.08*t3^4 + 83.14*t3^5 +
                    - 81.38*t3^6 + 43.24*t3^7 - 9.600*t3^8)}

"KMUdnA" <- function(t3) {
              return(0.1226 - 0.3206*t3 - 102.4*t3^2 - 4.753E4*t3^3 +
                     - 7.605E6*t3^4 - 5.244E8*t3^5 - 1.336E10*t3^6)}

"KMUdnB" <- function(t3) {
              return(0.09328 - 1.488*t3 + 16.29*t3^2 - 205.4*t3^3 +
                     + 1545*t3^4 - 5595*t3^5 + 7726*t3^6)}

"KMUdnC" <- function(t3) {
              return(0.07245 - 0.8631*t3 + 2.031*t3^2 - 0.01952*t3^3 +
                     - 0.7532*t3^4 + 0.7093*t3^5 - 0.2156*t3^6)}

"EMUup" <- function(t3) {
              return(0.1229 - 0.03548*t3 - 0.1835*t3^2 + 2.524*t3^3 +
                     - 2.954*t3^4 + 2.001*t3^5 - 0.4746*t3^6)}

# Here, we are drawing the trajectories of the tabulated parameters
# and L-moments within the internal storage of lmomco.
lines(.lmomcohash$EMU_lmompara_byeta$T3,
      .lmomcohash$EMU_lmompara_byeta$T4,   col=7, lwd=0.5)
lines(.lmomcohash$KMU_lmompara_bykappa$T3,
      .lmomcohash$KMU_lmompara_bykappa$T4, col=8, lwd=0.5)

# Draw the polynomials
lines(t3, KMUdnA(t3), lwd=4, col=2, lty=4)
lines(t3, KMUdnB(t3), lwd=4, col=3, lty=4)
lines(t3, KMUdnC(t3), lwd=4, col=4, lty=4)
lines(t3, EMUup(t3),  lwd=4, col=5, lty=4)
lines(t3, KMUup(t3),  lwd=4, col=6, lty=4)
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Eta-Mu}
