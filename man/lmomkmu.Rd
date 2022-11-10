\name{lmomkmu}
\alias{lmomkmu}
\title{L-moments of the Kappa-Mu Distribution}
\description{
This function estimates the L-moments of the Kappa-Mu (\eqn{\kappa:\mu})  distribution given the parameters (\eqn{\nu} and \eqn{\alpha}) from \code{\link{parkmu}}. The L-moments in terms of the parameters are complex. They are computed here by the \eqn{\alpha_r} probability-weighted moments in terms of the Marcum Q-function (see \code{\link{cdfkmu}}). The linear combination relating the L-moments to the \eqn{\beta_r} probability-weighted moments is
\deqn{
\lambda_{r+1} = \sum_{k=0}^{r} (-1)^{r-k} {r \choose k} { r + k \choose k } \beta_k
\mbox{,}}
for \eqn{r \ge 0} and the linear combination relating \eqn{\alpha_r} to \eqn{\beta_r} is
\deqn{
\alpha_r = \sum_{k=0}^r (-1)^k { r \choose k } \beta_k
\mbox{,}}
and by definition the \eqn{\alpha_r} are the expectations
\deqn{
\alpha_r \equiv E\{ X\,[1-F(X)]^r\}
\mbox{,}}
and thus
\deqn{
\alpha_r =  \int_{-\infty}^{\infty} x\, [1 - F(x)]^r f(x)\; \mathrm{d}x
\mbox{,}}
in terms of \eqn{x}, the PDF \eqn{f(x)}, and the CDF \eqn{F(x)}. Lastly, the \eqn{\alpha_r} for the Kappa-Mu distribution with substitutions of the Marcum Q-function are
\deqn{
\alpha_r = \int_{-\infty}^{\infty} Q_\mu\biggl(\sqrt{2\kappa\mu},\; x\sqrt{2(1+\kappa)\mu}\biggr)^r\,x\, f(x)\; \mathrm{d}x\mbox{.}
}
Although multiple methods for Marcum Q-function computation are in \code{\link{cdfkmu}} and discussed in that documentation, the \code{lmomkmu} presenting is built only using the \dQuote{chisq} approach.

Yacoub (2007, eq. 5) provides an expectation for the \eqn{j}th moment of the distribution as given by
\deqn{
\mathrm{E}(x^j) = \frac{\Gamma(\mu+j/2)\mathrm{exp}(-\kappa\mu)}{\Gamma(\mu)[(1+\kappa)\mu]^{j/2}} \times {}_1F_1(\mu+j/2; \mu; \kappa\mu)
\mbox{,}}
where \eqn{{}_1F_1(a;b;z)} is the confluent hypergeometric function of Abramowitz and Stegun (1972, eq. 13.1.2). The \code{lmomkmu} function optionally solves for the mean (\eqn{j=1}) using the above equation in conjunction with the mean as computed by the order statistic minimums. The \eqn{{}_1F_1(a;b;z)} is defined as
\deqn{
{}_1F_1(a;b;z) = \sum_{i=0}^\infty \frac{a^{(i)}}{b^{(i)}}\frac{z^i}{n!}
\mbox{,}}
where the notation \eqn{a^{(n)}} represents \dQuote{rising factorials} that are defined as \eqn{a^{(0)} = 1} and \eqn{a^{(n)} = a(a+1)(a+2)\ldots(a+n-1)}. The rising factorials are readily computed by \eqn{ a^{(n)} = \Gamma(n+1)/\Gamma(n)} without resorting to a series computation. Yacoub (2007, eq. 5) is used to compute the mean.
}
\usage{
lmomkmu(para, nmom=5, paracheck=TRUE, tol=1E-6, maxn=100)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{nmom}{The number of moments to compute.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity.}
  \item{tol}{An absolute tolerance term for series convergence of the confluent hypergeometric function when the Yacoub (2007) mean is to be computed.}
  \item{maxn}{The maximum number of interations in the series of the confluent hypergeometric function when the Yacoub (2007) mean is to be computed.}
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
source of the L-moments: \dQuote{lmomkmu}.}
   \item{yacoubsmean}{A list containing the mean, convergence error, and number of iterations in the series until convergence.}
}
\references{
Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81.
}
\author{W.H. Asquith}
\seealso{\code{\link{parkmu}}, \code{\link{cdfkmu}}, \code{\link{pdfkmu}}, \code{\link{quakmu}}
}
\examples{
kmu <- vec2par(c(1.19,2.3), type="kmu")
lmomkmu(kmu)
\dontrun{
par <- vec2par(c(1.67, .5), type="kmu")
lmomkmu(par)$lambdas
cdf2lmoms(par, nmom=4)$lambdas

system.time(lmomkmu(par))
system.time(cdf2lmoms(par, nmom=4))
}
# See the examples under lmomemu() so visualize L-moment
# relations on the L-skew and L-kurtosis diagram
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Kappa-Mu}
