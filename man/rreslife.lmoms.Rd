\name{rreslife.lmoms}
\alias{rreslife.lmoms}
\title{L-moments of Reversed Residual Life}
\description{
This function computes the L-moments of reversed residual life for a quantile function \eqn{x(F)} for an exceedance threshold in probabiliy of \eqn{u}. The L-moments of residual life are thoroughly described by Nair et al. (2013, p. 211).  These L-moments are define as
\deqn{{}_\mathrm{r}\lambda(u)_r = \sum_{k=0}^{r-1} (-1)^k {r-1 \choose k}^2 \int_0^u \left(\frac{p}{u}\right)^{r-k-1} \left(1 - \frac{p}{u}\right)^k \frac{x(p)}{u}\,\mathrm{d}p \mbox{,}}
where \eqn{{}_\mathrm{r}\lambda(u)_r} is the \eqn{r}th L-moment at residual life probability \eqn{u}.  The L-moment ratios \eqn{{}_\mathrm{r}\tau(u)_r} have the usual definitions. The implementation here exclusively uses the quantile function of the distribution. If \eqn{u=0}, then the usual L-moments of the quantile function are returned because the integration domain is the entire potential lifetime range.  If \eqn{u=0}, then \eqn{{}_\mathrm{r}\lambda(1)_1 = x(0)} is returned, which is the minimum lifetime of the distribution (the value for the lower support of the distribution), and the remaining \eqn{{}_\mathrm{r}\lambda(1)_r} for \eqn{r \ge 2} are set to \code{NA}. The reversal aspect is denoted by the prepended romanscript \eqn{\mathrm{r}} to the \eqn{\lambda}'s and \eqn{\tau}'s. Lastly, the notation \eqn{(u)} is neither super or subscripted to avoid confusion with L-moment order \eqn{r} or the TL-moments that indicate trimming level as a superscript (see \code{\link{TLmoms}}).
}
\usage{
rreslife.lmoms(f, para, nmom=5)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{{}_\mathrm{r}\lambda_1}, second element is \eqn{{}_\mathrm{r}\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{{}_\mathrm{r}\tau}, third element is \eqn{{}_\mathrm{r}\tau_3} and so on.}
  \item{life.notexceeds}{The value for \eqn{x(F)} for \eqn{F=} \code{f}.}
  \item{life.percentile}{The value \eqn{100\times}\code{f}.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{NULL} because no trimming theory for L-moments of residual life have been developed or researched.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL} because no trimming theory for L-moments of residual life have been developed or researched.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL} because no trimming theory for L-moments of residual life have been developed or researched.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{rreslife.lmoms}.}
}
\author{W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{rmlmomco}}, \code{\link{reslife.lmoms}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0
"afunc" <- function(p)        { return(par2qua(p,A,paracheck=FALSE)) }
"bfunc" <- function(p,u=NULL) { return((2*p - u)*par2qua(p,A,paracheck=FALSE)) }
f <- 0.35
rL1a <- integrate(afunc, lower=0, upper=f)$value      / f   # Nair et al. (2013, eq. 6.18)
rL2a <- integrate(bfunc, lower=0, upper=f, u=f)$value / f^2 # Nair et al. (2013, eq. 6.19)
rL <- rreslife.lmoms(f, A, nmom=2) # The data.frame shows equality of the two approaches.
rL1b <- rL$lambdas[1]; rL2b <- rL$lambdas[2]
print(data.frame(rL1a=rL1a, rL1b=rL1b, rL2b=rL2b, rL2b=rL2b))
\dontrun{
# 2nd Example, let us look at Tau3, each of the L-skews are the same.
T3    <- par2lmom(A)$ratios[3]
T3.0  <-  reslife.lmoms(0, A)$ratios[3]
rT3.1 <- rreslife.lmoms(1, A)$ratios[3]
}
\dontrun{
# Nair et al. (2013, p. 212), test shows rL2(u=0.77) = 12.6034
A <- vec2par(c(230, 269, 3.3), type="gpa"); F <- 0.77
"afunc" <- function(p) { return(p*rrmlmomco(p,A)) }
rL2u1 <- (F)^(-2)*integrate(afunc,0,F)$value
rL2u2 <- rreslife.lmoms(F,A)$lambdas[2]
}
}
\keyword{lifetime/reliability analysis}
\keyword{L-moment (residual life)}
\keyword{L-moment (theoretical)}
\concept{residual life}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
