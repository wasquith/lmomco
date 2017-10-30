\name{theoLmoms.max.ostat}
\alias{theoLmoms.max.ostat}
\alias{theoLmoms.min.ostat}
\title{Compute the Theoretical L-moments of a Distribution Distribution based on System of Maximum Order Statistic Expectations}
\description{
This function computes the theoretical L-moments of a distribution by the following
\deqn{
\lambda_r = (-1)^{r-1} \sum_{k=1}^r (-1)^{r-k}k^{-1}{r-1 \choose k-1}{r+k-2 \choose k-1}\mathrm{E}[X_{1:k}]
}
for the minima  (\code{\link{theoLmoms.min.ostat}}, theoretical L-moments from the minima of order statistics) or
\deqn{
\lambda_r = \sum_{k=1}^r (-1)^{r-k}k^{-1}{r-1 \choose k-1}{r+k-2 \choose k-1}\mathrm{E}[X_{k:k}]
}
for the maxima (\code{\link{theoLmoms.max.ostat}}, theoretical L-moments from the maxima of order statistics). The functions \code{\link{expect.min.ostat}} and \code{\link{expect.max.ostat}} compute the minima (\eqn{\mathrm{E}[X_{1:k}]}) and maxima (\eqn{\mathrm{E}[X_{k:k}]}), respectively.

If \code{qua != NULL}, then the first expectation equation shown under \code{\link{expect.max.ostat}} is used for the order statistic expectations and any function set in \code{cdf} and \code{pdf} is ignored.
}
\usage{
theoLmoms.max.ostat(para=NULL, cdf=NULL, pdf=NULL, qua=NULL,
                    nmom=4, switch2minostat=FALSE, showterms=FALSE, ...)
}
\arguments{
  \item{para}{A distribution parameter list from a function such as \code{lmom2par} or \code{\link{vec2par}}.}
  \item{cdf}{CDF of the distribution for the parameters.}
  \item{pdf}{PDF of the distribution for the parameters.}
  \item{qua}{Quantile function for the parameters.}
  \item{nmom}{The number of L-moments to compute.}
  \item{switch2minostat}{A logical in which a switch to the expectations of minimum order statistics will be used and \code{\link{expect.min.ostat}} instead of \code{\link{expect.max.ostat}} will be used with expected small change in overall numerics. The function \cr \code{\link{theoLmoms.min.ostat}} provides a direct interface for L-moment computation by minimum order statistics.}
  \item{showterms}{A logical controlling just a reference message that will show the multipliers on each of the order statistic minima or maxima that comprise the terms within the summations in the above formulae (see Asquith, 2011, p. 95).}
  \item{...}{Optional, but likely, arguments to pass to \code{\link{expect.min.ostat}} or \cr\code{\link{expect.max.ostat}}. Such arguments will likely tailor the integration limits that can be specific for the distribution in question. Further these arguments might be needed for the cumulative distribution function.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments: first element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal \code{NULL} until trimming support is made.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which will equal \code{NULL} until trimming support is made.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which will equal \code{NULL} until trimming support is made.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{theoLmoms.max.ostat}.}
}
\note{
Perhaps one of the neater capabilities that the \code{theoLmoms.max.ostat} and \code{theoLmoms.min.ostat} functions provide is for computing L-moments that are not analytically available from other authors or have no analytical solution.
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\author{ W.H. Asquith}
\seealso{\code{\link{theoLmoms}}, \code{\link{expect.min.ostat}}, \code{\link{expect.max.ostat}}}
\examples{
\dontrun{
para <- vec2par(c(40,20), type='nor')
A1 <- theoLmoms.max.ostat(para=para, cdf=cdfnor, pdf=pdfnor, switch2minostat=FALSE)
A2 <- theoLmoms.max.ostat(para=para, cdf=cdfnor, pdf=pdfnor, switch2minostat=TRUE)
B1 <- theoLmoms.max.ostat(para=para, qua=quanor, switch2minostat=FALSE)
B2 <- theoLmoms.max.ostat(para=para, qua=quanor, switch2minostat=TRUE)
print(A1$ratios[4]) # reports 0.1226017
print(A2$ratios[4]) # reports 0.1226017
print(B1$ratios[4]) # reports 0.1226012
print(B2$ratios[4]) # reports 0.1226012
# Theoretical value = 0.122601719540891.
# Confirm operational with native R-code being used inside lmomco functions
# Symmetrically correct on whether minima or maxima are used, but some
# Slight change when qnorm() used instead of dnorm() and pnorm().

para <- vec2par(c(40,20), type='exp')
A1 <- theoLmoms.max.ostat(para=para, cdf=cdfexp, pdf=pdfexp, switch2minostat=FALSE)
A2 <- theoLmoms.max.ostat(para=para, cdf=cdfexp, pdf=pdfexp, switch2minostat=TRUE)
B1 <- theoLmoms.max.ostat(para=para, qua=quaexp, switch2minostat=FALSE)
B2 <- theoLmoms.max.ostat(para=para, qua=quaexp, switch2minostat=TRUE)
print(A1$ratios[4]) # 0.1666089
print(A2$ratios[4]) # 0.1666209
print(B1$ratios[4]) # 0.1666667
print(B2$ratios[4]) # 0.1666646
# Theoretical value = 0.1666667

para <- vec2par(c(40,20), type='ray')
A1 <- theoLmoms.max.ostat(para=para, cdf=cdfray, pdf=pdfray, switch2minostat=FALSE)
A2 <- theoLmoms.max.ostat(para=para, cdf=cdfray, pdf=pdfray, switch2minostat=TRUE)
B1 <- theoLmoms.max.ostat(para=para, qua=quaray, switch2minostat=FALSE)
B2 <- theoLmoms.max.ostat(para=para, qua=quaray, switch2minostat=TRUE)
print(A1$ratios[4]) # 0.1053695
print(A2$ratios[4]) # 0.1053695
print(B1$ratios[4]) # 0.1053636
print(B2$ratios[4]) # 0.1053743
# Theoretical value = 0.1053695
}
\dontrun{
# The Rice distribution is complex and tailoring of the integration
# limits is needed to effectively trap errors, the limits for the
# Normal distribution above are infinite so no granular control is needed.
para <- vec2par(c(30,10), type="rice")
theoLmoms.max.ostat(para=para, cdf=cdfrice, pdf=pdfrice,
                    lower=0, upper=.Machine$double.max)
}
\dontrun{
para <- vec2par(c(0.6, 1.5), type="emu")
theoLmoms.min.ostat(para, cdf=cdfemu, pdf=pdfemu,
                    lower=0, upper=.Machine$double.max)
theoLmoms.min.ostat(para, cdf=cdfemu, pdf=pdfemu, yacoubsintegral = FALSE,
                    lower=0, upper=.Machine$double.max)

para <- vec2par(c(0.6, 1.5), type="kmu")
theoLmoms.min.ostat(para, cdf=cdfkmu, pdf=pdfkmu,
                    lower=0, upper=.Machine$double.max)
theoLmoms.min.ostat(para, cdf=cdfkmu, pdf=pdfkmu, marcumQ = FALSE,
                    lower=0, upper=.Machine$double.max)
}
\dontrun{
# The Normal distribution is used on the fly for the Rice for high to
# noise ratios (SNR=nu/alpha > some threshold). This example will error out.
nu <- 30; alpha <- 0.5
para <- vec2par(c(nu,alpha), type="rice")
theoLmoms.max.ostat(para=para, cdf=cdfrice, pdf=pdfrice,
                    lower=0, upper=.Machine$double.max)
}
}
\keyword{L-moment (theoretical)}
\keyword{order statistics (misc.)}
