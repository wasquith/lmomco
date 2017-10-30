\name{expect.max.ostat}
\alias{expect.max.ostat}
\alias{expect.min.ostat}
\alias{eostat}
\title{Compute the Expectation of a Maximum (or Minimum and others) Order Statistic}
\description{
The maximum (or minimum) expectation of an order statistic can be directly used for L-moment computation through either of the following two equations (Hosking, 2006) as dictated by using the maximum (\eqn{\mathrm{E}[X_{k:k}]}, \code{\link{expect.max.ostat}}) or minimum (\eqn{\mathrm{E}[X_{1:k}]}, \code{\link{expect.min.ostat}}):
\deqn{
\lambda_r = (-1)^{r-1} \sum_{k=1}^r (-1)^{r-k}k^{-1}{r-1 \choose k-1}{r+k-2 \choose k-1}\mathrm{E}[X_{1:k}]\mbox{,}
}
and
\deqn{
\lambda_r = \sum_{k=1}^r (-1)^{r-k}k^{-1}{r-1 \choose k-1}{r+k-2 \choose k-1}\mathrm{E}[X_{k:k}]\mbox{.}
}

In terms of the quantile function \code{\link{qlmomco}}, the expectation of an order statistic (Asquith, 2011, p. 49) is
\deqn{
  \mathrm{E}[X_{j:n}] = n {n-1 \choose j - 1}\int^1_0 \! x(F)\times F^{j-1} \times (1-F)^{n-j}\; \mathrm{d}F\mbox{,}
}
where \eqn{x(F)} is the quantile function, \eqn{F} is nonexceedance probability, \eqn{n} is sample size, and \eqn{j} is the \eqn{j}th order statistic.

In terms of the probability density function (PDF) \code{\link{dlmomco}} and cumulative density function (CDF) \code{\link{plmomco}}, the expectation of an order statistic (Asquith, 2011, p. 50) is
\deqn{
\mathrm{E}[X_{j:n}] = \frac{1}{\mathrm{B}(j,n-j+1)}\int_{-\infty}^{\infty} [F(x)]^{j-1}[1-F(x)]^{n-j} x\, f(x)\;\mathrm{d} x\mbox{,}
}
where \eqn{F(x)} is the CDF, \eqn{f(x)} is the PDF, and \eqn{\mathrm{B}(j,n-j+1)} is the complete Beta function, which in \R is \code{beta} with the same argument order as shown above.
}
\usage{
expect.max.ostat(n, para=NULL, cdf=NULL, pdf=NULL, qua=NULL,
                 j=NULL, lower=-Inf, upper=Inf, aslist=FALSE, ...)
}
\arguments{
  \item{n}{The sample size.}
  \item{para}{A distribution parameter list from a function such as \code{\link{vec2par}} or \code{\link{lmom2par}}.}
  \item{cdf}{cumulative distribution function of the distribution.}
  \item{pdf}{probability density function of the distribution.}
  \item{qua}{quantile function of the distribution. If this is defined, then \code{cdf} and \code{pdf} are ignored.}
  \item{j}{The \eqn{j}th value of the order statistic, which defaults to \var{n=j} (the maximum order statistic) if \code{j=NULL}.}
  \item{lower}{The lower limit for integration.}
  \item{upper}{The upper limit for integration.}
  \item{aslist}{A logically triggering whether an \R \code{list} is returned instead of just the expection.}
  \item{...}{Additional arguments to pass to the three distribution functions.}
}
\details{
If \code{qua != NULL}, then the first order-statistic expectation equation above is used, and any function that might have been set in \code{cdf} and \code{pdf} is \emph{ignored}. If the limits are infinite (default), then the limits of the integration will be set to \eqn{F\!\downarrow = 0} and  \eqn{F\!\uparrow = 1}. The user can replace these by setting the limits to something \dQuote{near} zero and(or) \dQuote{near} 1. Please consult the Note below concerning more information about the limits of integration.

If \code{qua == NULL}, then the second order-statistic expectation equation above is used and \code{cdf} and \code{pdf} must be set. The default \eqn{\pm\infty} limits are used unless the user \emph{knows} otherwise for the distribution or through supervision provides their meaning of \emph{small} and \emph{large}.

This function requires the user to provide either the \code{qua} or the \code{cdf} and \code{pdf} functions, which is somewhat divergent from the typical flow of logic of \pkg{lmomco}. This has been done so that \code{\link{expect.max.ostat}} can be used readily for experimental distribution functions. It is suggested that the parameter object be left in the \pkg{lmomco} style (see \code{\link{vec2par}}) even if the user is providing their own distribution functions.

Last comments: This function is built around the idea that either (1) the \code{cdf} and \code{pdf} ensemble or (2) \code{qua} exist in some clean analytical form and therefore the \code{qua=NULL} is the trigger on which order statistic expectation integral is used. This precludes an attempt to compute the support of the distribution  internally, and thus providing possibly superior (more refined) \code{lower} and \code{upper} limits. Here is a suggested re-implementation using the support of the Generalized Extreme Value distribution:
\preformatted{
para <- vec2par(c(100, 23, -0.5), type="gev")
lo <- quagev(0, para) # The value 54
hi <- quagev(1, para) # Infinity
E22 <- expect.max.ostat(2, para=para,cdf=cdfgev,pdf=pdfgev,
                           lower=lo, upper=hi)
E21 <- expect.min.ostat(2, para=para,cdf=cdfgev,pdf=pdfgev,
                           lower=lo, upper=hi)
L2 <- (E22 - E21)/2 # definition of L-scale
cat("L-scale: ",L2,"(integration)",
    lmomgev(para)$lambdas[2], "(theory)\n")
# The results show 33.77202 as L-scale.
}
The design intent makes it possible for some arbitrary and(or) new quantile function with difficult \code{cdf} and \code{pdf} expressions (or numerical approximations) to not be needed as the L-moments are explored.  Contrarily, perhaps some new \code{pdf} exists and simple integration of it is made to get the \code{cdf} but the \code{qua} would need more elaborate numerics to invert the \code{cdf}. The user could then still explore the L-moments with supervision on the integration limits or foreknowledge of the support of the distribution.
}
\value{
The expectation of the maximum order statistic, unless \eqn{j} is specified and then the expectation of that order statistic is returned. This similarly holds if the \code{\link{expect.min.ostat}} function is used except \dQuote{maximum} becomes the \dQuote{minimum}.

Alternatively, an \R \code{list} is returned.

\item{type}{The type of approach used: \dQuote{bypdfcdf} means the PDF and CDF of the distribution were used, and alternatively \dQuote{byqua} means that the quantile function was used.}
\item{value}{See previous discussion of value.}
\item{abs.error}{Estimate of the modulus of the absolute error from \R function \code{integrate}.}
\item{subdivisions}{The number of subintervals produced in the subdivision process from \R function \code{integrate}.}
\item{message}{\dQuote{OK} or a character string giving the error message.}
}
\note{
A function such as this might be helpful for computations involving distribution mixtures. Mixtures are readily made using the algebra of quantile functions (Gilchrist, 2000; Asquith, 2011, sec. 2.1.5 \dQuote{The Algebra of Quantile Functions}).

Last comments: Internally, judicious use of logarithms and exponents for the terms involving the \eqn{F} and \eqn{1-F} and the quantities to the left of the intergrals shown above are made in an attempt to maximize stability of the function without the user having to become too invested in the \code{lower} and \code{upper} limits. For example, \eqn{(1-F)^{n-j} \rightarrow \exp([n-j]\log(1-F))}. Testing indicates that this coding practice is quite useful. But there will undoubtedly be times for which the user needs to be informed enough about the expected value on return to identify that tweaking to the integration limits is needed. Also use of \R functions \code{lbeta} and \code{lchoose} is made to maximize operations in logarithmic space.

For \pkg{lmomco} v.2.1.+: Because of the extensive use of exponents and logarithms as described, enhanced deep tail estimation of the extrema for large \eqn{n} and large or small \eqn{j} results. This has come at the expense that expectations can be computed when the expectations actually do not exist. An error in the integration no longer occurs in \pkg{lmomco}. For example, the Cauchy distribution has infinite extrema but this function (for least for a selected parameter set and \code{n=10}) provides apparent values for the \eqn{\mathrm{E}[X_{1:n}]} and \eqn{\mathrm{E}[X_{n:n}]} when the \code{cdf} and \code{pdf} are used but not when the \code{qua} is used. Users are cautioned to not rely on \code{\link{expect.max.ostat}} \dQuote{knowing} that a given distribution has undefined order statistic extrema. Now for the Cauchy case just described, the extrema for \eqn{j = [1,n]} are hugely(!) greater in magnitude than for \eqn{j = [2,(n-1)]}, so some resemblance of \emph{infinity} remains.

The alias \code{\link{eostat}} is a shorter name dispatching to \code{\link{expect.max.ostat}} all of the arguments.
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Hosking, J.R.M., 2006, On the characterization of distributions by their L-moments: Journal of Statistical Planning and Inference, v. 136, no. 1, pp. 193--198.
}
\author{ W.H. Asquith}
\seealso{\code{\link{theoLmoms.max.ostat}}, \code{\link{expect.min.ostat}}, \code{\link{eostat}}}
\examples{
para <- vec2par(c(10,100), type="nor")
n <- 12
# The three outputted values from should be similar:
# (1) theoretical, (2) theoretical, and (3) simulation
expect.max.ostat(n, para=para, cdf=cdfnor, pdf=pdfnor)
expect.max.ostat(n, para=para, qua=quanor)
mean(sapply(1:1000, function(x) { max(rlmomco(n,para))}))

eostat(8, j=5, qua=quagum, para=vec2par(c(1670,1000), type="gum"))
\dontrun{
para <- vec2par(c(1280, 800), type="nor")
expect.max.ostat(10, j=9, para, qua=quanor)
[1] 2081.086      # SUCCESS ---------------------------
expect.max.ostat(10, j=9, para, pdf=pdfnor, cdf=cdfnor,
                                lower=-1E3, upper=1E6)
[1] 1.662701e-06  # FAILURE ---------------------------
expect.max.ostat(10, j=9, para, pdf=pdfnor, cdf=cdfnor,
                                lower=-1E3, upper=1E5)
[1] 2081.086      # SUCCESS ---------------------------
}
}
\keyword{order statistics (misc.)}

