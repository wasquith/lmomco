\name{lmoms.cov}
\alias{lmoms.cov}
\title{Distribution-Free Variance-Covariance Structure of Sample L-moments }
\description{
Compute the distribution-free, variance-covariance matrix (\eqn{\widehat{\mathrm{var}}(\lambda)}) of the sample L-moments (\eqn{\hat\lambda_r}) or alternatively the sample probability-weighted moments (\eqn{\hat\beta_k}, Elamir and Seheult, 2004, sec. 5). The \eqn{\widehat{\mathrm{var}}(\lambda)} is defined by the matrix product

\deqn{\widehat{\mathrm{var}}(\lambda) = \mathbf{C}\,\mathbf{\hat\Theta}\,\mathbf{C}^{\mathrm{T}}\mbox{,}}

where the \eqn{r \times r} matrix \eqn{\mathbf{C}} for number of moments \eqn{r} represents the coefficients of the linear combinations converting \eqn{\beta_k} to \eqn{\lambda_r} and the \eqn{r}th row in the matrix is defined as

\deqn{\mathbf{C}[r,]_{k{=}0:(r-1)} = (-1)^{(r-1-k)} {r-1 \choose k} {r-1+k \choose k}\mbox{,}}

where the row is padded from the right with zeros for \eqn{k < r} to form the required lower triangular structure.  Elamir and Seheult (2004) list the \eqn{\mathbf{C}} matrix for \eqn{r = 4}.

Letting the \emph{falling factorial} be defined (matching Elamir and Seheult's nomenclature) as
\deqn{a^{(b)} = \Gamma(b+1) {a \choose b}\mbox{,}}
and letting an entry in the \eqn{\mathbf{\hat\Theta}} matrix denoted as \eqn{\hat\theta_{kl}} be defined as

\deqn{\hat\theta_{kl} = \hat\beta_k\hat\beta_l - \frac{A}{n^{(k+l+2)}}\mbox{,}}

where \eqn{\hat\beta_k} are again the sample probability-weighted moments and are computed by \code{\link{pwm}}, and finally \eqn{A} is defined as
\deqn{A = \sum_{i=1}^{n-1}\sum_{j=i+1}^{n}
  \bigl[ (i-1)^{(k)} (j-k-2)^{(l)} +
          (i-1)^{(l)} (i-l-2)^{(k)}
  \bigr] X_{i:n}X_{j:n}\mbox{,}}
where \eqn{X_{i:n}} are the sample order statistics for a sample of size \eqn{n}.

Incidentally, the matrix \eqn{\mathbf{\hat\Theta}} is the variance-covariance structure (\eqn{\widehat{\mathrm{var}}}) of the \eqn{\hat\beta}, thus \eqn{\widehat{\mathrm{var}}(\beta) = \mathbf{\hat\Theta}}, which can be returned by a logical function argument (\code{as.pwm=TRUE}) instead of \eqn{\widehat{\mathrm{var}}(\lambda)}. The last example in \bold{Examples} provides a demonstration.
}
\usage{
lmoms.cov(x, nmom=5, as.pwm=FALSE, showC=FALSE,
             se=c("NA", "lamse", "lmrse", "pwmse"), ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{as.pwm}{A logical controlling whether the distribution-free, variance-covariance of sample probability-weighted moments (\eqn{\mathbf{\hat\Theta}}) is returned instead.}
  \item{showC}{A logical controlling whether the matrix \eqn{\mathbf{C}} is printed during function operation, and this matrix is not returned as a presumed safety feature.}
  \item{se}{Compute standard errors (\eqn{SE}) for the respective moments. The default of \code{"NA"} retains the return of either \eqn{\widehat{\mathrm{var}}(\beta)} or \eqn{\widehat{\mathrm{var}}(\lambda)} depending on setting of \code{as.pwm}. The \code{"lamse"} returns the square root of the diagonal of \eqn{\widehat{\mathrm{var}}(\lambda)}, and notationally these are \eqn{\lambda_r^{SE}}. Similarly, \code{"pwmse"} returns the square root of the diagonal of \eqn{\widehat{\mathrm{var}}(\beta)} by internally setting \code{as.pwm} to \code{TRUE}, and notationally these are \eqn{\beta_{r-1}^{SE}}. (Remember that \eqn{\beta_0 \equiv \lambda_1}---the indexing of the former starts at 0 and at the later at 1). The \code{"lmrse"} returns the square root of the first two terms of the \eqn{\widehat{\mathrm{var}}(\lambda)} diagonal (\eqn{\lambda_{1,2}^{SE}}) but computes \eqn{SE} for the L-moment ratios (\eqn{\tau_r^{SE}}) for \eqn{r \ge 3} using the Taylor-series-based approximation (see \bold{Note}) shown by Elamir and Seheult (2004, p. 348). (Remember that L-moment ratios are \eqn{\tau_r = \lambda_r/\lambda_2} for \eqn{r \ge 3} and that \eqn{\tau_2 = \lambda_2/\lambda_1} [coefficient of L-variation].)}
  \item{...}{Other arguments to pass should they be needed (none were at first implementation).}
}
\value{
  An \R \code{matrix} is returned. In small samples and substantially sized \eqn{r}, one or more \eqn{\hat\theta_{kl}} will be \code{NaN} starting from the lower right corner of the matrix.  The function does not test for this nor reduce the number of moments declared in \code{nmom} itself. To reiterate, the square roots along the \eqn{\widehat{\mathrm{var}}(\lambda)} diagonal are \eqn{SE} for the respective L-moments.
}
\note{
Function \code{lmoms.cov} was developed as a double check on the evidently separately developed \eqn{r \le 4} (\code{nmom}) implementations of \eqn{\widehat{\mathrm{var}}(\lambda)} in packages \pkg{Lmoments} and \pkg{nsRFA}. Also the internal structure closely matches the symbolic mathematics by Elamir and Seheult (2004), but this practice comes at the expense of more than an order of magnitude slower execution times than say either of the functions \code{Lmomcov()} (package \pkg{Lmoments}) or \code{varLmoments()} (package \pkg{nsRFA}). For a high speed and recommended implementation, please use the \pkg{Lmoments} package by Karvanen (2016)---Karvanen extended this implementation to larger \eqn{r} for the \pkg{lmomco} package.

For \code{se="lmrse"}, the Taylor-series-based approximation is suggested by Elamir and Seheult (2004, p. 348) to estimate the variance of an L-moment ratio (\eqn{\tau_r} for \eqn{r \ge 3}) is based on structure of the variance of the ratio of two uniform variables in which the numerator is the \eqn{r}th L-moment and the denominator is \eqn{\lambda_2}:
\deqn{\mathrm{var}(\tau_r) \cong
\biggl[ \frac{\mathrm{var}(\lambda_r)}{\mathrm{E}(\lambda_r)^2} + \frac{\mathrm{var}(\lambda_2)}{\mathrm{E}(\lambda_2)^2} - \frac{2\mathrm{cov}(\lambda_r,\lambda_2)}{\mathrm{E}(\lambda_r)\mathrm{E}(\lambda_2)} \biggr]
\biggl[\frac{\mathrm{E}(\lambda_r)}{\mathrm{E}(\lambda_2)} \biggr]^2\mbox{,}}
where \eqn{\mathrm{var}(\cdots)} are the along the diagonal of \eqn{\widehat{\mathrm{var}}(\lambda)} and \eqn{\mathrm{cov}(\cdots)} are the off-diagonal covariances. The expectations \eqn{\mathrm{E}(\cdots)} are replaced with the sample estimates. Only for \code{se="lmrse"} the \eqn{SE} of the coefficient of L-variation (\eqn{\tau_2^{SE}}) is computed but retained as an attribute (\code{attr()} function) of the returned vector and not housed within the vector---the \eqn{\lambda_2^{SE}} continues to be held in the 2nd position of the returned vector.
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2004, Exact variance structure of sample L-moments: Journal of Statistical Planning and Inference, v. 124, pp. 337--359.

Karvanen, Juha, 2016, Lmoments---L-moments and quantile mixtures: R package version 1.2-3, accessed February 22, 2016 at https://cran.r-project.org/web/packages/Lmoments/index.html
}
\author{W.H. Asquith}
\seealso{\code{\link{lmoms}}, \code{\link{pwm}}
}
\examples{
\dontrun{
nsim <- 1000; n <- 10 # Let us compute variance of lambda_3
VL3sample <- mean(replicate(nsim, { zz <- lmoms.cov(rexp(n),nmom=3); zz[3,3] }))
falling.factorial <- function(a, b) gamma(b+1)*choose(a,b)
VL3exact  <- ((4*n^2 - 3*n - 2)/30)/falling.factorial (10, 3) # Exact variance is from
print(c(VL3sample, VL3exact)) # Elamir and Seheult (2004, table 1, line 8)
#[1] 0.01755058 0.01703704  # the values obviously are consistent}
\dontrun{
# Data considered by Elamir and Seheult (2004, p. 348)
library(MASS); data(michelson); Light <- michelson$Speed
lmoms(Light, nmom=4)$lambdas # 852.4, 44.3, 0.83, 6.5 # matches those authors
lmoms.cov(Light) # [1, ] ==> 62.4267, 0.7116, 2.5912, -3.9847 # again matches
# The authors report standard error of L-kurtosis as 0.03695, which matches
lmoms.cov(Light, se="lmrse")[4] # 0.03695004 }
\dontrun{
D <- rnorm(100) # Check results of Lmoments package.
lmoms.cov(D, rmax=5)[,5]
#        lam1         lam2         lam3         lam4         lam5
#3.662721e-04 3.118812e-05 5.769509e-05 6.574662e-05 1.603578e-04
Lmoments::Lmomcov(D, rmax=5)[,5]
#          L1           L2           L3           L4           L5
#3.662721e-04 3.118812e-05 5.769509e-05 6.574662e-05 1.603578e-04}
}
\keyword{L-moment (sample)}
\keyword{PWM (sample)}
\keyword{L-moment (variance-covariance)}
\keyword{variance-covariance}

