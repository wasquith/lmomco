\encoding{utf8}
\name{lmrloc}
\alias{lmrloc}
\title{Line-of-Organic Correlation}
\description{
Compute the line-of-organic correlation (LOC) (Helsel and others, 2020, sec. 10.2.2, p. 280). The LOC is estimated by both L-moments and product moments. The LOC has other names in the literature including reduced major axis and line of diagonal correlation. When describing a functional relations between two variables without trying to predict one from the other, LOC is more appropriate than ordinary least squares (OLS).

The LOC is a regression line whose slope is computed by the ratio between respective variations of the predictor variable and the response variable. The intercept of the line is computed such that the line passes through the familiar arithmetic mean (first L-moment) (\eqn{\lambda_1}) each for the two variables. Relative variation is readily computed by the ratio of standard deviations or for more robust and less biased estimation by the ratio of the L-variations (second L-moment) (\eqn{\lambda_2}) of the two variables.

The \eqn{\lambda_2} is generically based on the so-called Gini mean difference statistic (GMD) (\eqn{\mathcal{G}}) by \eqn{\lambda_2 = \mathcal{G}/2} (\code{\link{gini.mean.diff}}). Incidentally for the normal distribution, the well-known standard deviation is the product \eqn{\lambda_2\sqrt{\pi}} (see also \code{\link{lmomnor}}). Mathematically, GMD is defined as the linear combination
\deqn{\mathcal{G} = \frac{2}{n(n-1)}\sum_{i=1}^n (2i - n - 1) x_{i:n}\mbox{,}}
where \eqn{x_{i:n}} are the sample ascending order statistics.

Returning to the need to estimate the LOC slope, algebra shows the slope is the ratio of the \eqn{\mathcal{G}} values as
\deqn{m = \mathrm{sign[} \rho \mathrm{]}\cdot\frac{\sum_{i=1}^n (2i - n - 1) X_{i:n}}{\sum_{i=1}^n (2i - n - 1) Y_{i:n}}\mbox{,}}
where \eqn{X_{i:n}} is an ordered (ascending) vector of random variable \eqn{X}, \eqn{Y_{i:n}} is an ordered (ascending) vector of random variable \eqn{Y}, and the slope sign can be computed by a correlation coefficient sign (Pearson R, Kendall Tau [computationally slowest], Spearman Rho would all work [implemented for the function, \eqn{\rho}]). For applications, it is critical that the correlation coefficient is computed using the original correlated-ordering of \eqn{X} and \eqn{Y} and not after individual vector sorting that is needed for the GMD (L-moments). A developer, therefore, must be cognizant of the placement in code when the two variables are sorted to the order statistics for \eqn{\mathcal{G}} computations.

The LOC intercept is given by algebra by
\deqn{b = \frac{1}{n}\biggl(\sum_{i=1}^n X_{i:n} - m \cdot \sum_{i=1}^n Y_{i:n}\biggr)\mbox{.}}

Helsel and others (2020, p. 281) enumerate some advantages to the use of the LOC: (1) it minimizes errors in both \code{x} and \code{y} directions, (2) it provides a single line regardless of which variable (x or y) is used as the response variable, and (3) its cumulative distribution function of the predictions, including the variance and probabilities, is correct (meaning not compressed as in OLS). The LOC is particularly useful for modeling the intrinsic functional relation between two variables, both of which are measured with error and (or) when neither variable is considered an independent variable appropriate to predict the other.
}
\usage{
lmrloc(x, y=NULL, terse=TRUE)
}
\arguments{
  \item{x}{A numeric vector, matrix or data frame.}
  \item{y}{\code{NULL} (default) or a vector of same length of \code{x}.}
  \item{terse}{A logical triggering only return of the coefficients of the two lines; otherwise, the intermediate computations are also returned.}
}
\value{
An \R \code{list} is returned with \code{terse=TRUE} with two vectors of the intercept and slope coefficients for the L-moment and the product moment versions. The names on the vectors, respectively, are \code{"LMR_Intercept", "LMR_Slope"} and \code{"PMR_Intercept", "PMR_Slope"} for LMR (L-moment ratio) and PMR (product moment ratio) are monikers for the two approaches. An expanded \R \code{list} is returned with \code{terse=FALSE} with the intermediate computations also provided.
  \item{loc_lmr}{The LOC by L-moments (L-variations or equivalently Gini Mean Differences).}
  \item{loc_pmr}{The LOC by product moments (standard deviations).}
  \item{srho}{The sign on Spearman Rho.}
  \item{mu_x}{The arithmetic mean of the \code{x} variable.}
  \item{mu_y}{The arithmetic mean of the \code{y} variable.}
  \item{gini_x}{The GMD of the \code{x} variable.}
  \item{gini_y}{The GMD of the \code{y} variable.}
  \item{sd_x}{The standard deviation of the \code{x} variable.}
  \item{sd_y}{The standard deviation of the \code{y} variable.}
}
\references{
Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A., and Gilroy, E.J., 2020, Statistical methods in water resources: U.S. Geological Survey Techniques and Methods, book 4, chap. A3, 458 p., \doi{10.3133/tm4a3}.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

\enc{Jurečková}{Jureckova}, J., and Picek, J., 2006, Robust statistical methods with R: Boca Raton, Fla., Chapman and Hall/CRC, ISBN 1--58488--454--1.
}
\author{W.H. Asquith}
\seealso{\code{\link{gini.mean.diff}}}
\examples{
n <- 100; x <- rnorm(n); y <- -0.4 * x + rnorm(n, sd=0.2)
y[x == min(x)] <- 2 * min(y) # throw in an outlier to help separate two lines
loc <- lmrloc(x, y, terse=FALSE)
plot(x, y)
abline(loc$loc_lmr, lty=1)
abline(loc$loc_pmr, lty=2)
legend("topright", c("LOC by L-moments", "LOC by product moments"), lty=c(1,2))

olsxy <- 1 / coefficients(stats::lm(x~y))[2] # yes inversion needed to show
olsyx <-     coefficients(stats::lm(y~x))[2] # geometric mean in proper way
mstar <- loc$srho * sqrt(abs(olsxy) * abs(olsyx)); names(mstar) <- NULL
m_pmr <- loc$loc_pmr[2]; names(m_pmr) <- NULL
m_lmr <- loc$loc_lmr[2]; names(m_lmr) <- NULL
message("Geometric mean OLS slopes = ", mstar) # see that these two are
message("           PMR LOC slope  = ", m_pmr) # equivalent by theory
message("           LMR LOC slope  = ", m_lmr) # this one is not
}
\keyword{line-of-organic correlation}
\keyword{reduced-major-axis line}
\concept{line of organic correlation}
\concept{organic correlation}
\concept{reduced major axis}
\keyword{ordinary least squares}
