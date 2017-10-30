\name{tau34sq.normtest}
\alias{tau34sq.normtest}
\title{The Tau34-squared Test: A Normality Test based on L-skew and L-kurtosis and an Elliptical Rejection Region on an L-moment Ratio Diagram}
\description{
This function performs highly intriguing test for normality using L-skew (\eqn{\tau_3}) and L-kurtosis (\eqn{\tau_4}) computed from an input vector of data. The test is simultaneously focused on L-skew and L-kurtosis. Harri and Coble (2011) presented two types of normality tests based on these two L-moment ratios. Their first test is dubbed the \eqn{\tau_3\tau_4} test. Those authors however conclude that a second test dubbed the \eqn{\tau^2_{3,4}} test \dQuote{in particular shows consistently high power against [sic] symmetric distributions and also against [sic] skewed distributions and is a powerful test that can be applied against a variety of distributions.}

A sample-size transformed quantity of the sample L-skew (\eqn{\hat\tau_3}) is
\deqn{Z(\tau_3) = \hat\tau_3 \times \frac{1}{\sqrt{0.1866/n + 0.8/n^2}}\mathrm{,}}
which has an approximate Standard Normal distribution. A sample-sized transformation of the sample L-kurtosis (\eqn{\hat\tau_4}) is
\deqn{Z(\tau_4)' = \hat\tau_4 \times \frac{1}{\sqrt{0.0883/n}}\mathrm{,}}
which also has an approximate Standard Normal distribution. A superior approximation for the variate of the Standard Normal distribution however is
\deqn{Z(\tau_4) = \hat\tau_4 \times \frac{1}{\sqrt{0.0883/n + 0.68/n^2 + 4.9/n^3}}\mathrm{,}}
and is highly preferred for the algorithms in \code{\link{tau34sq.normtest}}.


The \eqn{\tau_3\tau_4} test (not implemented in \code{\link{tau34sq.normtest}}) by Harri and Coble (2011) can be constructed from the \eqn{Z(\tau_3)} and \eqn{Z(\tau_4)} statistics as shown, and a square rejection region constructed on an L-moment ratio diagram of L-skew versus L-kurtosis. However, the preferred method is the \dQuote{Tau34-squared} test \eqn{\tau^2_{3,4}} that can be developed by expressing an ellipse on the L-moment ratio diagram of L-skew versus L-kurtosis. The \eqn{\tau^2_{3,4}} test statistic is defined as
\deqn{\tau^2_{3,4} = Z(\tau_3)^2 + Z(\tau_4)^2\mathrm{,}}
which is approximately distributed as a \eqn{\chi^2} distribution with two degrees of freedom. The \eqn{\tau^2_{3,4}} also is the expression of the ellipical region on the L-moment ratio diagram of L-skew versus L-kurtosis.
}
\usage{
tau34sq.normtest(x, alpha=0.05, pvalue.only=FALSE, getlist=TRUE,
                    useHoskingZt4=TRUE, verbose=FALSE, digits=4)
}
\arguments{
  \item{x}{A vector of values.}
  \item{alpha}{The \eqn{\alpha} significance level.}
  \item{pvalue.only}{Only return the p-value of the test and superceeds the \code{getlist} argument.}
  \item{getlist}{Return a list of salient parts of the computations.}
  \item{useHoskingZt4}{J.R.M. Hosking provided a better approximation \eqn{Z(\tau_4)} in personal correspondance to Harri and Coble (2011) than the one \eqn{Z(\tau_4)'} they first presented in their paper. This argument is a logical on whether this approximation should be used. It is highly recommended that \code{useHoskingZt4} be left at the default setting.}
  \item{verbose}{Print a nice summary of the test.}
  \item{digits}{How many digits to report in the summary.}
}
\value{
  An \R \code{list} is returned if \code{getlist} argument is true. The list contents are

  \item{SampleTau3}{The sample L-skew.}
  \item{SampleTau4}{The sample L-kurtosis.}
  \item{Ztau3}{The Z-value of \eqn{\tau_3}.}
  \item{Ztau4}{The Z-value of \eqn{\tau_4}.}
  \item{Tau34sq}{The \eqn{\tau^2_{3,4}} value.}
  \item{ChiSq.2df}{The Chi-squared distribution nonexceedance probability.}
  \item{pvalue}{The p-value of the test.}
  \item{isSig}{A logical on whether the p-value is \dQuote{statistically significant} based on the \eqn{\alpha} value.}
  \item{source}{The source of the parameters: \dQuote{tau34sq.normtest}.}
}
\references{
Harri, A., and Coble, K.H., 2011, Normality testing---Two new tests using L-moments: Journal of Applied Statistics, v. 38, no. 7, pp. 1369--1379.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfnor}}, \code{\link{plotlmrdia}} }
\examples{
HarriCoble <- tau34sq.normtest(rnorm(20), verbose=TRUE)
\dontrun{
# If this basic algorithm is run repeatedly with different arguments,
# then the first three rows of table 1 in Harri and Coble (2011) can
# basically be repeated. Testing by WHA indicates that even better
# empirical alphas will be computed compared to those reported in that table 1.
# R --vanilla --silent --args n 20 s 100 < t34.R
# Below is file t34.R
library(batch) # for command line argument parsing
a <- 0.05; n <- 50; s <- 5E5 # defaults
parseCommandArgs() # it will echo out those arguments on command line
sims <- sapply(1:s, function(i) {
          return(tau34sq.normtest(rnorm(n),
                 pvalue.only=TRUE)) })
p <- length(sims[sims <= a])
print("RESULTS(Alpha, SampleSize, EmpiricalAlpha)")
print(c(a, n, p/s))
}
}
\keyword{L-moment (hypothesis)}

