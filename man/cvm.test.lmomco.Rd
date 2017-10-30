\encoding{utf8}
\name{cvm.test.lmomco}
\alias{cvm.test.lmomco}
\title{ \enc{Cramér--von Mises}{Cramer--von Mises} Test for Goodness-of-Fit }
\description{
The \enc{Cramér--von Mises}{Cramer--von Mises} test for goodness-of-fit is implemented for the order statistics \eqn{x_{1:n} \le x_{i:n} \le x_{n:n}} of a sample of size \eqn{n}. Define the test statistic (\enc{Csörgő}{Csorgo} and Faraway, 1996) as
\deqn{\omega^2 = \frac{1}{12n} + \sum_{i=1}^n \biggl[\frac{2i-1}{2n} - F_\theta(x_i)\biggr]\mbox{,}}
where \eqn{F_\theta(x)} is the cumulative distribution function (continuous) for some distribution having parameters \eqn{\theta}. If the value for \eqn{\omega^2} is larger than some critical value, reject the null hypothesis. The null hypothesis is that \eqn{F} is the function specified by \eqn{\theta}, while the alternative hypothesis is that \eqn{F} is some other function.
}
\usage{
cvm.test.lmomco(x, para1, ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{para1}{The parameters of the distribution.}
  \item{...}{Additional arguments to pass to \code{\link{par2cdf}}.}
}
\value{
  An \R \code{list} is returned.
  \item{null.dist}{The null distribution, which is an echoing of the \code{para} argument, which recall for \pkg{lmomco} that is contains the distribution abbreviation.}
  \item{text}{The string \dQuote{Cramer-von Mises test of goodness-of-fit}.}
  \item{statistic}{The \eqn{\omega^2} as defined above (see \bold{Note}).}
  \item{p.value}{The p-value computed from the \code{pCvM()} function from the \pkg{goftest} package for the null distribution of the test statistic.}
  \item{source}{An attribute identifying the computational source of the L-moments:\cr \dQuote{cvm.test.lmomco}.}
}
\details{
The above definition for \eqn{\omega^2} as the \enc{Cramér--von Mises}{Cramer--von Mises} test statistic is consistent with the notation in \enc{Csörgő}{Csorgo} and Faraway (1996) as well as that in package \pkg{goftest}. Depending on how the null distribution is defined by other authors and attendant notation, the \enc{Cramér--von Mises}{Cramer--von Mises} statistic can be branded as \eqn{T = n\omega^2}. The null distribution herein requires just \eqn{\omega^2} and the sample size is delivered separately into the cumulative distribution function:
\preformatted{
  goftest::pCvM(omega.sq, n=n, lower.tail=FALSE)
}
}
\note{
An example of coverage probabilities demonstrating the differences in what the p-values mean on whether the parent is known or the "parent" is coming from the sample. The p-values are quite different and inference has subtle differences. In ensemble, comparing the test statistic amongst distribution choices might be more informative than a focus on p-values being below a critical alpha.
\preformatted{
  parent <- vec2par(c(20,120), type="gam"); nsim <- 10000
  pp <- nn <- ee <- rep(NA,nsim)
  for(i in 1:nsim) {
    x <- rlmomco(56, parent); lmr <- lmoms(x)
    pp[i] <- cvm.test.lmomco(x,          parent          )$p.value
    nn[i] <- cvm.test.lmomco(x, lmom2par(lmr, type="nor"))$p.value
    ee[i] <- cvm.test.lmomco(x, lmom2par(lmr, type="exp"))$p.value
  }
  message("GAMMA PARENT KNOWN     'rejection rate'=", sum(pp < 0.05)/nsim)
  message("ESTIMATED NORMAL       'rejection rate'=", sum(nn < 0.05)/nsim)
  message("ESTIMATED EXPONENTIAL  'rejection rate'=", sum(ee < 0.05)/nsim)
}
The rejection rate for the Gamma is about 5 percent, which matches the 0.05 specified in the conditional. The Normal is about zero, and the Exponential is about 21 percent. The fitted Normal almost always passes for the real parent, though Gamma, for the sample size and amount of L-skewness involved. The Exponential does not. This illustrates that the p-value can be misleading in the single-sample version of this test. Thus when fit by parameters from the sample, the test statistic is nearly always smaller than the one for a prespecified set of parameters. The significance level will be smaller than intended.
}
\references{
\enc{Csörgő}{Csorgo}, S., and Faraway, J.J., 1996, The exact and asymptotic distributions of \enc{Cramér--von Mises}{Cramer--von Mises} statistics: Journal of the Royal Statistical Society, Series B, v. 58, pp. 221--234.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lmrdia}} }
\examples{
# An example in which the test is conducted on a sample but the parent is known.
# This will lead to more precise inference than if the sample parameters are used.
mu <- 120; sd <- 25; para <- vec2par(c(120,25), type="nor")
x <- rnorm(56, mean=mu, sd=sd)
T1 <- cvm.test.lmomco(x, para)$statistic
T2 <- goftest::cvm.test(x, null="pnorm", mean=mu, sd=sd)$statistic
message("Cramer--von Mises: T1=", round(T1, digits=6), " and T2=", round(T2, digits=6))
}
\keyword{goodness-of-fit}
\keyword{Cramer--von Mises statistic}
