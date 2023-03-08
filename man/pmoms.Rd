\name{pmoms}
\alias{pmoms}
\title{The Sample Product Moments: Mean, Standard Deviation, Skew, and Excess Kurtosis}
\description{
Compute the first four sample product moments. Both classical (theoretical and biased) versions and unbiased (nearly) versions are produced. Readers are directed to the References and the source code for implementation details.
}
\usage{
pmoms(x)
}
\arguments{
  \item{x}{A real value vector.}
}
\value{
  An \R \code{list} is returned.

  \item{moments}{Vector of the product moments: first element is the mean (\code{mean} in \R), second is standard deviation, and the higher values typically are not used as these are not unbiased moments, but the ratios[3] and ratios[4] are nearly unbiased.}
  \item{ratios}{Vector of the product moment ratios. Second element is the coefficient of variation, ratios[3] is skew, and ratios[4] is kurtosis.}
  \item{sd}{Nearly unbiased standard deviation [well at least unbiased variance \cr (\code{unbiased.sd^2})] computed by \R function \code{sd}. }
  \item{umvu.sd}{Uniformly-minimum variance unbiased estimator of standard deviation.}
  \item{skew}{Nearly unbiased skew, same as ratios[3].}
  \item{kurt}{Nearly unbiased kurtosis, same as ratios[4].}
  \item{excesskurt}{Excess kurtosis from the Normal distribution: \code{kurt - 3}.}
  \item{classic.sd}{Classical (theoretical) definition of standard deviation.}
  \item{classic.skew}{Classical (theoretical) definition of skew.}
  \item{classic.kurt}{Classical (theoretical) definition of kurtosis}
  \item{classic.excesskurt}{Excess classical (theoretical) kurtosis from \cr Normal distribution: \code{classic.kurt - 3}.}
  \item{message}{The product moments are confusing in terms of definition because they are not naturally unbiased.  This characteristic is different from the L-moments. The author thinks that it is informative to show the biased versions within the \dQuote{classic} designations. Therefore, this \code{message} includes several clarifications of the output.}
  \item{source}{An attribute identifying the computational source (the function name) of the product moments: \dQuote{pmoms}.}
}
\note{
This function is primarily available for gamesmanship with the Pearson Type III distribution as its parameterization in \pkg{lmomco} returns the product moments as the very parameters of that distribution. This of course is like the Normal distribution in which the first two parameters are the first two product moments; the Pearson Type III just adds skew.  See the example below. Another reason for having this function in \pkg{lmomco} is that it demonstrates application of unbiased product moments and permits comparisons to the L-moments (see Asquith, 2011; figs. 12.13--12.16).

The \code{umvu.sd} is computed by
\deqn{\hat\sigma' = \frac{\Gamma[(n-1)/2]}{\Gamma(n/2)\sqrt{2}}\sqrt{\sum_{i=1}^{n} (x_i - \hat\mu)^2}\mbox{,}}

where \eqn{\hat\sigma'} is the estimate of standard deviation for the sample \eqn{x} of size \eqn{n}, \eqn{\Gamma(\cdots)} is the complete gamma function, and \eqn{\hat\mu} is the arthimetic mean.
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.

Joanes, D.N., Gill, C.A., 1998, Comparing measures of sample skewness and kurtosis: The Statistician, v. 47, no. 1, pp. 183--189.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmoms}}}
\examples{
# A simple example
PM <- pmoms(rnorm(1000)) # n standard normal values as a fake data set.
cat(c(PM$moments[1],PM$moments[2],PM$ratios[3],PM$ratios[4],"\n"))
# As sample size gets very large the four values returned should be
# 0,1,0,0 by definition of the standard normal distribution.

# A more complex example
para <- vec2par(c(100,500,3),type='pe3') # mean=100, sd=500, skew=3
# The Pearson type III distribution is implemented here such that
# the "parameters" are equal to the mean, standard deviation, and skew.
simDATA <- rlmomco(100,para) # simulate 100 observations
PM <- pmoms(simDATA) # compute the product moments

p.tmp <- c(PM$moments[1],PM$moments[2],PM$ratios[3])
cat(c("Sample P-moments:",p.tmp,"\n"))
# This distribution has considerable variation and large skew. Stability
# of the sample product moments requires LARGE sample sizes (too large
# for a builtin example)

# Continue the example through the L-moments
lmr <- lmoms(simDATA) # compute the L-moments
epara <- parpe3(lmr) # estimate the Pearson III parameters. This is a
# hack to back into comparative estimates of the product moments. This
# can only be done because we know that the parent distribution is a
# Pearson Type III

l.tmp <- c(epara$para[1],epara$para[2],epara$para[3])
cat(c("PearsonIII by L-moments:",l.tmp,"\n"))
# The first values are the means and will be identical and close to 100.
# The second values are the standard deviations and the L-moment to
#   PearsonIII will be closer to 500 than the product moment (this
#   shows the raw power of L-moment based analysis---they work).
# The third values are the skew. Almost certainly the L-moment estimate
#   of skew will be closer to 3 than the product moment.
}
\keyword{product moment}

