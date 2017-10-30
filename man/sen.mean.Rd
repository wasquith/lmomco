\encoding{UTF-8}
\name{sen.mean}
\alias{sen.mean}
\title{ Sen's Weighted Mean Statistic }
\description{
The Sen weighted mean statistic \eqn{\mathcal{S}_{n,k}} is a robust estimator of the mean of a distribution
\deqn{\mathcal{S}_{n,k} = {n \choose 2k+1}^{-1} \sum_{i=1}^n {i - 1 \choose k} {n - i \choose k } x_{i:n}\mbox{,}}
where \eqn{x_{i:n}} are the sample order statistics and \eqn{k} is a weighting or trimming parameter. If \eqn{k = 2}, then the \eqn{\mathcal{S}_{n,2}} is the first symmetrical TL-moment (trim = 1).

Note that \eqn{\mathcal{S}_{n,0} = \mu = \overline{X}_n} or the arithmetic \code{mean} and \eqn{\mathcal{S}_{n,k}} is the sample \code{median} if either \eqn{n} is even and \eqn{k = (n/2) - 1} or \eqn{n} is odd and \eqn{k = (n-1)/2}.
}
\usage{
sen.mean(x,k=0)
}
\arguments{
  \item{x}{A vector of data values that will be reduced to non-missing values.}
  \item{k}{A weighting or trimming parameter \eqn{0 < k < (n-1)/2}.}
}
\value{
  An \R \code{list} is returned.

  \item{sen}{The sen mean \eqn{\mathcal{S}_{n,k}}}
  \item{source}{An attribute identifying the computational source: \dQuote{sen.mean}.}
}
\references{
\enc{Jurečková}{Jureckova}, J., and Picek, J., 2006, Robust statistical methods with R: Boca Raton, Fla., Chapman and Hall/CRC, ISBN 1--58488--454--1, 197 p.

Sen, P.K., 1964, On some properties of the rank-weighted means: Journal Indian Society of Agricultural Statistics: v. 16, pp. 51--61.
}
\author{ W.H. Asquith}
\seealso{\code{\link{TLmoms}}, \code{\link{gini.mean.diff}}}
\examples{
fake.dat <- c(123,34,4,654,37,78)
sen.mean(fake.dat); mean(fake.dat) # These should be the same values

sen.mean(fake.dat, k=(length(fake.dat)/2) - 1); median(fake.dat)
# Again, same values

# Finally, the sen.mean() is like a symmetrically trimmed TL-moment
# Let us demonstrate by computed a two sample trimming for each side
# for a Normal distribution having a mean of 100.
fake.dat <- rnorm(20, mean=100)
lmr <- TLmoms(fake.dat, trim=2)
sen <- sen.mean(fake.dat, k=2)

print(abs(lmr$lambdas[1] - sen$sen)) # zero is returned
}
\keyword{special function}

