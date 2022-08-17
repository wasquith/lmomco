\encoding{utf8}
\name{gini.mean.diff}
\alias{gini.mean.diff}
\title{Gini Mean Difference Statistic }
\description{
The Gini mean difference statistic \eqn{\mathcal{G}} is a robust estimator of distribution scale and is closely related to the second L-moment \eqn{\lambda_2 = \mathcal{G}/2}.
\deqn{\mathcal{G} = \frac{2}{n(n-1)}\sum_{i=1}^n (2i - n - 1) x_{i:n}\mbox{,}}
where \eqn{x_{i:n}} are the sample order statistics.
}
\usage{
gini.mean.diff(x)
}
\arguments{
  \item{x}{A vector of data values that will be reduced to non-missing values.}
}
\value{
  An \R \code{list} is returned.
  \item{gini}{The gini mean difference \eqn{\mathcal{G}}.}
  \item{L2}{The L-scale (second L-moment) because \eqn{\lambda_2 = 0.5\times\mathcal{G}} (see \code{\link{lmom.ub}}).}
  \item{source}{An attribute identifying the computational source of the Gini's Mean Difference: \dQuote{gini.mean.diff}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

\enc{Jurečková}{Jureckova}, J., and Picek, J., 2006, Robust statistical methods with R: Boca Raton, Fla., Chapman and Hall/CRC, ISBN 1--58488--454--1.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmoms}}}
\examples{
fake.dat <- c(123, 34, 4, 654, 37, 78)
gini <- gini.mean.diff(fake.dat)
lmr <- lmoms(fake.dat)
str(gini)
print(abs(gini$L2 - lmr$lambdas[2]))
}
\keyword{special function}
