\name{qua.ostat}
\alias{qua.ostat}
\title{Compute the Quantiles of the Distribution of an Order Statistic}
\description{
This function computes a specified quantile by nonexceedance probability \eqn{F} for the \eqn{j}th-order statistic of a sample of size \eqn{n} for a given distribution. Let the quantile function (inverse distribution) of the Beta distribution be

\deqn{ \mathrm{B}^{(-1)}(F,j,n-j+1) \mbox{,} }

and let \eqn{x(F,\Theta)} represent the quantile function of the given distribution and \eqn{\Theta} represents a vector of distribution parameters. The quantile function of the distribution of the \eqn{j}th-order statistic is

\deqn{ x\bigl(\mathrm{B}^{(-1)}(F,j,n-j+1),\Theta\bigr) \mbox{.} }
}
\usage{
qua.ostat(f, j, n, para=NULL)
}
\arguments{
  \item{f}{The nonexceedance probability \eqn{F} for the quantile.}
  \item{j}{The \eqn{j}th-order statistic \eqn{x_{1:n} \le x_{2:n} \le \ldots \le x_{j:n} \le x_{n:n}.}}
  \item{n}{The sample size.}
  \item{para}{A distribution parameter list from a function such as \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  The quantile of the distribution of the \eqn{j}th-order statistic is returned.
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton, Fla.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmom2par}}, \code{\link{vec2par}}}
\examples{
gpa <- vec2par(c(100, 500, 0.5), type="gpa")
n <- 20   # the sample size
j <- 15   # the 15th order statistic
F <- 0.99 # the 99th percentile
theoOstat <- qua.ostat(F, j, n, gpa)

\dontrun{
# Let us test this value against a brute force estimate.
Jth <- vector(mode="numeric")
for(i in seq_len(50000)) {
  Q <- sort( rlmomco(n, gpa) )
  Jth[i] <- Q[j]
}
bruteOstat <- quantile(Jth, F) # estimate by built-in function
theoOstat  <- signif( theoOstat, digits=5)
bruteOstat <- signif(bruteOstat, digits=5)
cat(c("Theoretical=", theoOstat, "  Simulated=", bruteOstat, "\n")) # }
}
\keyword{order statistics (misc.)}
\keyword{order statistics (quantile function of)}
