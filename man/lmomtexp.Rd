\name{lmomtexp}
\alias{lmomtexp}
\title{L-moments of the Truncated Exponential Distribution}
\description{
This function estimates the L-moments of the Truncated Exponential distribution. The parameter \eqn{\psi} is the right truncation of the distribution and \eqn{\alpha} is a scale parameter, letting \eqn{\beta = 1/\alpha} to match nomenclature of Vogel and others (2008), the L-moments in terms of the parameters, letting \eqn{\eta = \mathrm{exp}(-\alpha\psi)}, are
\deqn{\lambda_1 = \frac{1}{\beta} - \frac{\psi\eta}{1-\eta} \mbox{,}}
\deqn{\lambda_2 = \frac{1}{1-\eta}\biggl[\frac{1+\eta}{2\beta} -
                                         \frac{\psi\eta}{1-\eta}\biggr] \mbox{,}}
\deqn{\lambda_3 = \frac{1}{(1-\eta)^2}\biggl[\frac{1+10\eta+\eta^2}{6\alpha} -
                                         \frac{\psi\eta(1+\eta)}{1-\eta}\biggr] \mbox{, and}}
\deqn{\lambda_4 = \frac{1}{(1-\eta)^3}\biggl[\frac{1+29\eta+29\eta^2+\eta^3}{12\alpha} -
                                         \frac{\psi\eta(1+3\eta+\eta^2)}{1-\eta}\biggr] \mbox{.}}

The distribution is restricted to a narrow range of L-CV (\eqn{\tau_2 = \lambda_2/\lambda_1}). If \eqn{\tau_2 = 1/3}, the process represented is a stationary Poisson for which the probability density function is simply the uniform distribution and \eqn{f(x) = 1/\psi}. If \eqn{\tau_2 = 1/2}, then the distribution is represented as the usual exponential distribution with a location parameter of zero and a scale parameter \eqn{1/\beta}. Both of these limiting conditions are supported.

If the distribution shows to be Uniform (\eqn{\tau_2 = 1/3}), then \eqn{\lambda_1 = \psi/2}, \eqn{\lambda_2 = \psi/6}, \eqn{\tau_3 = 0}, and \eqn{\tau_4 = 0}. If the distribution shows to be Exponential (\eqn{\tau_2 = 1/2}), then \eqn{\lambda_1 = \alpha}, \eqn{\lambda_2 = \alpha/2}, \eqn{\tau_3 = 1/3} and \eqn{\tau_4 = 1/6}.
}
\usage{
lmomtexp(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
\item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomtexp}.}
}
\references{
Vogel, R.M., Hosking, J.R.M., Elphick, C.S., Roberts, D.L., and Reed, J.M., 2008, Goodness of fit of probability distributions for sightings as species approach extinction: Bulletin of Mathematical Biology, DOI 10.1007/s11538-008-9377-3, 19 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{partexp}}, \code{\link{cdftexp}}, \code{\link{pdftexp}}, \code{\link{quatexp}} }
\examples{
set.seed(1) # to get a suitable L-CV
X <- rexp(1000, rate=.001) + 100
Y <- X[X <= 2000]
lmr <- lmoms(Y)

print(lmr$lambdas)
print(lmomtexp(partexp(lmr))$lambdas)

print(lmr$ratios)
print(lmomtexp(partexp(lmr))$ratios)
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Exponential (trimmed)}
