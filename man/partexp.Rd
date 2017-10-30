\name{partexp}
\alias{partexp}
\title{Estimate the Parameters of the Truncated Exponential Distribution}
\description{
This function estimates the parameters of the Truncated Exponential distribution given
the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The parameter \eqn{\psi} is the right truncation of the distribution, and \eqn{\alpha} is a scale parameter, letting \eqn{\beta = 1/\alpha} to match nomenclature of Vogel and others (2008), and recalling the L-moments in terms of the parameters and letting \eqn{\eta = \exp(-\beta\psi)} are
\deqn{\lambda_1 = \frac{1 - \eta + \eta\log(\eta)}{\beta(1-\eta)}\mbox{,}}
\deqn{\lambda_2 = \frac{1 + 2\eta\log(\eta) - \eta^2}{2\beta(1-\eta)^2}\mbox{, and}}
\deqn{\tau_2 = \frac{\lambda_2}{\lambda_1} = \frac{1 + 2\eta\log(\eta) - \eta^2}{2(1-\eta)[1-\eta+\eta\log(\eta)]}\mbox{,}}
and \eqn{\tau_2} is a monotonic function of \eqn{\eta} is decreasing from \eqn{\tau_2 = 1/2} at \eqn{\eta = 0} to \eqn{\tau_2 = 1/3} at \eqn{\eta = 1} the parameters are readily solved given \eqn{\tau_2 = [1/3, 1/2]}, the \R function \code{uniroot}  can be used to solve for \eqn{\eta} with a starting interval of \eqn{(0, 1)}, then the parameters in terms of the parameters are
\deqn{\alpha = \frac{1 - \eta + \eta\log(\eta)}{(1 - \eta)\lambda_1}\mbox{, and}}
\deqn{\psi = -\log(\eta)/\alpha\mbox{.}}

If the \eqn{\eta} is rooted as equaling zero, then it is assumed that \eqn{\hat\tau_2 \equiv \tau_2} and the exponential distribution triggered, or if the \eqn{\eta} is rooted as equaling unity, then it is assumed that \eqn{\hat\tau_2 \equiv \tau_2} and the uniform distribution triggered (see below).

The distribution is restricted to a narrow range of L-CV (\eqn{\tau_2 = \lambda_2/\lambda_1}). If \eqn{\tau_2 = 1/3}, the process represented is a stationary Poisson for which the probability density function is simply the uniform distribution and \eqn{f(x) = 1/\psi}. If \eqn{\tau_2 = 1/2}, then the distribution is represented as the usual exponential distribution with a location parameter of zero and a scale parameter \eqn{1/\beta}. Both of these limiting conditions are supported.

If the distribution shows to be uniform (\eqn{\tau_2 = 1/3}), then the third element in the returned parameter vector is used as the \eqn{\psi} parameter for the uniform distribution, and the first and second elements are \code{NA} of the returned parameter vector.

If the distribution shows to be exponential (\eqn{\tau_2 = 1/2}), then the second element in the returned parameter vector is the inverse of the rate parameter for the exponential distribution, and the first element is \code{NA} and the third element is \code{0} (a numeric \code{FALSE}) of the returned parameter vector.
}
\usage{
partexp(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or
\code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{texp}.}
  \item{para}{The parameters of the distribution.}
  \item{ifail}{A logical value expressed in numeric form indicating the failure or success state of the parameter estimation. A value of two indicates that the \eqn{\tau_2 < 1/3} whereas a value of three indicates that the \eqn{\tau_2 > 1/2}; for each of these inequalities a fuzzy tolerance of one part in one million is used. Successful parameter estimation, which includes the uniform and exponential boundaries, is indicated by a value of zero.}
  \item{ifail.message}{Various messages for successful and failed parameter estimations are reported. In particular, there are two conditions for which each distributional boundary (uniform or exponential) can be obtained. First, for the uniform distribution, one message would indicate if the \eqn{\tau_2 = 1/3} is assumed within a one part in one million will be identified or if \eqn{\eta} is rooted to 1. Second, for the exponential distribution, one message would indicate if the \eqn{\tau_2 = 1/2} is assumed within a one part in one million will be identified or if \eqn{\eta} is rooted to 0.}
  \item{eta}{The value for \eqn{\eta}. The value is set to either unity or zero if the \eqn{\tau_2} fuzzy tests as being equal to 1/3 or 1/2, respectively. The value is set to the rooted value of \eqn{\eta} for all other valid solutions. The value is set to \code{NA} if \eqn{\tau_2} tests as being outside the 1/3 and 1/2 limits.}
  \item{source}{The source of the parameters: \dQuote{partexp}.}
}
\references{
Vogel, R.M., Hosking, J.R.M., Elphick, C.S., Roberts, D.L., and Reed, J.M., 2008, Goodness of fit of probability distributions for sightings as species approach extinction: Bulletin of Mathematical Biology, DOI 10.1007/s11538-008-9377-3, 19 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomtexp}}, \code{\link{cdftexp}}, \code{\link{pdftexp}}, \code{\link{quatexp}}
}
\examples{
# truncated exponential is a nonstationary poisson process
A  <- partexp(vec2lmom(c(100, 1/2),   lscale=FALSE)) # pure exponential
B  <- partexp(vec2lmom(c(100, 0.499), lscale=FALSE)) # almost exponential
BB <- partexp(vec2lmom(c(100, 0.45),  lscale=FALSE)) # truncated exponential
C  <- partexp(vec2lmom(c(100, 1/3),   lscale=FALSE)) # stationary poisson process
D  <- partexp(vec2lmom(c(100, 40))) # truncated exponential
}
\keyword{distribution (parameters)}
\keyword{Distribution: Exponential (trimmed)}
