\name{partri}
\alias{partri}
\title{Estimate the Parameters of the Asymmetric Triangular Distribution}
\description{
This function estimates the parameters of the Asymmetric Triangular distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomtri}}.

The estimtion by the \code{partri} function is built around simultaneous numerical optimization of an objective function defined as
\deqn{\epsilon = \biggl(\frac{\lambda_1 - \hat\lambda_1}{\hat\lambda_1}\biggr)^2 + \biggl(\frac{\lambda_2 - \hat\lambda_2}{\hat\lambda_2}\biggr)^2 + \biggl(\frac{\tau_3 - \hat\tau_3}{1}\biggr)^2}
for estimation of the three parameters (\eqn{\nu}, minimum; \eqn{\omega}, mode; and \eqn{\psi}, maximum) from the sample L-moments (\eqn{\hat\lambda_1}, \eqn{\hat\lambda_2}, \eqn{\hat\tau_3}). The divisions shown in the objective function are used for scale removal to help make each L-moment order somewhat similar in its relative contribution to the solution. The coefficient of L-variation is not used because the distribution implementation by the \pkg{lmomco} package supports entire real number line and the loss of definition of \eqn{\tau_2} at \eqn{x = 0}, in particular, causes untidiness in coding.

The function is designed to support both left- or right-hand right triangular shapes because of (1) \code{paracheck} argument availability in \code{\link{lmomtri}}, (2) the sorting of the numerical estimates if the mode is no compatable with either of the limits, and (3) the snapping of \eqn{\nu = \omega \equiv (\nu^\star + \omega^\star)/2} when \eqn{\hat\tau_3 > 0.142857} or \eqn{\psi = \omega \equiv (\psi^\star + \omega^\star)/2} when \eqn{\hat\tau_3 < 0.142857} where the \eqn{\star} versions are the optimized values if the \eqn{\tau_3} is very near to its numerical bounds.
}
\usage{
partri(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}

\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{tri}.}
  \item{para}{The parameters of the distribution.}
  \item{obj.val}{The value of the objective function, which is the error of the optimization.}
  \item{source}{The source of the parameters: \dQuote{partri}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomtri}},
\code{\link{cdftri}}, \code{\link{pdftri}}, \code{\link{quatri}}
}
\examples{
lmr <- lmomtri(vec2par(c(10,90,100), type="tri"))
partri(lmr)

partri(lmomtri(vec2par(c(-11, 67,67), type="tri")))$para
partri(lmomtri(vec2par(c(-11,-11,67), type="tri")))$para
}
\keyword{distribution (parameters)}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}

