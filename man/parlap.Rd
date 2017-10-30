\name{parlap}
\alias{parlap}
\title{Estimate the Parameters of the Laplace Distribution}
\description{
This function estimates the parameters of the Laplace distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and sample L-moments are simple, but there are two methods. The first method, which is the only one implemented in \pkg{lmomco}, jointly uses \eqn{\lambda_1, \lambda_2, \lambda_3}, and \eqn{\lambda_4}. The mathematical expressions are
\deqn{\xi = \lambda_1 - 50/31\times\lambda_3 \mbox{and}}
\deqn{\alpha = 1.4741\lambda_2 - 0.5960\lambda_4 \mbox{.}}

The alternative and even simpler method only uses \eqn{\lambda_1} and \eqn{\lambda_2}. The mathematical expressions are
\deqn{\xi = \lambda_1\mbox{\, and}}
\deqn{\alpha = \frac{4}{3}\lambda_2\mbox{.}}
The user could easily estimate the parameters from the L-moments and use \code{\link{vec2par}} to create a parameter object.
}
\usage{
parlap(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{lap}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parlap}.}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: IBM Research Report RC12210, T.J. Watson Research Center, Yorktown Heights, New York.
}
\note{
The decision to use only one of the two systems of equations for Laplace fitting is largely arbitrary, but it seems most fitting to use four L-moments instead of two.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomlap}},
\code{\link{cdflap}}, \code{\link{pdflap}}, \code{\link{qualap}}
}
\examples{
lmr <- lmoms(rnorm(20))
parlap(lmr)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Laplace}
