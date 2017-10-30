\name{lmom2vec}
\alias{lmom2vec}
\title{ Convert an L-moment object to a Vector of L-moments}
\description{
This function converts an L-moment object in the structure used by \pkg{lmomco} into a simple vector. The precise operation of this function is dependent on the L-moment object argument. The \code{\link{lmorph}} function is not used. This function is useful if one needs to use certain functions in the \pkg{lmoms} package that are built around vectors of L-moments and L-moment ratios as arguments.
}
\usage{
lmom2vec(lmom, ...)
}
\arguments{
  \item{lmom}{L-moment object as from functions such as \code{lmoms}, \code{lmom.ub}, and \code{vec2lmom}.}
  \item{...}{Not presently used.}
}
\value{
A vector of the L-moments (\eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\tau_3}, \eqn{\tau_4}, \eqn{\tau_5}, \dots, \eqn{\tau_r}).
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lmom.ub}}, \code{\link{lmoms}}, \code{\link{lmorph}}, \code{\link{vec2lmom}}, \code{\link{pwm2vec}}
}
\examples{
  lmr <- lmoms(rnorm(40))
  lmom2vec(lmr)
  lmr <- vec2lmom(c(140,150,.3,.2,-.1))
  lmom2vec(lmr)
}
\keyword{utility (vector)}
\keyword{L-moment (utility)}
\keyword{utility (L-moments to vector of)}

