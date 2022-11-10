\name{par2cdf2}
\alias{par2cdf2}
\title{Equivalent Cumulative Distribution Function of Two Distributions}
\description{
This function computes the nonexceedance probability of a given quantile from a linear weighted combination of two quantile functions but accomplishes this from the perspective of cumulative distribitution functions (see \code{\link{par2qua2}}).  For the current implementation simply \code{uniroot}'ing of a internally declared function and \code{\link{par2qua2}} is made. Mathematical details are provided under \code{\link{par2qua2}}.
}
\usage{
par2cdf2(x, para1, para2, weight=NULL, ...)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para1}{The first distribution parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{para2}{The second distribution parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{weight}{An optional weighting argument to use in lieu of the \code{F}. Consult the documentation for \code{\link{par2qua2}} for the implementation details when \code{weight} is \code{NULL}.}
  \item{...}{The additional arguments are passed to the quantile function.}
}
\value{
  Nonexceedance probabilities (\eqn{0 \le F \le 1}) for \code{x} from the two distributions.
}
\author{W.H. Asquith}
\seealso{\code{\link{par2cdf}}, \code{\link{lmom2par}}, \code{\link{par2qua2}}
}
\examples{
lmr <- lmoms(rnorm(20)); left <- parnor(lmr); right <- pargev(lmr)
mixed.median    <- par2qua2(0.5,          left, right)
mixed.nonexceed <- par2cdf2(mixed.median, left, right)
}
\keyword{distribution (mixture)}
\keyword{distribution (mixture)}
\keyword{mixed distribution}
