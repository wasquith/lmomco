\name{vec2lmom}
\alias{vec2lmom}
\title{Convert a Vector of L-moments to a L-moment Object}
\description{
This function converts a vector of L-moments to a L-moment object of \pkg{lmomco}. The object is an \R \code{list}. This function is intended to facilitate the use of L-moments (and TL-moments) that the user might have from other sources. L-moments and L-moment ratios of arbitrary length are supported.

Because in typical practice, the \eqn{k \ge 3} order L-moments are dimensionless ratios (\eqn{\tau_3}, \eqn{\tau_4}, and \eqn{\tau_5}), this function computes \eqn{\lambda_3}, \eqn{\lambda_4}, \eqn{\lambda_5} from \eqn{\lambda_2} from the ratios. However, typical practice is not set on the use of \eqn{\lambda_2} or \eqn{\tau} as measure of dispersion. Therefore, this function takes an \code{lscale} optional logical (\code{TRUE|FALSE}) argument---if \eqn{\lambda_2} is provided and \code{lscale=TRUE}, then \eqn{\tau} is computed by the function and if \eqn{\tau} is provided, then \eqn{\lambda_2} is computed by the function.
}
\usage{
vec2lmom(vec, lscale=TRUE,
         trim=NULL, leftrim=NULL, rightrim=NULL, checklmom=TRUE)
}
\arguments{
  \item{vec}{A vector of L-moment values in \eqn{\lambda_1}, \eqn{\lambda_2} or \eqn{\tau}, \eqn{\tau_3}, \eqn{\tau_4}, and \eqn{\tau_5} order.}
  \item{lscale}{A logical switch on the type of the second value of first argument.
                L-scale (\eqn{\lambda_2}) or LCV (\eqn{\tau}). Default is \code{TRUE},
                the second value in the first argument is \eqn{\lambda_2}.}
  \item{trim}{Level of symmetrical trimming, which should equal \code{NULL} if asymmetrical trimming is used.}
  \item{leftrim}{Level of trimming of the left-tail of the sample, which will equal \code{NULL} even if \code{trim = 1} if the trimming is symmetrical.}
  \item{rightrim}{Level of trimming of the right-tail of the sample, which will equal \code{NULL} even if \code{trim = 1} if the trimming is symmetrical.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default unless TL-moments are being constructed in lieu of using \code{\link{vec2TLmom}}.}
}
\value{
  An \R \code{list} is returned.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmoms}}, \code{\link{vec2pwm}}}
\examples{
lmr <- vec2lmom(c(12,0.6,0.34,0.20,0.05),lscale=FALSE)
}
\keyword{utility (vector)}

