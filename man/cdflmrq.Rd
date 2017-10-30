\name{cdflmrq}
\alias{cdflmrq}
\title{Cumulative Distribution Function of the  Linear Mean Residual Quantile Function Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the \dQuote{Linear Mean Residual Quantile Function} distribution given parameters computed by \code{\link{parlmrq}}.  The cumulative distribution function  has no explicit form and requires numerical methods. The \R function \code{uniroot} is used to root the quantile function \code{\link{qualmrq}} to compute the nonexceedance probability. The function returns 0 or 1 if the \code{x} argument is at or beyond the limits of the distribution as specified by the parameters. The \code{\link{cdflmrq}} function is also used with numerical methods to solve the \code{\link{pdflmrq}}.
}
\usage{
cdflmrq(x, para, paracheck=FALSE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parlmrq}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}. This switch is made so that the root solution needed for \code{\link{cdflmrq}} exhibits an extreme speed increase because of the repeated calls to \code{qualmrq}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Midhu, N.N., Sankaran, P.G., and Nair, N.U., 2013, A class of distributions with linear mean residual quantile function and it's generalizations: Statistical Methodology, v. 15, pp. 1--24.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdflmrq}}, \code{\link{qualmrq}}, \code{\link{lmomlmrq}}, \code{\link{parlmrq}}
}
\examples{
  lmr <- lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2))
  cdflmrq(2,parlmrq(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Linear Mean Residual Quantile Function}
