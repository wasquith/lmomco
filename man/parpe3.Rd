\name{parpe3}
\alias{parpe3}
\title{Estimate the Parameters of the Pearson Type III Distribution}
\description{
This function estimates the parameters of the Pearson Type III distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The L-moments in terms of the parameters are complicated and solved numerically. For the implementation in \pkg{lmomco}, the three parameters are \eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma} for the mean, standard deviation, and skew, respectively.
}
\usage{
parpe3(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{pe3}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parpe3}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{
\code{\link{lmompe3}}, \code{\link{cdfpe3}}, \code{\link{pdfpe3}}, \code{\link{quape3}}
}
\examples{
lmr <- lmoms(rnorm(20))
parpe3(lmr)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Pearson Type III}
