\encoding{utf8}
\name{parpdq3}
\alias{parpdq3}
\title{Estimate the Parameters of the Polynomial Density-Quantile3 Distribution}
\description{
This function estimates the parameters of the Polynomial Density-Quantile3 distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between the distribution parameters and L-moments are seen under \code{\link{lmompdq3}}.
}
\usage{
parpdq3(lmom, checklmom=TRUE)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is unlikely that the L-moments will not be viable. However, for some circumstances or large simulation exercises then one might want to bypass this check.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{pdq3}.}
  \item{para}{The parameters of the distribution.}
  \item{ifail}{A numeric field connected to the \code{ifailtext}; a value of 0 indicates fully successful operation of the function.}
  \item{ifailtext}{A message, instead of a warning, about the internal operations or operational limits of the function.}
  \item{source}{The source of the parameters: \dQuote{parpdq3}.}
}
\note{
The following is a study of the performance of \code{parpdq3} as the upper limit of the shape parameter \eqn{\kappa} is approached. The algorithms have the ability to estimate the \eqn{\kappa} reliabily, it is the scale parameter \eqn{\alpha} that breaks down and hence there is a hard-wired setting of \eqn{|\kappa| > 0.98} in which a warning is issue in \code{parpdq3} about \eqn{\alpha} reliability:
\preformatted{
  A <- 10
  K <- seq(0.8, 1, by=0.0001)
  K <- sort(c(-K, K))
  As <- Ks <- rep(NA, length(K))
  for(i in 1:length(K)) {
    para <- list(para=c(0, A, K[i]), type="pdq3")
    As[i] <- parpdq3( lmompdq3(para) )$para[2]
    Ks[i] <- parpdq3( lmompdq3(para) )$para[3]
  }
  plot( K, (As-A)/A, type="l", col="red")
  abline(v=c(-0.98, +0.98)) # heuristically determined threshold
}
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmompdq3}}, \code{\link{cdfpdq3}}, \code{\link{pdfpdq3}}, \code{\link{quapdq3}}
}
\examples{
para <- list(para=c(0, 0.4332, -0.7029), type="pdq3")
parpdq3(lmompdq3(para))$para

para <- list(para=c(0, 0.4332, 0.7029), type="pdq3")
parpdq3(lmompdq3(para))$para

para <- list(para=c(0, 0.4332, 1-sqrt(.Machine$double.eps)), type="pdq3")
parpdq3(lmompdq3(para))$para

para <- list(para=c(0, 0.4332, -1+sqrt(.Machine$double.eps)), type="pdq3")
parpdq3(lmompdq3(para))$para

para <- list(para=c(0, 0.4332, +0.0001), type="pdq3")
parpdq3(lmompdq3(para))$para

para <- list(para=c(0, 0.4332, -0.0001), type="pdq3")
parpdq3(lmompdq3(para))$para

para <- list(para=c(0, 0.4332, 0), type="pdq3")
parpdq3(lmompdq3(para))$para
}
\keyword{distribution (parameters)}
\keyword{Distribution: Polynomial Density-Quantile3}
