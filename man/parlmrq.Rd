\name{parlmrq}
\alias{parlmrq}
\title{Estimate the Parameters of the Linear Mean Residual Quantile Function Distribution
}
\description{This function estimates the parameters of the Linear Mean Residual Quantile Function distribution given the L-moments of the data in an L-moment object such as that returned by \code{lmoms}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomlmrq}}.
}
\usage{
parlmrq(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{lmrq}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parlmrq}.}
}
\references{
Midhu, N.N., Sankaran, P.G., and Nair, N.U., 2013, A class of distributions with linear mean residual quantile function and it's generalizations: Statistical Methodology, v. 15, pp. 1--24.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomlmrq}}, \code{\link{cdflmrq}}, \code{\link{pdflmrq}}, \code{\link{qualmrq}}
}
\examples{
lmr <- lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2))
parlmrq(lmr)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Linear Mean Residual Quantile Function}
