\name{parsmd}
\alias{parsmd}
\title{Estimate the Parameters of the Singh--Maddala Distribution}
\description{
This function estimates the parameters of the Singh--Maddala (Burr Type XII) distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The L-moments in terms of the parameters are complicated and solved numerically.
}
\usage{
parsmd(lmom, checklmom=TRUE, emplims=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_3} and \eqn{\tau_4} inequality, \code{\link{are.lmom.valid}}). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{emplims}{A logical to trigger the application of the empirical limits of the distribution in terms of \eqn{\tau_3} and \eqn{\tau_5} wherein parameter estimation appears numerically possible and such parameters return the given values of these L-moment ratios.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{smd}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parsmd}.}
}
\references{
Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomsmd}}, \code{\link{cdfsmd}}, \code{\link{pdfsmd}}, \code{\link{quasmd}}
}
\examples{
lmr <- lmoms(rnorm(20))
parsmd(lmr)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
