\name{parwei}
\alias{parwei}
\title{Estimate the Parameters of the Weibull Distribution}
\description{
This function estimates the parameters of the Weibull distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The Weibull distribution is a reverse Generalized Extreme Value distribution.  As result, the Generalized Extreme Value algorithms are used for implementation of the Weibull in this package. The relations between the Generalized Extreme Value parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) and the Weibull parameters are
\deqn{\kappa = 1/\delta \mbox{,}}
\deqn{\alpha = \beta/\delta \mbox{, and}}
\deqn{\xi = \zeta - \beta \mbox{.}}
These relations are taken from Hosking and Wallis (1997). The relations between the distribution parameters and L-moments are seen under \code{\link{lmomgev}}.
}
\usage{
parwei(lmom,checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}

\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{wei}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parwei}.}
}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomwei}},
\code{\link{cdfwei}}, \code{\link{pdfwei}}, \code{\link{quawei}}
}
\examples{
parwei(lmoms(rnorm(20)))
\dontrun{
str(parwei(lmoms(rweibull(3000,1.3, scale=340)-1200))) #}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Weibull}
