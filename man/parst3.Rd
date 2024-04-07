\name{parst3}
\alias{parst3}
\title{Estimate the Parameters of the 3-Parameter Student t Distribution}
\description{
This function estimates the parameters of the 3-parameter Student t distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomst3}}. The largest value of \eqn{\nu} recognized is \eqn{10^5.5}, which is the near the Normal distribution and the smallest value recognized is \eqn{1.001}, which is near the Cauchy. As \eqn{\nu \rightarrow 1} the distribution limits to the Cauchy, but the implementation here does not switch over to the Cauchy. Therefore in \pkg{lmomco} \eqn{1.001 \le \nu \le 10^5.5}. The \eqn{\nu} is the \dQuote{degrees of freedom} parameter that is well-known with the 1-parameter Student t distribution. The \eqn{nu} limits are studied in the \code{inst/doc/t4t6/studyST3.R} script and the \code{\link{theoTLmoms}} function and its performance on the quantile function of the distribution provide the guidance including range of numerically computable \eqn{\tau_6}. The \eqn{\tau_4} value can be set as low as \eqn{0.1226} as short-hand for the true lower L-kurtosis limit, which is that of the Normal (\eqn{30/\pi \times \mathrm{atan}(\sqrt{2}) - 9 = 0.1226017} and additional decimals). Internally, the given \eqn{0.1226 \le \tau_4 \le 0.1226017} is snapped to that of the Normal with an internal small positive nudge up. The \eqn{\tau_4 > 0.998} are set to \eqn{\tau_4 = 0.998}.
}
\usage{
parst3(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{st3}.}
  \item{para}{The parameters of the distribution.}
  \item{rt}{The returned list of the \code{uniroot()} call to estimate \eqn{\nu}.}
  \item{source}{The source of the parameters: \dQuote{parst3}.}
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{lmomst3}}, \code{\link{cdfst3}}, \code{\link{pdfst3}}, \code{\link{quast3}}
}
\examples{
  parst3(vec2lmom(c(10, 2, 0, 0.1226)))$para
  parst3(vec2lmom(c(10, 2, 0, 0.14  )))$para
  parst3(vec2lmom(c(10, 2, 0, 0.4   )))$para
  parst3(vec2lmom(c(10, 2, 0, 0.9   )))$para
  parst3(vec2lmom(c(10, 2, 0, 0.998 )))$para
}
\keyword{distribution (parameters)}
\keyword{Distribution: Student t (3-parameter)}
