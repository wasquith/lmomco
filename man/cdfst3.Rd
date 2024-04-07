\name{cdfst3}
\alias{cdfst3}
\title{Cumulative Distribution Function of the 3-Parameter Student t Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the  3-parameter Student t distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\nu}) computed by \code{\link{parst3}}.  There is no explicit solution for the cumulative distribution function for value \code{X} but built-in \R functions can be used. For \code{U} = \eqn{\xi} and \code{A} = \eqn{\alpha} and  for \eqn{1.001 \le \nu \le 10^5.5}, one can use \code{pt((X-U)/A, N)} for \code{N} = \eqn{\nu}. The \R function \code{pt} is used for the 1-parameter Student t cumulative distribution function. The limits for \eqn{\nu} stem from study of ability for theoretical integration of the quantile function to produce viable \eqn{\tau_4} and \eqn{\tau_6} (see \code{inst/doc/t4t6/studyST3.R}).
}
\usage{
cdfst3(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parst3}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical on whether the parameter should be check for validity.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{pdfst3}}, \code{\link{quast3}}, \code{\link{lmomst3}}, \code{\link{parst3}}
}
\examples{
  lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
  cdfst3(191.5143, parst3(lmr)) # 75th percentile
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Student t (3-parameter)}
