\name{quagdd}
\alias{quagdd}
\title{Quantile Function of the Gamma Difference Distribution}
\description{
This function computes the quantiles of the Gamma Difference distribution (Klar, 2015) given parameters (\eqn{\alpha_1 > 0}, \eqn{\beta_1 > 0}, \eqn{\alpha_2 > 0}, \eqn{\beta_2 > 0}) computed by \code{\link{pargdd}}. The quantile function requires numerical rooting of the cumulative distribution function \code{\link{cdfgdd}}.
}
\usage{
quagdd(f, para, paracheck=TRUE, silent=TRUE, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargdd}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity.}
  \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}.}
  \item{...}{Additional arguments to pass.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\author{W.H. Asquith}
\references{
Klar, B., 2015, A note on gamma difference distributions: Journal of Statistical Computation and Simulation v. 85, no. 18, pp. 1--8, \doi{10.1080/00949655.2014.996566}.
}
\seealso{\code{\link{cdfgdd}}, \code{\link{pdfgdd}}, \code{\link{lmomgdd}}, \code{\link{pargdd}}}
\examples{
para <- list(para=c(3, 0.1, 0.1, 4), type="gdd")
quagdd(0.5, para) # [1] 26.71568
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Gamma Difference}
