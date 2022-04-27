\name{quagpa}
\alias{quagpa}
\title{Quantile Function of the Generalized Pareto Distribution}
\description{
This function computes the quantiles of the Generalized Pareto distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) computed by \code{\link{pargpa}}. The quantile function is
\deqn{x(F) = \xi + \frac{\alpha}{\kappa} \left( 1-(1-F)^\kappa \right)\mbox{,}}
for \eqn{\kappa \ne 0}, and
\deqn{x(F) = \xi - \alpha\log(1-F)\mbox{,}}
for \eqn{\kappa = 0}, where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and
\eqn{\kappa} is a shape parameter. The range of \eqn{x} is \eqn{\xi \le x \le \xi + \alpha/\kappa} if \eqn{k > 0}; \eqn{\xi \le x < \infty} if \eqn{\kappa \le 0}. Note that the shape parameter \eqn{\kappa} parameterization of the distribution herein follows that in tradition by the greater L-moment community and others use a sign reversal on \eqn{\kappa}. (The \pkg{evd} package is one example.)
}
\usage{
quagpa(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargpa}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124, \doi{10.1111/j.2517-6161.1990.tb01775.x}.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgpa}}, \code{\link{pdfgpa}}, \code{\link{lmomgpa}}, \code{\link{pargpa}}}
\examples{
  lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
  quagpa(0.5,pargpa(lmr))

\dontrun{
  # Let us compare L-moments, parameters, and 90th percentile for a simulated
  # GPA distibution of sample size 100 having the following parameters between
  # lmomco and lmom packages in R. The answers are the same.
  gpa.par <- lmomco::vec2par(c(1.02787, 4.54603, 0.07234), type="gpa")
  X <- lmomco::rlmomco(100, gpa.par)
   lmom::samlmu(X)
  lmomco::lmoms(X)
    lmom::pelgpa( lmom::samlmu(X))
  lmomco::pargpa(lmomco::lmoms(X))
    lmom::quagpa(0.90,   lmom::pelgpa(  lmom::samlmu(X)))
  lmomco::quagpa(0.90, lmomco::pargpa(lmomco::lmoms( X))) # }
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Generalized Pareto}
