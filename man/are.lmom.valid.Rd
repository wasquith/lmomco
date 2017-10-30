\name{are.lmom.valid}
\alias{are.lmom.valid}
\title{Are the L-moments valid}
\description{
The L-moments have particular constraints on magnitudes and relation to each other.  This function evaluates and L-moment object whether the bounds for \eqn{\lambda_2 > 0} (L-scale), \eqn{|\tau_3| < 1} (L-skew), \eqn{\tau_4 < 1} (L-kurtosis), and \eqn{|\tau_5| < 1} are satisfied. An optional check on \eqn{\tau_4 \ge (5\tau_3^2 - 1)/4} is made. Also for further protection, the finitenesses of the mean (\eqn{\lambda_1}) and \eqn{\lambda_2} are also checked. These checks provide protection against say L-moments being computed on the logarithms of some data but the data themselves have values less than or equal to zero.

The TL-moments as implemented by the \code{TL} functions (\code{\link{TLmoms}}) are not applicable to the boundaries (well finiteness of course). The \code{are.lmom.valid} function should not be consulted on the TL-moments.
}
\usage{
are.lmom.valid(lmom, checkt3t4=TRUE)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}}, \code{\link{lmom.ub}}, \code{\link{pwm2lmom}}; and}
  \item{checkt3t4}{A logical triggering the above test on L-skew to L-kurtosis. This bounds in very small samples can be violated---usually the user will want this set and until (first release in 2017, v2.2.6) this bounds check was standard in \pkg{lmomco} for over a decade.}
}
\value{
  \item{TRUE}{L-moments are valid.}
  \item{FALSE}{L-moments are not valid.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M. and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmom.ub}}, \code{\link{lmoms}}, \code{\link{pwm2lmom}} }
\examples{
lmr <- lmoms(rnorm(20))
if(are.lmom.valid(lmr)) print("They are.")
\dontrun{
X <- c(1.7106278,  1.7598761,  1.2111335,  0.3447490,  1.8312889,
       1.3938445, -0.5376054, -0.2341009, -0.4333601, -0.2545229)
are.lmom.valid(lmoms(X))
are.lmom.valid(pwm2lmom(pwm.pp(X, a=0.5)))

# Prior to version 2.2.6, the next line could leak through as TRUE. This was a problem.
# Nonfiniteness of the mean or L-scale should have been checked; they are for v2.2.6+
are.lmom.valid(lmoms(log10(c(1,23,235,652,0)), nmom=1)) # of other nmom
}
}
\keyword{L-moment (theoretical)}
\keyword{L-moment (validation)}

