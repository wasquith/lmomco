\name{lmomgam}
\alias{lmomgam}
\title{L-moments of the Gamma Distribution}
\description{
This function estimates the L-moments of the Gamma distribution given the parameters (\eqn{\alpha} and \eqn{\beta}) from \code{\link{pargam}}. The L-moments in terms of the parameters are complicated and solved numerically. This function is adaptive to the 2-parameter and 3-parameter Gamma versions supported by this package.  For legacy reasons, \pkg{lmomco} continues to use a port of Hosking's FORTRAN into R if the 2-parameter distribution is used but the 3-parameter generalized Gamma distribution calls upon \code{\link{theoLmoms.max.ostat}}. Alternatively, the \code{\link{theoTLmoms}} could be used: \code{theoTLmoms(para)} is conceptually equivalent to the internal calls to \code{\link{theoLmoms.max.ostat}} made for the \code{lmomgam} implementation.
}
\usage{
lmomgam(para, ...)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{...}{Additional arguments to pass to \code{\link{theoLmoms.max.ostat}}.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomgam}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, p. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pargam}}, \code{\link{cdfgam}}, \code{\link{pdfgam}}, \code{\link{quagam}} }
\examples{
lmomgam(pargam(lmoms(c(123,34,4,654,37,78))))

\dontrun{
# 3-p Generalized Gamma Distribution and comparisons of 3-p Gam parameterization.
#     1st parameter A[lmomco] = A[gamlss] =  exp(A[flexsurv])
#     2nd parameter B[lmomco] = B[gamlss] =      B[flexsurv]
#     3rd parameter C[lmomco] = C[gamlss] -->    C[flexsurv] = B[lmomco]/C[lmomco]
lmomgam(vec2par(c(7.4, 0.2, 14), type="gam"), nmom=5)$lambdas      # numerics
lmoms(gamlss.dist::rGG(50000, mu=7.4, sigma=0.2, nu=14))$lambdas   # simulation
lmoms(flexsurv::rgengamma(50000, log(7.4), 0.2, Q=0.2*14))$lambdas # simulation
#[1]  5.364557537  1.207492689 -0.110129217  0.067007941 -0.006747895
#[1]  5.366707749  1.209455502 -0.108354729  0.066360223 -0.006716783
#[1]  5.356166684  1.197942329 -0.106745364  0.069102821 -0.008293398#}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Gamma}
\keyword{Distribution: Generalized Gamma}
