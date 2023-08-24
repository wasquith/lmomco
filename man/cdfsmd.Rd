\name{cdfsmd}
\alias{cdfsmd}
\title{Cumulative Distribution Function of the Singh--Maddala Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Singh--Maddala (Burr Type XII) distribution given parameters (\eqn{a}, \eqn{b}, and \eqn{q}) of the distribution computed by \code{\link{parsmd}}. The cumulative distribution function is
\deqn{F(x) = 1 - \biggl(1 + \bigl[ (x - \xi) / a \bigr]^b \biggl)^{-q}\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x} with \eqn{0 \le x \le \infty}, \eqn{\xi} is a location parameter, \eqn{a} is a scale parameter (\eqn{a > 0}), \eqn{b} is a shape parameter (\eqn{b > 0}), and \eqn{q} is another shape parameter (\eqn{q > 0}).
}
\usage{
cdfsmd(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parsmd}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Kumar, D., 2017, The Singh--Maddala distribution---Properties and estimation: International Journal of System Assurance Engineering and Management, v. 8, no. S2, 15 p., \doi{10.1007/s13198-017-0600-1}.

Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfsmd}}, \code{\link{quasmd}}, \code{\link{lmomsmd}}, \code{\link{parsmd}}}
\examples{
# The SMD approximating the normal and use x=0
tau4_of_normal <- 30 * pi^-1 * atan(sqrt(2)) - 9 # from theory
cdfsmd(0, parsmd( vec2lmom( c( -pi, pi, 0, tau4_of_normal ) ) ) ) # 0.7138779
pnorm( 0, mean=-pi, sd=pi*sqrt(pi))                               # 0.7136874

\dontrun{
t3 <- 0.6
t4 <- (t3 * (1 + 5 * t3))/(5 + t3) # L-kurtosis of GPA from lmrdia()
paraA <- parsmd( vec2lmom( c( -1000, 200, t3, t4-0.02 ) ) )
paraB <- parsmd( vec2lmom( c( -1000, 200, t3, t4      ) ) )
paraC <- parsmd( vec2lmom( c( -1000, 200, t3, t4+0.02 ) ) )
FF <- nonexceeds(); x <- quasmd(FF, paraA)
plot( x, prob2grv(cdfsmd(x, paraA)), col="red", type="l",
      xlab="Quantile", ylab="Gumbel Reduced Variate, prob2grv()")
lines(x, prob2grv(cdfsmd(x, paraB)), col="green")
lines(x, prob2grv(cdfsmd(x, paraC)), col="blue" ) # }
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
