\name{lmomsmd}
\alias{lmomsmd}
\title{L-moments of the Singh--Maddala Distribution}
\description{
This function computes the L-moments of the Singh--Maddala (Burr Type XII) distribution given the parameters (\eqn{\xi}, \eqn{a}, \eqn{b}, and \eqn{q}) from \code{\link{parsmd}}. The first L-moment (\eqn{\lambda_1}) for \eqn{b' = 1/b} and \eqn{R = a\Gamma(1 + b')} is
\deqn{\lambda_1 = R\times\biggl[\frac{a\Gamma(1q-b')}{\Gamma(1q)}\biggr] + \xi\mbox{.}}

The second L-moment (\eqn{\lambda_2}) is
\deqn{\lambda_2 = R\times\biggl[\frac{1\Gamma(1q - b')}{\Gamma(1q)} -
                                \frac{1\Gamma(2q - b')}{\Gamma(2q)}\biggr]\mbox{.}}

The third L-moment (\eqn{\lambda_3}) is
\deqn{\lambda_3 = R\times\biggl[\frac{1\Gamma(1q - b')}{\Gamma(1q)} -
                                \frac{3\Gamma(2q - b')}{\Gamma(2q)} +
                                \frac{2\Gamma(3q - b')}{\Gamma(3q)}\biggr]\mbox{.}}

The fourth L-moment (\eqn{\lambda_4}) is
\deqn{\lambda_4 = R\times\biggl[\frac{ 1\Gamma(1q - b')}{\Gamma(1q)} -
                                \frac{ 6\Gamma(2q - b')}{\Gamma(2q)} +
                                \frac{10\Gamma(3q - b')}{\Gamma(3q)} -
                                \frac{ 5\Gamma(4q - b')}{\Gamma(4q)}\biggr]\mbox{.}}

The fifth L-moment (\eqn{\lambda_5}) (unique to \pkg{lmomco} development) is
\deqn{\lambda_5 = R\times\biggl[\frac{ 1\Gamma(1q - b')}{\Gamma(1q)} -
                                \frac{10\Gamma(2q - b')}{\Gamma(2q)} +
                                \frac{30\Gamma(3q - b')}{\Gamma(3q)} -
                                \frac{35\Gamma(4q - b')}{\Gamma(4q)} +
                                \frac{14\Gamma(5q - b')}{\Gamma(5q)}\biggr]\mbox{.}}

The sixth L-moment (\eqn{\lambda_6}) (unique to \pkg{lmomco} development) is
\deqn{\lambda_6 = R\times\biggl[\frac{  1\Gamma(1q - b')}{\Gamma(1q)} -
                                \frac{ 15\Gamma(2q - b')}{\Gamma(2q)} +
                                \frac{ 70\Gamma(3q - b')}{\Gamma(3q)} -
                                \frac{140\Gamma(4q - b')}{\Gamma(4q)} +}
\deqn{\frac{126\Gamma(5q - b')}{\Gamma(5q)} -
      \frac{ 42\Gamma(6q - b')}{\Gamma(6q)}\biggr]\mbox{.}}
}
\usage{
lmomsmd(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
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
source of the L-moments: \dQuote{lmomsmd}.}
}
\references{
Bhatti, F.A., Hamedani, G.G., Korkmaz, M.\enc{C}{C}., and Munir Ahmad, M., 2019, New modified Singh--Maddala distribution---Development, properties, characterizations, and applications: Journal of Data Science, v. 17, no. 3, pp. 551--574, \doi{10.6339/JDS.201907_17(3).0006}.

Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\author{W.H. Asquith}
\seealso{\code{\link{parsmd}}, \code{\link{cdfsmd}}, \code{\link{pdfsmd}}, \code{\link{quasmd}} }
\examples{
lmr <- lmoms(c(123, 34, 4, 654, 37, 78), nmom=6)
lmr$source <- lmr$trim <- lmr$leftrim <- lmr$rightrim <-NULL
# The parsmd() reports Tau4 is too big and snaps it to an empirical boundary.
# "Tau4(~Tau3) snapped to upper limit, Tau4=0.65483 for Tau3=0.75126"
bmr <- lmomsmd(parsmd(lmr, snap.tau4=TRUE))
dmr <- data.frame(bmr$lambdas, bmr$ratios)
cbind(as.data.frame(lmr), dmr) # See in table that row 4 has different Tau4s
#  lambdas    ratios bmr.lambdas bmr.ratios
# 1   155.0        NA   155.00000         NA
# 2   118.6 0.7651613   118.60000  0.7651613
# 3    89.1 0.7512648    89.18739  0.7520016
# 4    82.1 0.6922428    77.59904  0.6542921 # see different Tau4s (snapping)
# 5    69.5 0.5860034    68.40150  0.5767411 # We are not fitting to these
# 6   102.5 0.8642496    62.58792  0.5277228 # higher L-moments.

# T3 and T4 of the Gumbel distribution, which is inside the SMD domain.
gumt3t4 <- c(log(9/8)/log(2), (16 * log(2) - 10 * log(3))/log(2))
lmr <- theoLmoms(pargum(vec2lmom(c(155, 118.6, gumt3t4))), nmom=6)
lmr$source <- lmr$trim <- lmr$leftrim <- lmr$rightrim <-NULL
bmr <- lmomsmd(parsmd(lmr, snap.tau4=TRUE))
dmr <- data.frame(bmr$lambdas, bmr$ratios)
cbind(as.data.frame(lmr), dmr)
#      lambdas     ratios bmr.lambdas bmr.ratios
# 1 155.000000         NA  155.000000         NA
# 2 118.600005 0.76516132  118.600005  0.7651613
# 3  20.153103 0.16992498   20.153104  0.1699250
# 4  17.834464 0.15037490   17.834464  0.1503749 # see same Tau4s (no snapping)
# 5   6.625972 0.05586823    7.688957  0.0648310 # We are not fitting to these
# 6   6.891842 0.05810997    7.213039  0.0608182 # higher L-moments.

\dontrun{
  # T3 and T4 of the Gumbel distribution, which is inside the SMD domain.
  gumt3t4 <- c(log(9/8)/log(2), (16 * log(2) - 10 * log(3))/log(2))
  FF <- nonexceeds(); qFF <- qnorm(FF)
  gumx <- qlmomco(FF, pargum(vec2lmom(c(155, 118.6, gumt3t4))))
  smdx <- qlmomco(FF, parsmd(lmr, snap.tau4=TRUE))
  plot( qFF, gumx, col="blue", type="l",
       xlab="Standard normal variate", ylab="Quantile")
  lines(qFF, smdx, col="red") # }
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
