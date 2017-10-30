\name{lmomgov}
\alias{lmomgov}
\title{L-moments of the Govindarajulu Distribution}
\description{
This function estimates the L-moments of the Govindarajulu distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\beta}) from \code{\link{pargov}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi + \frac{2\alpha}{\beta+2} \mbox{,}}
\deqn{\lambda_2 = \frac{2\alpha\beta}{(\beta+2)(\beta+3)} \mbox{,}}
\deqn{\tau_3 = \frac{\beta-2}{\beta+4} \mbox{, and}}
\deqn{\tau_4 = \frac{(\beta-5)(\beta-1)}{(\beta+4)(\beta+5)} \mbox{.}}

The limits of \eqn{\tau_3} are \eqn{(-1/2, 1)} for \eqn{\beta \rightarrow 0} and \eqn{\beta \rightarrow \infty}.
}
\usage{
lmomgov(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R list is returned.

    \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational source
of the L-moments: \dQuote{lmomgov}.}
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Nair, N.U., Sankaran, P.G., Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.

Nair, N.U., Sankaran, P.G., and Vineshkumar, B., 2012, The Govindarajulu distribution---Some Properties and applications: Communications in Statistics, Theory and Methods, v. 41, no. 24, pp. 4391--4406.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pargov}}, \code{\link{cdfgov}}, \code{\link{pdfgov}}, \code{\link{quagov}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmorph(lmr)
lmomgov(pargov(lmr))
\dontrun{
Bs <- exp(seq(log(.01),log(10000),by=.05))
T3 <- (Bs-2)/(Bs+4)
T4 <- (Bs-5)*(Bs-1)/((Bs+4)*(Bs+5))
plotlmrdia(lmrdia(), autolegend=TRUE)
points(T3, T4)
T3s <- c(-0.5,T3,1)
T4s  <- c(0.25,T4,1)
the.lm <- lm(T4s~T3s+I(T3s^2)+I(T3s^3)+I(T3s^4)+I(T3s^5))
lines(T3s, predict(the.lm), col=2)
max(residuals(the.lm))
summary(the.lm)
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Govindarajulu}
