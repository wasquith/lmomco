\name{hlmomco}
\alias{hlmomco}
\title{Hazard Functions of the Distributions}
\description{
This function acts as a front end to \code{\link{dlmomco}} and \code{\link{plmomco}} to compute the hazard function \eqn{h(x)} or conditional failure rate. The function is defined by
\deqn{h(x) = \frac{f(x)}{1 - F(x)}\mbox{,}}
where \eqn{f(x)} is a probability density function and \eqn{F(x)} is the cumulative distribution function.

To help with intuitive understanding of what \eqn{h(x)} means (Ugarte and others, 2008), let \eqn{\mathrm{d}x} represent a small unit of measurement. Then the quantity \eqn{h(x)\,\mathrm{d}x} can be conceptualized as the approximate probability that random variable \eqn{X} takes on a value in the interval \eqn{[x, x+\mathrm{d}x]}.

Ugarte and others (2008) continue by stating that \eqn{h(x)} represents the instantaneous rate of death or failure at time \eqn{x}, given the survival to time \eqn{x} has occurred. Emphasis is needed that \eqn{h(x)} is a rate of probability change and not a probability itself.
}
\usage{
hlmomco(x,para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
}
\value{
  Hazard rate for \code{x}.
}
\note{
The hazard function is numerically solved for the given cumulative distribution and probability density functions and not analytical expressions for the hazard function that do exist for many distributions.
}
\references{
Ugarte, M.D., Militino, A.F., and Arnholt, A.T., 2008, Probability and statistics with R: CRC Press, Boca Raton, FL.
}
\author{W.H. Asquith}
\seealso{\code{\link{plmomco}}, \code{\link{dlmomco}}}
\examples{
my.lambda <- 100
para <- vec2par(c(0,my.lambda), type="exp")

x <- seq(40:60)
hlmomco(x,para) # returns vector of 0.01
# because the exponential distribution has a constant
# failure rate equal to 1/scale or 1/100 as in this example.
}
\keyword{hazard function}
\keyword{lifetime/reliability analysis}
\concept{hazard}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
\keyword{The lmomco function mimics of R nomenclature}
