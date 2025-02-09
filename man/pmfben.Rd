\encoding{utf8}
\name{pmfben}
\alias{pmfben}
\title{Probability Density Function of the Benford Distribution}
\description{
This function computes the probability mass function of the Benford distribution (Benford's Law) given parameters defining the number of first M-significant digits and the numeric base. The mass function has the simple expression
\deqn{
P(d) = \mathrm{log}_b\biggl(1 + \frac{1}{d}\biggr)\mbox{.}
}
for any base \eqn{b \ge 2} and digits \eqn{d}. The first significant digits in decimal are \eqn{d \in 1, \cdots, 9}, the first two-significant digits similarly are \eqn{d \in 10, \cdots, 99}, and the first three-significant digits similarly are \eqn{d \in 100, \cdots, 999}.
}
\usage{
pmfben(d, para=list(para=c(1, 10)), ...)
}
\arguments{
  \item{d}{A integer value vector of M-significant digits.}
  \item{para}{The number of first M-significant digits followed by the numerical base (only base10 supported) and the list structure mimics similar uses of the \pkg{lmomco} list structure. Default are the first significant digits and hence the digits 1 through 9.}
  \item{...}{Additional arguments to pass (not likely to be needed but changes in base handling might need this).}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Benford, F., 1938, The law of anomalous numbers: Proceedings of the American Philosophical Society, v. 78, no. 4, pp. 551--572, \url{https://www.jstor.org/stable/984802}.

Goodman, W., 2016, The promises and pitfalls of Benford’s law: Significance (Magazine), June 2015, pp. 38--41, \doi{10.1111/j.1740-9713.2016.00919.x}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfben}}, \code{\link{quaben}} }
\examples{
# probability masses matching values in authoritative texts
pmfben(1:9, para=list(para=c(1, 10)))
# [1] 0.30103000 0.17609126 0.12493874 0.09691001
# [5] 0.07918125 0.06694679 0.05799195 0.05115252
# [9] 0.04575749
cumsum( pmfben(1:9, para=list(para=c(1, 10))) ) # should end in unity
# [1] 0.3010300 0.4771213 0.6020600 0.6989700 0.7781513
# [6] 0.8450980 0.9030900 0.9542425 1.0000000
}
\keyword{distribution (discrete)}
\keyword{probability mass function}
\keyword{Distribution: Benford}
\concept{Benford Law}
