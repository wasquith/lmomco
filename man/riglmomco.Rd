\name{riglmomco}
\alias{riglmomco}
\title{Income Gap Ratio Quantile Function for the Distributions}
\description{
This function computes the Income Gap Ratio for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair et al. (2013, p. 230) as
\deqn{G(u) = 1 - \frac{{}_\mathrm{r}\lambda_1(u)}{x(u)}\mbox{,}}
where \eqn{G(u)} is the  income gap quantile for nonexceedance probability \eqn{u}, \eqn{x(u)} is a constant for \eqn{x(F = u)} is the quantile for \eqn{u}, and \eqn{{}_\mathrm{r}\lambda_1(u)} is the 1st reversed residual life L-moment (\code{\link{rreslife.lmoms}}).
}
\usage{
riglmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
   Income gap ratio quantile value for \eqn{F}.
}
\author{W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{rreslife.lmoms}}
}
\examples{
# Let us parametize some "income" distribution.
A <- vec2par(c(123, 264, 2.11), type="gov")
riglmomco(0.5, A)
\dontrun{
F <- nonexceeds(f01=TRUE)
plot(F, riglmomco(F,A), type="l",
     xlab="NONEXCEEDANCE PROBABILITY", ylab="INCOME GAP RATIO")}
}
\keyword{quantile function}
\keyword{income gap ratio quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
