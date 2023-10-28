\name{f2flo}
\alias{f2flo}
\title{Conversion of Annual Nonexceedance Probability to Conditional Probability Nonexceedance Probabilities}
\description{
This function converts the cumulative distribution function of \eqn{F(x)} to a conditional cumulative distribution function \eqn{P(x)} based on the probability level of the left-hand threshold. It is recommended that this threshold (as expressed as a probability) be that value returned from \code{\link{x2xlo}} in element \code{pp}. The conversion is simple
\deqn{
P(x) <- (F(x) - pp)/(1-pp)\mathrm{,}
}
where the term \eqn{\mathrm{pp}} corresponds to the estimated probability or plotting position of the left-hand threshold.

This function is particularly useful for applications in which zero values in the data set require truncation so that logarithms of the data may be used. But also this function contributes to the isolation of the right-hand tail of the distribution for analysis. Finally, \code{f <- f[f >= pp]} for \code{type="ge"} or \code{f <- f[f > pp]} for \code{type="gt"} is used internally for probability subsetting, so the user does not have to do that with the nonexceedance probability before calling this function. The function \code{\link{f2f}} does similar subsetting without converting \eqn{F(x)} to \eqn{P(x)}. Users are directed to \bold{Examples} under \code{\link{par2qua2lo}} and carefully note how \code{f2flo} and \code{\link{f2f}} are used.
}
\usage{
f2flo(f, pp=NA, xlo=NULL, type=c("ge", "gt"))
}
\arguments{
  \item{f}{A vector of nonexceedance probabilities.}
  \item{pp}{The plotting position of the left-hand threshold and recommended to come from \code{\link{x2xlo}}.}
  \item{xlo}{An optional result from \code{\link{x2xlo}} from which the \code{pp} will be take instead of from the argument \code{pp}.}
  \item{type}{The type of the logical construction \code{gt} means greater than the \code{pp} and \code{ge} means greater than or equal to the \code{pp} for the computations. There can be subtle variations in conceptualization of the truncation need or purpose and hence this argument is provided for flexibility.}
}
\value{
  A vector of conditional nonexceedance probabilities.
}
\author{W.H. Asquith}
\seealso{\code{\link{x2xlo}}, \code{\link{flo2f}}, \code{\link{f2f}}, \code{\link{xlo2qua}} }
\examples{
# See examples for x2xlo().
}
\keyword{utility (probability)}
\keyword{conditional probability}
\keyword{probability (conditional)}
