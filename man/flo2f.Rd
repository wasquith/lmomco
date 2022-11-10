\name{flo2f}
\alias{flo2f}
\title{Conversion of Conditional Nonexceedance Probability to Annual Nonexceedance Probability}
\description{
This function converts  the conditional cumulative distribution function of \eqn{P(x)} to a cumulative distribution function \eqn{F(x)} based on the probability level of the left-hand threshold. It is recommended that this threshold (as expressed as a probability) be that value returned from \code{\link{x2xlo}} in attribute \code{pp}. The conversion is simple
\deqn{
F(x) = pp + (1 - pp)P(x)\mathrm{,}
}
where the term \eqn{pp} corresponds to the estimated probability or plotting position of the left-hand threshold.

This function is particularly useful for applications in which zero values in the data set require truncation so that logarithms of the data may be used. But also this function  contributes to the isolation of the right-hand tail of the distribution for analysis by conditionally trimming out the left-hand tail at the analyst's discretion.
}
\usage{
flo2f(f, pp=NA, xlo=NULL)
}
\arguments{
  \item{f}{A vector of nonexceedance probabilities.}
  \item{pp}{The plotting position of the left-hand threshold and recommended to come from \code{\link{x2xlo}}.}
  \item{xlo}{An optional result from \code{\link{x2xlo}} from which the \code{pp} will be take instead of from the argument \code{pp}.}
}
\value{
  A vector of converted nonexceedance probabilities.
}
\author{W.H. Asquith}
\seealso{\code{\link{x2xlo}}, \code{\link{f2flo}} }
\examples{
flo2f(f2flo(.73,pp=.1),pp=.1)
# Also see examples for x2xlo().
}
\keyword{utility (probability)}
\keyword{conditional probability}
\keyword{probability (conditional)}
