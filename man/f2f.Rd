\name{f2f}
\alias{f2f}
\title{Subsetting of Nonexceedance Probabilities Related to Conditional Probability Adjustment}
\description{
This function subsetting nonexceedance probability according to
\deqn{
F(x) <- F(x | F(x) [>,\ge] p)\mathrm{,}
}
where \eqn{F} is nonexceedance probability for \eqn{x} and \code{pp} is the probability of a threshold. In \R logic, this is simply \code{f <- f[f > pp]} for \code{type == "gt"} or \code{f <- f[f >= pp]} for \code{type == "ge"}.

This function is particularly useful to shorten a commonly needed code logic related such as \code{FF[FF >= XloALL$pp]}, which would be needed in conditional probability adjustements and \code{XloALL} is from \code{\link{x2xlo}}. This could be replaced by syntax such as \code{f2f(FF, xlo=XloALL)}.  This function is very similar to \code{\link{f2flo}} with the only exception that the conditional probability adjustment is not made.
}
\usage{
f2f(f, pp=NA, xlo=NULL, type=c("ge", "gt"))
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
\seealso{\code{\link{x2xlo}}, \code{\link{f2flo}} }
\examples{
# See examples for x2xlo().
}
\keyword{utility (probability)}
\keyword{conditional probability}
\keyword{probability (conditional)}
