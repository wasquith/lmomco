\name{par2pdf}
\alias{par2pdf}
\title{Probability Density Function of the Distributions}
\description{
This function acts as a frontend or dispatcher to the distribution-specific probability density functions.
}
\usage{
par2pdf(f, para, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
  \item{...}{The additional arguments are passed to the quantile function such as \cr \code{paracheck = FALSE} for the Generalized Lambda distribution (\code{\link{quagld}}).}
}
\value{
  Quantile value for \eqn{F}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{par2cdf}}, \code{\link{par2qua}} }
\examples{
lmr     <- lmoms(rnorm(20))
para    <- parnor(lmr)
median  <- par2qua(0.5,para)
}
\keyword{probability density function}
