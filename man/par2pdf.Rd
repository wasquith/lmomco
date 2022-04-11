\name{par2pdf}
\alias{par2pdf}
\title{Probability Density Function of the Distributions}
\description{
This function acts as a frontend or dispatcher to the distribution-specific probability density functions.
}
\usage{
par2pdf(x, para, ...)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
  \item{...}{The additional arguments are passed to the quantile function such as \cr \code{paracheck = FALSE} for the Generalized Lambda distribution (\code{\link{quagld}}).}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\author{W.H. Asquith}
\seealso{\code{\link{par2cdf}}, \code{\link{par2qua}} }
\examples{
para    <- parnor(lmoms(rnorm(20)))
density <- par2pdf(par2qua(0.5, para), para)
}
\keyword{probability density function}
