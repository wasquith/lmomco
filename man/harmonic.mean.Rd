\name{harmonic.mean}
\alias{harmonic.mean}
\title{The Harmonic Mean with Zero-Value Correction }
\description{
Compute the harmonic mean of a vector with a zero-value correction.
\deqn{\check{\mu} =
           \biggl(\frac{\sum^{N_T - N_0}_{i=1} 1/x_i}
                       {N_T - N_0}\biggr)^{-1} \times \frac{N_T - N_0}
                 {N_T} \mbox{,}}
where \eqn{\check{\mu}} is harmonic mean, \eqn{x_i} is a nonzero value of the data vector, \eqn{N_T} is the (total) sample size, \eqn{N_0} is the number of zero values.

}
\usage{
harmonic.mean(x)
}
\arguments{
  \item{x}{A vector of data values that will be reduced to non-missing values.}
}
\value{
  An \R \code{list} is returned.

  \item{harmean}{The harmonic mean with zero-value correction, \eqn{\check{\mu}}.}
  \item{correction}{The zero-value correction, \eqn{(N_T - N_0)/N_T}.}
  \item{source}{An attribute identifying the computational source of the harmonic mean: \dQuote{harmonic.mean}.}
}
\note{
The harmonic mean can not be computed when zero values are present. This situation is common in surface-water hydrology. As stated in the reference below, in order to calculate water-quality-based effluent limits (WQBELs) for human health protection, a harmonic mean flow is determined for all perennial streams and for streams that are intermittent with perennial pools. Sometimes these streams have days on which measured flow is zero. Because a zero flow cannot be used in the calculation of harmonic mean flow, the second term in the harmonic mean equation is an adjustment factor used to lower the harmonic mean to compensate for days on which the flow was zero. The zero-value correction is the same correction used by the EPA computer program DFLOW.
}
\references{
Texas Commission on Environmental Quality, 2003, Procedures to implement the Texas surface-water-quality standards: TCEQ RG--194, p. 47
}
\author{W.H. Asquith}

\seealso{\code{\link{pmoms}}}
\examples{
Q <- c(0,0,5,6,7)
harmonic.mean(Q)
}
\keyword{special function}

