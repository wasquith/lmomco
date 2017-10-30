\name{qlmomco}
\alias{qlmomco}
\title{Quantile Function of the Distributions}
\description{
This function acts as an alternative front end to \code{\link{par2qua}}. The nomenclature of the \code{\link{qlmomco}} function is to mimic that of built-in \R functions that interface with distributions.
}
\usage{
qlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
}
\value{
  Quantile value for \eqn{F} for the specified parameters.
}
\author{ W.H. Asquith}
\seealso{\code{\link{dlmomco}}, \code{\link{plmomco}}, \code{\link{rlmomco}}, \code{\link{slmomco}}, \code{\link{add.lmomco.axis}}, \code{\link{supdist}}
}
\examples{
para <- vec2par(c(0,1),type='nor') # standard normal parameters
p75  <- qlmomco(.75,para) # 75th percentile of one standard deviation
}
\keyword{quantile function}
\concept{quantile}
\keyword{The lmomco functions}
\keyword{The lmomco function mimics of R nomenclature}
