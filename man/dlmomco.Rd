\name{dlmomco}
\alias{dlmomco}
\title{Probability Density Function of the Distributions}
\description{
This function acts as an alternative front end to \code{\link{par2pdf}}. The nomenclature of the \code{\link{dlmomco}} function is to mimic that of built-in \R functions that interface with distributions.
}
\usage{
dlmomco(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
}
\value{
  Probability density for \code{x}.
}
\author{W.H. Asquith}
\seealso{\code{\link{plmomco}}, \code{\link{qlmomco}}, \code{\link{rlmomco}}, \code{\link{slmomco}} }
\examples{
para <- vec2par(c(0,1),type="nor") # standard normal parameters
nonexceed <- dlmomco(1,para) # percentile of one standard deviation
}
\keyword{probability density function}
\concept{probability density}
\keyword{The lmomco functions}
\keyword{The lmomco function mimics of R nomenclature}

