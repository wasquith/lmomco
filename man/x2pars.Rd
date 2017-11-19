\name{x2pars}
\alias{x2pars}
\title{Estimate an Ensemble of Parameters from Three Different Methods}
\description{
This function acts as a frontend to estimate an ensemble of parameters from the methods of L-moments (\code{\link{lmr2par}}), maximum likelihood (MLE, \code{\link{mle2par}}), and maximum product of spacings (MPS, \code{\link{mps2par}}). The parameters estimated by the L-moments are used as the \code{para.int} for the subsequent calls to MLE and MPS.
}
\usage{
x2pars(x, verbose=TRUE, ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{verbose}{A logical to control a sequential message ahead of each method.}
  \item{...}{The additional arguments, if ever used.}
}
\value{
  A \code{list} having
  \item{lmr}{Parameters from method of L-moments. This is expected to be \code{NULL} if the method fails, and the \code{NULL} is tested for in \code{\link{pars2x}}.}
  \item{mle}{Parameters from MLE. This is expected to be \code{NULL} if the method fails, and the \code{NULL} is tested for in \code{\link{pars2x}}.}
  \item{mps}{Parameters from MPS. This is expected to be \code{NULL} if the method fails, and the \code{NULL} is tested for in \code{\link{pars2x}}.}
}
\author{ W.H. Asquith}
\seealso{\code{\link{pars2x}} }
\examples{
\dontrun{
# Simulate from GLO and refit it. Occasionally, the simulated data
# will result in MLE or MPS failing to converge, just a note to users.
set.seed(3237)
x <- rlmomco(126, vec2par(c(2.5, 0.7, 0.3), type="glo"))
three.para.est <- x2pars(x, type="glo")
print(three.para.est$lmr$para) # 2.5598083 0.6282518 0.1819538
print(three.para.est$mle$para) # 2.5887340 0.6340132 0.2424734
print(three.para.est$mps$para) # 2.5843058 0.6501916 0.2364034}
}
\keyword{parameter (L-moments to)}
\keyword{parameter (maximum likelihood to)}
\keyword{parameter (maximum product spacing to)}
\keyword{parameter (estimation)}
