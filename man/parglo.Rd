\name{parglo}
\alias{parglo}
\title{Estimate the Parameters of the Generalized Logistic Distribution}
\description{
This function estimates the parameters of the Generalized Logistic distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomglo}}.
}
\usage{
parglo(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{glo}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parglo}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomglo}}, \code{\link{cdfglo}},  \code{\link{pdfglo}}, \code{\link{quaglo}}
}
\examples{
lmr <- lmoms(rnorm(20))
parglo(lmr)
\dontrun{
# A then Ph.D. student, L. Read inquired in February 2014 about the relation between
# GLO and the "Log-Logistic" distributions:
par.glo  <- vec2par(c(10, .56, 0), type="glo")         # Define GLO parameters
par.lnlo <- c(exp(par.glo$para[1]), 1/par.glo$para[2]) # Equivalent LN-LO parameters
F <- nonexceeds(); qF <- qnorm(F) # use a real probability axis to show features
plot(qF, exp(quaglo(F, par.glo)), type="l", lwd=5, xaxt="n", log="y",
     xlab="", ylab="QUANTILE") # notice the exp() wrapper on the GLO quantiles
lines(qF, par.lnlo[1]*(F/(1-F))^(1/par.lnlo[2]), col=2, lwd=2) # eq. for LN-LO
add.lmomco.axis(las=2, tcl=0.5, side.type="RI", otherside.type="NPP")
}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Generalized Logistic}
\keyword{Distribution: Log-Logistic}
