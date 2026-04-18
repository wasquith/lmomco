\name{pargno}
\alias{pargno}
\title{Estimate the Parameters of the Generalized Normal Distribution}
\description{
This function estimates the parameters of the Generalized Normal (Log-Normal3) distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomgno}}.

\bold{Remark} In spring 2026, a study was conducted using using \code{inst/}\code{newcoes/}\code{makeGNOcoes.R} of \pkg{lmomco} sources for computation of better coefficients used in approximation of the relation between \eqn{\kappa} and L-skew (\eqn{\tau_3}). A vector of \eqn{\kappa} values was made for \eqn{\mu=0} and \eqn{\lambda_2 = 0}, without a loss in generality, and the theoretical \eqn{\tau_3} computed by \code{\link{theoLmoms}} and \code{optim()} for fitting. These new coefficients relative to those taken from Hosking (1996) were implemented in spring 2026 by new settings in the \code{parpe3} source code. Slightly more accurate conversion of \eqn{\tau_3} to \eqn{\gamma} was attained. The structural formulation by Hosking was retained and only the fitting coefficients altered.
}
\usage{
pargno(lmom, checklmom=TRUE, useHosking=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{useHosking}{A logical to use the old Hosking (1996) FORTRAN coefficients or use the new optimized coefficients of Spring 2026 custom to the \pkg{lmomco} package.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gno}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{pargno}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomgno}}, \code{\link{cdfgno}}, \code{\link{pdfgno}}, \code{\link{quagno}}, \code{\link{parln3}}
}
\examples{
lmr <- lmoms(rnorm(20))
pargno(lmr)

\dontrun{
x   <- c(2.4, 2.7, 2.3, 2.5, 2.2, 62.4, 3.8, 3.1)
lmr <- lmoms(x, nmom=3)
gno <- pargno(lmr, useHosking=TRUE)  # triggers warning
# Hosking's limit is Tau3=+-0.95
gno <- pargno(lmr, useHosking=FALSE)
tmr <- theoLmoms(gno, nmom=3)
print( lmr$lambdas ) # matching
print( tmr$lambdas ) # matching
print( lmr$ratios  ) # matching
print( tmr$ratios  ) # matching }

\dontrun{
  Tau3 <- runif(500000, min=-1, max=1)
  system.time(H <- sapply(Tau3, function(t) {
              pargno(vec2lmom(c(0,1,t)), useHosking=TRUE )$para[3] }))
  #   user  system elapsed
  #  9.072   0.027   9.200
  system.time(W <- sapply(Tau3, function(t) {
              pargno(vec2lmom(c(0,1,t)), useHosking=FALSE)$para[3] }))
  #   user  system elapsed
  #  7.261   0.024   7.389
  # By useHosking=FALSE is about 1.25 times slower but still fast.
  plot(sort(H), sort(W), type="p")
  abline(0,1, col="red") # Hosking's Tau3=0.95 truncation evident.}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Generalized Normal}
