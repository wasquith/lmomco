\name{pargev}
\alias{pargev}
\title{Estimate the Parameters of the Generalized Extreme Value Distribution}
\description{
This function estimates the parameters of the Generalized Extreme Value distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomgev}}.
}
\usage{
pargev(lmom, checklmom=TRUE, useHosking=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{useHosking}{A logical to use the old Hosking (1996) FORTRAN coefficients or use direct computation by \code{uniroot()} wrapping the \eqn{\tau_3} and \eqn{k} relation shown under \code{\link{lmomgev}}.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gev}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{pargev}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomgev}},
\code{\link{cdfgev}}, \code{\link{pdfgev}}, \code{\link{quagev}}
}
\examples{
lmr <- lmoms(rnorm(20))
pargev(lmr)

\dontrun{
  set.seed(1)
  Tau3 <- runif(500000, min=-1, max=1)
  system.time(H <- sapply(Tau3, function(t) {
              pargev(vec2lmom(c(0,1,t)), useHosking=TRUE )$para[3] }))
  #   user  system elapsed
  #  8.036   0.033   8.198
  system.time(W <- sapply(Tau3, function(t) {
              pargev(vec2lmom(c(0,1,t)), useHosking=FALSE)$para[3] }))
  #   user  system elapsed
  # 19.348   0.062  19.737
  # By useHosking=FALSE is about 2.4 times slower but still fast.
  max(abs(H-W)) # 3.131861e-07 # Viability of either algorithm is proven.}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Generalized Extreme Value}
