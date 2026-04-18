\name{parpe3}
\alias{parpe3}
\title{Estimate the Parameters of the Pearson Type III Distribution}
\description{
This function estimates the parameters of the Pearson Type III distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The L-moments in terms of the parameters are complicated and solved numerically. For the implementation in \pkg{lmomco}, the three parameters are \eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma} for the mean, standard deviation, and skew, respectively.

\bold{Remark} In spring 2026, a study was conducted using \code{inst/}\code{newcoes/}\code{makePE3coes.R} of \pkg{lmomco} sources for computation of better coefficients used in approximation of the relation between \eqn{\gamma} and L-skew (\eqn{\tau_3}). A vector of \eqn{\gamma} values was made for \eqn{\mu=0} and \eqn{\sigma = 0}, without a loss in generality, and the theoretical \eqn{\tau_3} computed by \code{\link{theoLmoms}} and \code{optim()} for fitting. These new coefficients relative to those taken from Hosking (1996) were implemented in spring 2026 by new settings in the \code{parpe3} source code. Slightly more accurate conversion of \eqn{\tau_3} to \eqn{\gamma} was attained. The structural formulation by Hosking was retained and only the fitting coefficients altered.
}
\usage{
parpe3(lmom, checklmom=TRUE, useHosking=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{useHosking}{A logical to use the old Hosking (1996) FORTRAN coefficients or use the new optimized coefficients of Spring 2026 custom to the \pkg{lmomco} package.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{pe3}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parpe3}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmompe3}}, \code{\link{cdfpe3}}, \code{\link{pdfpe3}}, \code{\link{quape3}}
}
\examples{
lmr <- lmoms(rnorm(20))
parpe3(lmr)

\dontrun{
  Tau3 <- runif(500000, min=-1, max=1)
  system.time(H <- sapply(Tau3, function(t) {
              parpe3(vec2lmom(c(0,1,t)), useHosking=TRUE )$para[3] }))
  #   user  system elapsed
  #  7.005   0.016   7.080
  system.time(W <- sapply(Tau3, function(t) {
              parpe3(vec2lmom(c(0,1,t)), useHosking=FALSE)$para[3] }))
  #   user  system elapsed
  #  7.084   0.023   7.186
  err <- H - W # error between the two methods
  plot(qnorm(pp(err)), abs(sort(err)), type="l", log="y",
       xlab="Standard normal variate of gamma parameter error",
       ylab="| Hosking 1996 minus lmomco 2026 |") #}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Pearson Type III}
