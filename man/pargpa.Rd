\name{pargpa}
\alias{pargpa}
\title{Estimate the Parameters of the Generalized Pareto Distribution}
\description{
This function estimates the parameters of the Generalized Pareto distribution given the L-moments of the data in an ordinary L-moment object (\code{\link{lmoms}}) or a trimmed L-moment object (\code{\link{TLmoms}} for \code{t=1}). The relations between distribution parameters and L-moments are seen under \code{\link{lmomgpa}} or \code{\link{lmomTLgpa}}.
}
\usage{
pargpa(lmom, zeta=1, xi=NULL, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}}, \code{\link{TLmoms}} with \code{trim=0}, or \code{\link{vec2lmom}}.}
  \item{zeta}{The right censoring fraction. If less than unity then a dispatch to the \code{\link{pargpaRC}} is made and the \code{lmom} argument must contain the B-type L-moments. If the data are not right censored, then this value must be left alone to the default of unity.}
  \item{xi}{The lower limit of the distribution. If \eqn{\xi} is known, then alternative algorithms are used.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gpa}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{pargpa}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomgpa}},
\code{\link{cdfgpa}}, \code{\link{pdfgpa}}, \code{\link{quagpa}}
}
\examples{
X   <- rexp(200)
lmr <- lmoms(X)
P1  <- pargpa(lmr)
P2  <- pargpa(lmr, xi=0.25)

\dontrun{
F <- nonexceeds()
plot(pp(X), sort(X))
lines(F,quagpa(F,P1))  # black line
lines(F, quagpa(F,P2), col=2)   # red line
}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Generalized Pareto}
