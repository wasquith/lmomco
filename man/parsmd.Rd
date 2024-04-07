\name{parsmd}
\alias{parsmd}
\title{Estimate the Parameters of the Singh--Maddala Distribution}
\description{
This function estimates the parameters of the Singh--Maddala (Burr Type XII) distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The L-moments in terms of the parameters are complicated and solved numerically. Extensive study of the computational limits of the \R implementation are incorporated within the source code of the function. The file \code{lmomco/inst/doc/domain_of_smd.R} contains the algorithmic sweep used to compute the L-skew and L-kurtosis attainable domain of the distribution.
}
\usage{
parsmd(lmom, checklmom=TRUE, checkbounds=TRUE, snap.tau4=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_3} and \eqn{\tau_4} inequality, \cr \code{\link{are.lmom.valid}}). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{checkbounds}{Should the lower bounds of \eqn{\tau_4} be verified and if sample \eqn{\hat\tau_3} and \eqn{\hat\tau_4} are outside of these bounds, then \code{NA} are returned for the solutions.}
  \item{snap.tau4}{A logical to trigger the application of the empirical limits of the distribution in terms of \eqn{\tau_3} and \eqn{\tau_4} wherein parameter estimation appears numerically possible and such parameters return the given values of these L-moment ratios. The lower and upper limits of \eqn{\tau_4} are defined by separate polynomials as functions of \eqn{\tau_3}. If the logical is true, then \eqn{\tau_4} in excess of the upper bounds are assigned to the upper bounds and \eqn{\tau_4} in deficit of the lower bounds are assigned to the lower bounds. Messages within the returned parameter object are provided if the snapping occurs.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{smd}.}
  \item{para}{The parameters of the distribution.}
  \item{last_para}{The last or final iteration of the parameters that are the same as \code{para} if \code{ifail} is zero. This provides a way to preserve where the parameter left off or gave up.}
  \item{source}{The source of the parameters: \dQuote{parsmd}.}
  \item{iter}{The number of iteration attempts looping on the \code{optim()} call.}
  \item{rt}{The output of the \code{optim()} call.}
  \item{message}{A message from \code{parsmd}, which generally involves \code{checkbounds=TRUE} and \cr\code{snap.tau4=TRUE} on the resetting or snapping of the \eqn{\tau_3} and \eqn{\tau_4} to the computational bounds for the distribution.}
  \item{ifail}{A interger flag to status of the operations: -1 means that the L-moments are invalid (if they are checked), 0 means that the parameter estimation appears successful, and 1 means that the parameter estimation appears to have failed.}
}
\references{
Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomsmd}}, \code{\link{cdfsmd}}, \code{\link{pdfsmd}}, \code{\link{quasmd}}
}
\examples{
lmr <- lmoms(rnorm(20))
parsmd(lmr, snap.tau4=TRUE)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
