\name{pwmLC}
\alias{pwmLC}
\title{Sample Probability-Weighted Moments for Left-Tail Censoring }
\description{
Compute the sample probability-weighted moments (PWMs) for left-tail censored data set---that is a data set censored from below. The censoring threshold is denoted as \eqn{T}.
}
\usage{
pwmLC(x, threshold=NULL, nmom=5, sort=TRUE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{threshold}{The left-tail censoring (lower) threshold.}
  \item{nmom}{Number of PWMs to return.}
  \item{sort}{Do the data need sorting? Note that convention is the have a \eqn{\beta'_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
}
\value{
  An \R \code{list} is returned.

  \item{Aprimebetas}{The A'-type PWMs. These should be same as \code{pwm()} returns if there is no censoring. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{Bprimebetas}{The B'-type PWMs. These should be \code{NA} if there is no censoring. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i = 1} of the \code{betas} vector.}
  \item{source}{Source of the PWMs: \dQuote{pwmLC}.}
  \item{threshold}{The upper censoring threshold.}
  \item{zeta}{The left censoring fraction: \code{numbelowthreshold/samplesize}.}
  \item{numbelowthreshold}{Number of data points equal to or above the threshold.}
  \item{observedsize}{Number of real data points in the sample (above the threshold).}
  \item{samplesize}{Number of actual sample values.}
}
\details{
There is some ambiguity if the threshold also numerically equals valid data in the data set. In the data for the examples below, which are taken from elsewhere, there are real observations at the censoring level. One can see how a hack is made to marginally decrease or increase the data or the threshold for the computations. This is needed because the code uses
\preformatted{
sapply(x, function(v) { if(v >= T) return(T); return(v) } )
}
to reset the data vector \code{x}. By operating on the data in this fashion one can toy with various levels of the threshold for experimental purposes; this seemed a more natural way for general implementation. The code sets \eqn{n} = \code{length(x)} and \eqn{m} = \code{n - length(x[x == T])}, which also seems natural. The \eqn{\beta^A_r} are computed by dispatching to \code{\link{pwm}}.
}
\references{
Zafirakou-Koulouris, A., Vogel, R.M., Craig, S.M., and Habermeier, J., 1998, L-moment diagrams for censored observations: Water Resources Research, v. 34, no. 5, pp. 1241--1249.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmoms}}, \code{\link{pwm2lmom}}, \code{\link{pwm}}, \code{\link{pwmRC}} }
\examples{
#
}
\keyword{probability-weighted moment (sample)}
\keyword{data censoring}
