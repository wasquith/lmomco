\name{lmomRCmark}
\alias{lmomRCmark}
\title{ Sample L-moment for Right-Tail Censoring by a Marking Variable }
\description{
Compute the sample L-moments for right-tail censored data set in which censored data values are identified by a marking variable.
}
\usage{
lmomRCmark(x, rcmark=NULL, r=1, sort=TRUE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{rcmark}{The right-tail censoring (upper) marking variable for unknown threshold: 1 is uncensored, 0 is censored.}
  \item{r}{The L-moment order to return, default is the mean.}
  \item{sort}{Do the data need sorting? The availability of this option is to avoid unnecessary overhead of sorting on each call to this function by the primary higher-level function \code{\link{lmomsRCmark}}.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\hat{\lambda}^{(0,0)}_1}, second element is \eqn{\hat{\lambda}^{(0,0)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\hat{\tau}^{(0,0)}}, third element is \eqn{\hat{\tau}^{(0,0)}_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal \code{NULL} if asymmetrical trimming was used.}
  \item{leftrim}{Level of left-tail trimming used in the computation.}
  \item{rightrim}{Level of right-tail trimming used in the computation.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomsRCmark}.}
}
\references{
Wang, Dongliang, Hutson, A.D., Miecznikowski, J.C., 2010, L-moment estimation for parametric survival models given censored data: Statistical Methodology, v. 7, no. 6, pp. 655--667.
}
\author{W.H. Asquith}
\seealso{ \code{\link{lmomsRCmark}} }
\examples{
# See example under lmomsRCmark
}
\keyword{L-moment (sample)}
\keyword{data censoring}
