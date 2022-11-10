\name{vec2TLmom}
\alias{vec2TLmom}
\title{Convert a Vector of TL-moments to a TL-moment Object}
\description{
This function converts a vector of trimmed L-moments (TL-moments) to a TL-moment object of \pkg{lmomco} by dispatch to \code{\link{vec2lmom}}. The object is an \R \code{list}. This function is intended to facilitate the use of TL-moments that the user might have from other sources. The trimming on the left-tail is denoted by \eqn{t} and the trimming on the right-tail is denoted as \eqn{s}. The first five TL-moments are \eqn{\lambda^{(t,s)}_1}, \eqn{\lambda^{(t,s)}_2}, \eqn{\lambda^{(t,s)}_3}, \eqn{\lambda^{(t,s)}_4}, \eqn{\lambda^{(t,s)}_5}, \eqn{\tau^{(t,s)}}, \eqn{\tau^{(t,s)}_3}, \eqn{\tau^{(t,s)}_4}, and \eqn{\tau^{(t,s)}_5}. The function supports TL-moments and TL-moment ratios of arbitrary length. Because in typical practice the \eqn{k \ge 3} order L-moments are dimensionless ratios (\eqn{\tau^{(t,s)}_3}, \eqn{\tau^{(t,s)}_4}, and \eqn{\tau^{(t,s)}_5}), this function computes \eqn{\lambda^{(t,s)}_3}, \eqn{\lambda^{(t,s)}_4}, \eqn{\lambda^{(t,s)}_5} from \eqn{\lambda^{(t,s)}_2} and the ratios. However, typical practice is not set on the use of \eqn{\lambda^{(t,s)}_2} or \eqn{\tau^{(t,s)}} as measure of dispersion. Therefore, this function takes an \code{lscale} optional logical argument---if \eqn{\lambda^{(t,s)}_2} is provided and \code{lscale=TRUE}, then \eqn{\tau} is computed by the function and if \eqn{\tau} is provided, then \eqn{\lambda^{(t,s)}_2} is computed by the function. The trim level of the TL-moment is required. Lastly, it might be common for \eqn{t=s} and hence symmetrical trimming is used.
}
\usage{
vec2TLmom(vec, ...)
}
\arguments{
  \item{vec}{A vector of L-moment values in \eqn{\lambda^{(t,s)}_1}, \eqn{\lambda^{(t,s)}_2} or \eqn{\tau^{(t,s)}}, \eqn{\tau^{(t,s)}_3}, \eqn{\tau^{(t,s)}_4}, and \eqn{\tau^{(t,s)}_5} order.}
  \item{...}{The arguments used by \code{\link{vec2lmom}}.}
}
\value{
  An \R \code{list} is returned where \eqn{t} represents the \code{trim} level.

  \item{lambdas}{Vector of the TL-moments. First element is \eqn{\lambda^{(t,s)}_1}, second element is \eqn{\lambda^{(t,s)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is \eqn{\tau^{(t,s)}}, third element is \eqn{\tau^{(t,s)}_3} and so on.}
  \item{trim}{Level of symmetrical trimming, which should equal \code{NULL} if asymmetrical trimming is used.}
  \item{leftrim}{Level of trimming of the left-tail of the sample.}
  \item{rightrim}{Level of trimming of the right-tail of the sample.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{TLmoms}.}
}
\note{
The motiviation for this function that arrange trivial arguments for \code{\link{vec2lmom}} is that it is uncertain how TL-moments will grow in the research community and there might someday be a needed for alternative support without having to touch \code{\link{vec2lmom}}. Plus there is nice function name parallelism in having a dedicated function for the TL-moments as there is for L-moments and probability-weighted moments.
}
\author{W.H. Asquith}
\seealso{\code{\link{TLmoms}}, \code{\link{vec2lmom}}}
\examples{
TL <- vec2TLmom(c(12,0.6,0.34,0.20,0.05),lscale=FALSE,trim=1)
}
\keyword{utility (vector)}
