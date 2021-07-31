\name{tlmr2par}
\alias{tlmr2par}
\title{Sample Trimmed L-moments to Fitted Distribution}
\description{
Parameter estimation of a distribution given initial estimate of the parameters of the distribution to the sample trimmed L-moment (TL-moment) using numerical optimization.
}
\usage{
tlmr2par(x, type, para.int=NULL, trim=NULL, leftrim=NULL, rightrim=NULL, ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{type}{Three character (minimum) distribution type (for example, \code{type="gev"}, see \code{\link{dist.list}}.}
  \item{para.int}{Initial parameters as a vector \eqn{\Theta} or as an \pkg{lmomco} parameter \dQuote{object} from say \code{\link{vec2par}}. If a vector is given, then internally \code{\link{vec2par}} is called with distribution equal to \code{type}.}
  \item{trim}{Level of symmetrical trimming to use in the computations. Although \code{NULL} is in the argument list, the default is 0---the usual L-moment is returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
  \item{rightrim}{Level of trimming of the right-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
  \item{...}{Other arguments to pass to the \code{optim()} function.}
}
\value{
  An \R \code{list} is returned.  This list should contain at least the following items, but some distributions such as the \code{revgum} have extra.

  \item{type}{The type of distribution in three character (minimum) format.}
  \item{para}{The parameters of the distribution.}
  \item{source}{Attribute specifying source of the parameters.}
  \item{rt}{The list from the \code{optim()} function.}
  \item{para.int}{A copy of the initial parameters given.}
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{ W.H. Asquith}
\seealso{\code{\link{theoTLmoms}}, \code{\link{TLmoms}}}
\examples{
FF <- lmomco::nonexceeds(); qFF <- qnorm(FF); type <- "gev"
the.data <- c(3.064458, 3.139879, 3.167317, 3.225309, 3.324282, 3.330414,
             3.3304140, 3.340444, 3.357935, 3.376577, 3.378398, 3.392697,
             3.4149730, 3.421604, 3.424882, 3.434569, 3.448706, 3.451786,
             3.4517860, 3.462398, 3.465383, 3.469822, 3.491362, 3.501059,
             3.5224440, 3.523746, 3.527630, 3.527630, 3.531479, 3.546543,
             3.5932860, 3.597695, 3.600973, 3.614897, 3.620136, 3.660865,
             3.6848450, 3.820858, 4.708421)
# visually, looks like 3 outliers to the left and one outlier to the right
whole.para <- lmom2par(lmoms(the.data), type=type)
trim.para  <- tlmr2par(the.data, type, para.int=whole.para, leftrim=3, rightrim=1)

plot( qFF, qlmomco(FF, whole.para), type="l", lwd=2, ylim=c(3.1,4.8),
           xlab="Standard normal variate", ylab="Some phenomena, log10(cfs)")
lines(qFF, qlmomco(FF, trim.para), col=4, lwd=3)
points(qnorm(pp(the.data)), sort(the.data), pch=1, col=grey(.5))
legend("topleft", c("L-moments", "TL-moments(3,2)"), bty="n",
                  lty=c(1,1), lwd=c(2,3), col=c(1,4))
# see the massive change from the whole sample to the trim(3,1) curve
}
\keyword{L-moment (sample)}
\keyword{L-moment (trimmed sample)}
