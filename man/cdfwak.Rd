\name{cdfwak}
\alias{cdfwak}
\title{Cumulative Distribution Function of the Wakeby Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Wakeby distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{\delta}) computed by \code{\link{parwak}}. The cumulative distribution function has no explicit form, but the \code{\link{pdfwak}} (density) and \code{\link{quawak}} (quantiles) do.
}
\usage{
cdfwak(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parwak}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfwak}}, \code{\link{quawak}}, \code{\link{lmomwak}}, \code{\link{parwak}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfwak(50,parwak(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Wakeby}
