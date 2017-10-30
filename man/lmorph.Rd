\name{lmorph}
\alias{lmorph}
\title{ Morph an L-moment Object }
\description{
Morph or change one L-moment object type into another. The first L-moment object created for \pkg{lmomco} used an \R \code{list} with named L-moment values (\code{\link{lmom.ub}}) such as \code{L1} or \code{TAU3}. This object was bounded for L-moment orders less than or equal to five. However, subsequent \pkg{lmomco} development in early 2006 that was related to the trimmed L-moments suggested that an alternative L-moment object structure be used that utilized two vectors for the L-moments and the L-moment ratios (\code{\link{lmorph}}). This second object type is not bounded by L-moment order. In turn it became important to seemlessly morph from one object structure to the other and back again. The canonical structure of the first L-moment object type is documented under \code{\link{lmom.ub}}; whereas, the canonical structure for the second L-moment object type is documented under \code{\link{lmoms}} (actually through \code{\link{TLmoms}}). Because the first L-moment object is bounded by five, L-moment order larger than this will be ignored in the morphing process.
}
\usage{
lmorph(lmom)
}
\arguments{
  \item{lmom}{An L-moment object of type like \code{lmom.ub} or \code{lmoms}.}
}
\value{
A two different \R \code{list}s (L-moment objects), which are the opposite of the argument type---see the documentation for \code{\link{lmom.ub}} and \code{\link{lmoms}}.
}
\note{
If any of the trimming characteristics of the second type of L-moment object (\code{trim}, \code{leftrim}, or \code{rightrim}) have a greater than zero value, then conversion to the L-moment object with named values will not be performed. A message will be provided that the conversion was not performed.  In April 2014, it was decided that all \code{lmomCCC()} functions, such as \code{\link{lmomgev}} or \code{\link{lmomnor}}, would be standardized to the less limited and easier to maintain vector output style of \code{\link{lmoms}}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmom.ub}}, \code{\link{lmoms}}, \code{\link{TLmoms}}}
\examples{
lmr <- lmom.ub(c(123,34,4,654,37,78))
lmorph(lmr)
lmorph(lmorph(lmr))
}
\keyword{L-moment (utility)}
