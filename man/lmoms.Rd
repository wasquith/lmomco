\name{lmoms}
\alias{lmoms}
\title{ The Sample L-moments and L-moment Ratios }
\description{
Compute the sample L-moments. The mathematical expression for sample L-moment computation is shown under \code{\link{TLmoms}}. The formula jointly handles sample L-moment computation and sample TL-moment (Elamir and Seheult, 2003) computation. A description of the most common L-moments is provided under \code{\link{lmom.ub}}.
}
\usage{
lmoms(x, nmom=5, no.stop=FALSE, vecit=FALSE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{no.stop}{A logical to return \code{NULL} instead of issuing a \code{stop()} if \code{nmom} is greater than the sample size or if all the values are equal. This is a very late change (decade+) to the foundational function in the package. Auxiliary coding to above this function to avoid the internal \code{stop()} became non-ignorable in large data mining exercises. It was a design mistake to have the \code{stop()} and not a \code{warning()} instead.}
  \item{vecit}{A logical to return the first two \eqn{\lambda_i \in 1,2} and then the \eqn{\tau_i \in 3,\cdots} where the length of the returned vector is controlled by the \code{nmom} argument. This argument will store the trims (see \code{\link{TLmoms}}) as \code{NULL} used (see the \bold{Example} that follows).}
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
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmoms}.}
}
\note{
This function computes the L-moments through the generalization of the TL-moments (\code{\link{TLmoms}}). In fact, this function calls the default TL-moments with no trimming of the sample. This function is equivalent to \code{\link{lmom.ub}}, but returns a different data structure. The \code{\link{lmoms}} function is preferred by the author.
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational statistics and data analysis, vol. 43, pp. 299-314.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{ W.H. Asquith}

\seealso{\code{\link{lmom.ub}}, \code{\link{TLmoms}}, \code{\link{lmorph}}, \code{\link{lmoms.bernstein}}, \code{\link{vec2lmom}}}
\examples{
lmoms(rnorm(30),nmom=4)

vec2lmom(lmoms(rexp(30), nmom=3, vecit=TRUE)) # re-vector
}
\keyword{L-moment (sample)}
