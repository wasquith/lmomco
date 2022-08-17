\name{Lcomoment.matrix}
\alias{Lcomoment.matrix}
\title{ Compute Sample L-comoment Matrix }
\description{
Compute the L-comoments from a rectangular \code{data.frame} containing arrays of random variables. The order of the L-comoments is specified.
}
\usage{
Lcomoment.matrix(DATAFRAME, k=1)
}
\arguments{
  \item{DATAFRAME}{A convential \code{data.frame} that is rectangular.}
  \item{k}{The order of the L-comoments to compute. Default is \eqn{k = 1}.}
}
\details{
L-comoments are computed for each item in the \code{data.frame}. L-comoments of order \eqn{k = 1} are means and co-means. L-comoments of order \eqn{k = 2} are L-scale and L-coscale values. L-comoments of order \eqn{k = 3} are L-skew and L-coskews. L-comoments of order \eqn{k = 4} are L-kurtosis and L-cokurtosis, and so on. The usual univariate L-moments of order \eqn{k} as seen from \code{\link{lmom.ub}} or \code{\link{lmoms}} are along the diagonal. This function does not make use of \code{\link{lmom.ub}} or \code{\link{lmoms}}. The function \code{\link{Lcomoment.matrix}} calls \code{\link{Lcomoment.Lk12}} for each cell in the matrix. The L-comoment matrix for \eqn{d}-random variables is
  \deqn{ \mathbf{\Lambda}_k = (\hat{\lambda}_{k[ij]})
       }
computed over the pairs (\eqn{X^{(i)},X^{(j)}}) where \eqn{1 \le i \le j \le d}.
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of L-comoment representation in the matrix: \dQuote{Lcomoments}.}
  \item{order}{The order of the matrix---specified by k in the argument list.}
  \item{matrix}{A kth order L-comoment matrix.}
}
\note{
The function begins with a capital letter. This is intentionally done so that lower case namespace is preserved. By using a capital letter now, then \code{lcomoment.matrix} remains an available name in future releases.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Serfling, R., and Xiao, P., 2007, A contribution to multivariate L-moments---L-comoment matrices: Journal of Multivariate Analysis, v. 98, pp. 1765--1781.
}
\seealso{\code{\link{Lcomoment.Lk12}}, \code{\link{Lcomoment.coefficients}}}
\examples{
D  <- data.frame(X1=rnorm(30), X2=rnorm(30), X3=rnorm(30))
L1 <- Lcomoment.matrix(D, k=1)
L2 <- Lcomoment.matrix(D, k=2)
}
\keyword{multivariate}
\keyword{L-comoment}
