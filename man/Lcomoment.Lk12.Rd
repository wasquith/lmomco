\name{Lcomoment.Lk12}
\alias{Lcomoment.Lk12}
\title{Compute a Single Sample L-comoment }
\description{
Compute the L-comoment (\eqn{\lambda_{k[12]}}) for a given pair of sample of \eqn{n} random variates \eqn{\{(X_i^{(1)}, X_i^{(1)}), 1 \le i \le n \}} from a joint distribution \eqn{H(x^{(1)}, x^{(2)})} with marginal distribution functions \eqn{F_1} and \eqn{F_2}. When the \eqn{X^{(2)}} are sorted to form the sample order statistics \eqn{X^{(2)}_{1:n} \le X^{(2)}_{2:n} \le \cdots \le X^{(2)}_{n:n}}, then the element of \eqn{X^{(1)}} of the unordered (at leasted expected to be) but shuffled set \eqn{\{X^{(1)}_1, \ldots, X^{(1)}_n\}} that is paired with \eqn{X^{(2)}_{r:n}} the \emph{concomitant} \eqn{X^{(12)}_{[r:n]}} of \eqn{X^{(2)}_{r:n}}. (The shuffling occurs by the sorting of \eqn{X^{(2)}}.) The \eqn{k \ge 1}-order L-comoments are defined (Serfling and Xiao, 2007, eq. 26) as
\deqn{\hat\lambda_{k[12]} = \frac{1}{n}\sum_{r=1}^n w^{(k)}_{r:n} X^{(12)}_{[r:n]}\mbox{,}}
where \eqn{w^{(k)}_{r:n}} is defined under \code{\link{Lcomoment.Wk}}. (The author is aware that \eqn{k \ge 1} is \eqn{k \ge 2} in Serfling and Xiao (2007) but \eqn{k=1} returns sample means. This matters only in that the \pkg{lmomco} package returns matrices for \eqn{k \ge 1} by \code{\link{Lcomoment.matrix}} even though the off diagnonals are \code{NAs}.)
}
\usage{
Lcomoment.Lk12(X1,X2,k=1)
}
\arguments{
  \item{X1}{A vector of random variables (a sample of random variable 1).}
  \item{X2}{Another vector of random variables (a sample of random variable 2).}
  \item{k}{The order of the L-comoment to compute. The default is 1.}
}
\details{
Now directing explanation of L-comoments with some reference heading into \R code. L-comoments of random variable \code{X1} (a vector) are computed from the concomitants of \code{X2} (another vector). That is, \eqn{X2} is sorted in ascending order to create the order statistics of \code{X2}. During the sorting process, \code{X1} is reshuffled to the order of \code{X2} to form the concomitants of \code{X2} (denoted as \code{X12}). So the trailing \code{2} is the sorted variable and the leading \code{1} is the variable that is shuffled. The \code{X12} in turn are used in a weighted summation and expectation calculation to compute the L-comoment of \code{X1} with respect to \code{X2} such as by \code{Lk3.12 <-} \code{Lcomoment.Lk12(X1,X2,k=3)}. The notation of \code{Lk12} is to read \dQuote{Lambda for kth order L-comoment}, where the \code{12} portion of the notation reflects that of Serfling and Xiao (2007) and then Asquith (2011). The weights for the computation are derived from calls made by \code{\link{Lcomoment.Lk12}} to the weight function \code{\link{Lcomoment.Wk}}. The L-comoments of \code{X2} are computed from the concomitants of \code{X1}, and the \code{X21} are formed by sorting \code{X1} in ascending order and in turn shuffling \code{X2} by the order of \code{X1}. The often asymmetrical L-comoment of \code{X2} with respect to \code{X1} is readily done (\code{Lk3.21 <-} \code{Lcomoment.Lk12(X2,X1,k=3)}) and is not necessarily equal to (\code{Lk3.12 <-} \code{Lcomoment.Lk12(X1,X2,k=3)}).
}
\value{A single L-comoment.}
\note{
The function begins with a capital letter. This is intentionally done so that lower case namespace is preserved. By using a capital letter now, then \code{lcomoment.Lk12} or similar remains an available name in future releases.
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Serfling, R., and Xiao, P., 2007, A contribution to multivariate L-moments---L-comoment matrices: Journal of Multivariate Analysis, v. 98, pp. 1765--1781.
}
\author{W.H. Asquith}
\seealso{\code{\link{Lcomoment.matrix}}, \code{\link{Lcomoment.Wk}}}
\examples{
X1 <- rnorm(101); X2 <- rnorm(101) + X1
Lcoskew12 <- Lcomoment.Lk12(X1,X2, k=3)
Lcorr12 <- Lcomoment.Lk12(X1,X2,k=2)/Lcomoment.Lk12(X1,X1,k=2)
rhop12 <- cor(X1, X2, method="pearson")
print(Lcorr12 - rhop12) # smallish number
}
\keyword{multivariate}
\keyword{L-comoment}
