\name{Lcomoment.correlation}
\alias{Lcomoment.correlation}
\title{ L-correlation Matrix (L-correlation through Sample L-comoments) }
\description{
Compute the L-correlation from an L-comoment matrix of order \eqn{k = 2}. This function assumes that the 2nd order matrix is already computed by the function \code{\link{Lcomoment.matrix}}.
}
\usage{
Lcomoment.correlation(L2)
}
\arguments{
  \item{L2}{A \eqn{k = 2} L-comoment matrix from \code{Lcomoment.matrix(Dataframe,k=2)}.}
}
\details{
L-correlation is computed by \code{Lcomoment.coefficients(L2,L2)} where \code{L2} is second order L-comoment matrix. The usual L-scale values as seen from \code{\link{lmom.ub}} or \code{\link{lmoms}} are along the diagonal. This function does not make use of \code{\link{lmom.ub}} or \code{\link{lmoms}} and can be used to verify  computation of \eqn{\tau} (coefficient of L-variation).
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of L-comoment representation in the matrix: \dQuote{Lcomoment.coefficients}.}
  \item{order}{The order of the matrix---extracted from the first matrix in arguments.}
  \item{matrix}{A \eqn{k \ge 2} L-comoment coefficient matrix.}
}
\note{
The function begins with a capital letter. This is intentionally done so that lower case namespace is preserved. By using a capital letter now, then \code{lcomoment.correlation} remains an available name in future releases.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Serfling, R., and Xiao, P., 2007, A contribution to multivariate L-moments---L-comoment matrices: Journal of Multivariate Analysis, v. 98, pp. 1765--1781.
}
\seealso{\code{\link{Lcomoment.matrix}}, \code{\link{Lcomoment.correlation}}}
\examples{
D   <- data.frame(X1=rnorm(30), X2=rnorm(30), X3=rnorm(30))
L2  <- Lcomoment.matrix(D,k=2)
RHO <- Lcomoment.correlation(L2)
\dontrun{
"SerfXiao.eq17" <-
 function(n=25, A=10, B=2, k=4,
          method=c("pearson","lcorr"), wrt=c("12", "21")) {
   method <- match.arg(method); wrt <- match.arg(wrt)
   # X1 is a linear regression on X2
   X2 <- rnorm(n); X1 <- A + B*X2 + rnorm(n)
   r12p <- cor(X1,X2) # Pearson's product moment correlation
   XX <- data.frame(X1=X1, X2=X2) # for the L-comoments
   T2 <- Lcomoment.correlation(Lcomoment.matrix(XX, k=2))$matrix
   LAMk <- Lcomoment.matrix(XX, k=k)$matrix # L-comoments of order k
   if(wrt == "12") { # is X2 the sorted variable?
      lmr <- lmoms(X1, nmom=k); Lamk <- LAMk[1,2]; Lcor <- T2[1,2]
   } else {          # no X1 is the sorted variable (21)
      lmr <- lmoms(X2, nmom=k); Lamk <- LAMk[2,1]; Lcor <- T2[2,1]
   }
   # Serfling and Xiao (2007, eq. 17) state that
   # L-comoment_k[12] = corr.coeff * Lmoment_k[1] or
   # L-comoment_k[21] = corr.coeff * Lmoment_k[2]
   # And with the X1, X2 setup above, Pearson corr. == L-corr.
   # There will be some numerical differences for any given sample.
   ifelse(method == "pearson",
             return(lmr$lambdas[k]*r12p - Lamk),
             return(lmr$lambdas[k]*Lcor - Lamk))
   # If the above returns a expected value near zero then, their eq.
   # is numerically shown to be correct and the estimators are unbiased.
}

# The means should be near zero.
nrep <- 2000; seed <- rnorm(1); set.seed(seed)
mean(replicate(n=nrep, SerfXiao.eq17(method="pearson", k=4)))
set.seed(seed)
mean(replicate(n=nrep, SerfXiao.eq17(method="lcorr", k=4)))
# The variances should nearly be equal.
seed <- rnorm(1); set.seed(seed)
var(replicate(n=nrep, SerfXiao.eq17(method="pearson", k=6)))
set.seed(seed)
var(replicate(n=nrep, SerfXiao.eq17(method="lcorr", k=6)))
}
}
\keyword{multivariate}
\keyword{L-comoment}
