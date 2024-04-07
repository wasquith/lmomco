\name{lmomcau}
\alias{lmomcau}
\title{Trimmed L-moments of the Cauchy Distribution}
\description{
This function estimates the trimmed L-moments of the Cauchy distribution given the parameters (\eqn{\xi} and \eqn{\alpha}) from \code{\link{parcau}}. The trimmed L-moments in terms of the parameters are \eqn{\lambda^{(1)}_1 = \xi},
\eqn{\lambda^{(1)}_2 = 0.69782723\alpha}, \eqn{\tau^{(1)}_{3, 5, \cdots} = 0}, \eqn{\tau^{(1)}_4 = 0.34280842}, and \eqn{\tau^{(1)}_6 = 0.20274358}. These TL-moments (trim=1) are symmetrical for the first L-moments defined because \eqn{\mathrm{E}[X_{1:n}]} and \eqn{\mathrm{E}[X_{n:n}]} undefined expectations for the Cauchy.
}
\usage{
lmomcau(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the trimmed L-moments. First element is
\eqn{\lambda^{(1)}_1}, second element is \eqn{\lambda^{(1)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau^{(1)}}, third element is \eqn{\tau^{(1)}_3} and so on. }
  \item{trim}{Level of symmetrical trimming used in the computation, which is unity.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is unity.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is unity.}
  \item{source}{An attribute identifying the computational source
of the L-moments: \dQuote{lmomcau}.}
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\seealso{\code{\link{parcau}}, \code{\link{cdfcau}}, \code{\link{pdfcau}}, \code{\link{quacau}} }
\examples{
X1 <- rcauchy(20)
lmomcau( parcau( TLmoms(X1, trim=1) ) )

alpha <- 30
tlmr <- theoTLmoms(vec2par(c(100, alpha), type="cau"), nmom=6, trim=1)
print( c(tlmr$lambdas[2] / alpha, tlmr$ratios[c(4,6)]), 8 )
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Cauchy}
