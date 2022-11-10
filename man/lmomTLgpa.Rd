\name{lmomTLgpa}
\alias{lmomTLgpa}
\title{Trimmed L-moments of the Generalized Pareto Distribution}
\description{
This function estimates the symmetrical trimmed L-moments (TL-moments) for \eqn{t=1} of the Generalized Pareto distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) from \code{\link{parTLgpa}}.
The TL-moments in terms of the parameters are
\deqn{\lambda^{(1)}_1 = \xi + \frac{\alpha(\kappa+5)}{(\kappa+3)(\kappa+2)} \mbox{,}}
\deqn{\lambda^{(1)}_2 = \frac{6\alpha}{(\kappa+4)(\kappa+3)(\kappa+2)} \mbox{,}}
\deqn{\tau^{(1)}_3 = \frac{10(1-\kappa)}{9(\kappa+5)} \mbox{, and}}
\deqn{\tau^{(1)}_4 = \frac{5(\kappa-1)(\kappa-2)}{4(\kappa+6)(\kappa+5)}  \mbox{.}}
}
\usage{
lmomTLgpa(para)
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
  \item{source}{An attribute identifying the computational source of the TL-moments: \dQuote{lmomTLgpa}.}
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomgpa}}, \code{\link{parTLgpa}}, \code{\link{cdfgpa}}, \code{\link{pdfgpa}}, \code{\link{quagpa}} }
\examples{
TL <- TLmoms(c(123,34,4,654,37,78,21,3400),trim=1)
TL
lmomTLgpa(parTLgpa(TL))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Generalized Pareto}
