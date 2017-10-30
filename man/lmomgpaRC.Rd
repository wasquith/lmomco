\name{lmomgpaRC}
\alias{lmomgpaRC}
\title{B-type L-moments of the Generalized Pareto Distribution with Right-Tail Censoring}
\description{
This function computes the \dQuote{B}-type L-moments of the Generalized Pareto distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) from \code{\link{pargpaRC}} and the right-tail censoring fraction \eqn{\zeta}. The B-type L-moments in terms of the parameters are
\deqn{\lambda^B_1 = \xi + \alpha m_1 \mbox{,}}
\deqn{\lambda^B_2 = \alpha (m_1 - m_2) \mbox{,}}
\deqn{\lambda^B_3 = \alpha (m_1 - 3m_2 + 2m_3)\mbox{,}}
\deqn{\lambda^B_4 = \alpha (m_1 - 6m_2 + 10m_3 - 5m_4)\mbox{, and}}
\deqn{\lambda^B_5 = \alpha (m_1 - 10m_2 + 30m_3 - 35m_4 + 14m_5)\mbox{,}}
where \eqn{m_r = \lbrace 1-(1-\zeta)^{r+\kappa}\rbrace/(r+\kappa)} and \eqn{\zeta} is the right-tail censor fraction or the probability \eqn{\mathrm{Pr}\lbrace \rbrace} that \eqn{x} is less than the quantile at \eqn{\zeta} nonexceedance probability: (\eqn{\mathrm{Pr}\lbrace x < X(\zeta) \rbrace}). In other words, if \eqn{\zeta = 1}, then there is no right-tail censoring. Finally, the \code{RC} in the function name is to denote \code{R}ight-tail \code{C}ensoring.
}
\usage{
lmomgpaRC(para)
}
\arguments{
  \item{para}{The parameters of the distribution. Note that if the \eqn{\zeta} part of the parameters (see \code{\link{pargpaRC}}) is not present then \code{zeta=1} (no right-tail censoring) is assumed.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomgpaRC}.}
  \item{message}{For clarity, this function adds the unusual message to an L-moment object that the \code{lambdas} and \code{ratios} are B-type L-moments.}
  \item{zeta}{The censoring fraction. Assumed equal to unity if not present in the \code{gpa} parameter object.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data, in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan, chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pargpa}}, \code{\link{pargpaRC}}, \code{\link{lmomgpa}}, \code{\link{cdfgpa}}, \code{\link{pdfgpa}}, \code{\link{quagpa}} }
\examples{
para <- vec2par(c(1500,160,.3),type="gpa") # build a GPA parameter set
lmorph(lmomgpa(para))
lmomgpaRC(para) # zeta = 1 is internally assumed if not available
# The previous two commands should output the same parameter values from
# independent code bases.
# Now assume that we have the sample parameters, but the zeta is nonunity.
para$zeta = .8
lmomgpaRC(para) # The B-type L-moments.
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Generalized Pareto}
\keyword{data censoring}
