\name{pargov}
\alias{pargov}
\title{Estimate the Parameters of the Govindarajulu Distribution}
\description{
This function estimates the parameters of the Govindarajulu distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments also are seen under \code{\link{lmomgov}}. The \eqn{\beta} is estimated as
\deqn{\beta = -\frac{(4\tau_3 + 2)}{(\tau_3 - 1)}\mbox{,}}
and \eqn{\alpha} then \eqn{\xi} are estimated for \emph{unknown} \eqn{\xi} as
\deqn{\alpha = \lambda_2\frac{(\beta+2)(\beta+3)}{2\beta}\mbox{, and}}
\deqn{\xi = \lambda_1 - \frac{2\alpha}{(\beta+2)}\mbox{,}}
and \eqn{\alpha} is estimated for \emph{known} \eqn{\xi} as 
\deqn{\alpha = (\lambda_1 - \xi)\frac{(\beta + 2)}{2}\mbox{.}}
The shape preservation for this distribution is an ad hoc decision. It could be that for given \eqn{\xi}, that solutions could fall back to estimating \eqn{\xi} and \eqn{\alpha} from \eqn{\lambda_1} and \eqn{\lambda_2} only. Such as solution would rely on \eqn{\tau_2 = \lambda_2/\lambda_1} with \eqn{\beta} estimated as
\deqn{\beta = \frac{3\tau_2}{(1-\tau_2)}\mbox{, and}}
\deqn{\alpha = \lambda_1\frac{(\beta+2)}{2}\mbox{,}}
but such a practice yields remarkable changes in shape for this distribution even if the provided \eqn{\xi} precisely matches that from a previous parameter estimation for which the \eqn{\xi} was treated as unknown.
}
\usage{
pargov(lmom, xi=NULL, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{xi}{An optional lower limit of the distribution. If not \code{NULL}, the \eqn{B} is still uniquely determined by \eqn{\tau_3}, the \eqn{\alpha} is adjusted so that the given lower bounds is honored. It is generally accepted to let the distribution fitting process determine its own lower bounds so \code{xi=NULL} should suffice in many circumstances.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gov}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{pargov}.}
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Nair, N.U., Sankaran, P.G., Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.

Nair, N.U., Sankaran, P.G., and Vineshkumar, B., 2012, The Govindarajulu distribution---Some Properties and applications: Communications in Statistics, Theory and Methods, v. 41, no. 24, pp. 4391--4406.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomgov}},
\code{\link{cdfgov}}, \code{\link{pdfgov}}, \code{\link{quagov}}
}
\examples{
lmr <- lmoms(rnorm(20))
pargov(lmr)

lmr <- vec2lmom(c(1391.8, 215.68, 0.01655, 0.09628))
pargov(lmr)$para             # see below
#         xi       alpha        beta 
# 868.148125 1073.740595    2.100971 
pargov(lmr, xi=868)$para     # see below
#         xi       alpha        beta 
# 868.000000 1074.044324    2.100971 
pargov(lmr, xi=100)$para     # see below
#         xi       alpha        beta 
# 100.000000 2648.817215    2.100971 
}
\keyword{distribution (parameters)}
\keyword{Distribution: Govindarajulu}
