\name{parkur}
\alias{parkur}
\title{Estimate the Parameters of the Kumaraswamy Distribution}
\description{
This function estimates the parameters of the Kumaraswamy distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomkur}}.
}
\usage{
parkur(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{kur}.}
  \item{para}{The parameters of the distribution.}
  \item{err}{The convergence error.}
  \item{convergence}{Logical showing whether error convergence occurred.}
  \item{source}{The source of the parameters: \dQuote{parkur}.}
}
\references{
Jones, M.C., 2009, Kumaraswamy's distribution---A beta-type distribution with some tractability advantages: Statistical Methodology, v. 6, pp. 70--81.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomkur}},
\code{\link{cdfkur}}, \code{\link{pdfkur}}, \code{\link{quakur}}
}
\examples{
lmr <- lmoms(runif(20)^2)
parkur(lmr)

kurpar <- list(para=c(1,1), type="kur");
lmr <- lmomkur(kurpar)
parkur(lmr)

kurpar <- list(para=c(0.1,1), type="kur");
lmr <- lmomkur(kurpar)
parkur(lmr)

kurpar <- list(para=c(1,0.1), type="kur");
lmr <- lmomkur(kurpar)
parkur(lmr)

kurpar <- list(para=c(0.1,0.1), type="kur");
lmr <- lmomkur(kurpar)
parkur(lmr)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Kumaraswamy}
