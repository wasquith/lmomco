\name{lmom2par}
\alias{lmom2par}
\alias{lmr2par}
\title{Convert L-moments to the Parameters of a Distribution}
\description{
This function converts L-moments to the parameters of a distribution. The type of distribution is specified in the argument list:
\code{aep4}, \code{cau}, \code{emu}, \code{exp}, \code{gam}, \code{gep},
\code{gev}, \code{gld}, \code{glo}, \code{gno}, \code{gov},
\code{gpa}, \code{gum}, \code{kap}, \code{kmu}, \code{kur},
\code{lap}, \code{lmrq}, \code{ln3}, \code{nor}, \code{pe3},
\code{ray}, \code{revgum}, \code{rice}, \code{sla}, \code{st3},
\code{texp}, \code{wak}, or \code{wei}.
}
\usage{
lmom2par(lmom, type, ...)
lmr2par(x, type, ...)
}
\arguments{
  \item{lmom}{An L-moment object such as that returned by \code{\link{lmoms}} or \code{\link{pwm2lmom}}.}
  \item{type}{Three character (minimum) distribution type (for example, \code{type="gev"}).}
  \item{\ldots}{Additional arguments for the \code{parCCC} functions.}
  \item{x}{In the \code{lmr2par} call the L-moments are computed from the \eqn{x} values. This function is a parallel to \code{\link{mle2par}} and \code{\link{mps2par}}.}
}
\value{
  An \R \code{list} is returned.  This list should contain at least the following items, but some distributions such as the \code{revgum} have extra.

  \item{type}{The type of distribution in three character (minimum) format.}
  \item{para}{The parameters of the distribution.}
  \item{source}{Attribute specifying source of the parameters.}
}
\author{ W.H. Asquith}
\seealso{\code{\link{par2lmom}}}
\examples{
lmr  <- lmoms(rnorm(20))
para <- lmom2par(lmr,type="nor")

# The lmom2par() calls will error if trim != 1.
X <- rcauchy(20)
cauchy <- lmom2par(TLmoms(X, trim=1), type="cau")
slash  <- lmom2par(TLmoms(X, trim=1), type="sla")
\dontrun{
plot(pp(X), sort(X), xlab="PROBABILITY", ylab="CAUCHY")
lines(nonexceeds(), par2qua(nonexceeds(), cauchy))
lines(nonexceeds(), par2qua(nonexceeds(), slash), col=2)
}
}
\keyword{distribution (parameters)}
\keyword{parameter (L-moments to)}
\keyword{distribution (utility)}
\keyword{parameter (estimation)}

