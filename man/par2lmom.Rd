\name{par2lmom}
\alias{par2lmom}
\title{Convert the Parameters of a Distribution to the L-moments}
\description{
This function acts as a frontend or dispatcher to the distribution-specific L-moments of the parameter values. This function dispatches to \code{lmomCCC} where CCC represents the three character (minimum) distribution identifier:
\code{aep4}, \code{cau}, \code{emu}, \code{exp}, \code{gam}, \code{gev}, \code{gld}, \code{glo}, \code{gno}, \code{gov}, \code{gpa}, \code{gum}, \code{kap}, \code{kmu}, \code{kur}, \code{lap}, \code{lmrq}, \code{ln3}, \code{nor}, \code{pe3}, \code{ray}, \code{revgum}, \code{rice}, \code{sla}, \code{st3}, \code{texp}, \code{wak}, and \code{wei}.

The conversion of parameters to TL-moments (\code{\link{TLmoms}}) is not supported. Specific use of functions such as \code{\link{lmomTLgld}} and \code{\link{lmomTLgpa}} for the TL-moments of the Generalized Lambda and Generalized Pareto distributions is required.
}
\usage{
par2lmom(para, ...)
}
\arguments{
  \item{para}{A parameter object of a distribution.}
  \item{...}{Other arguments to pass.}
}
\value{
  An L-moment object (an \R \code{list}) is returned.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmom2par}}}
\examples{
lmr      <- lmoms(rnorm(20))
para     <- parnor(lmr)
frompara <- par2lmom(para)
}
\keyword{L-moment (parameters to)}
\keyword{L-moment (distribution)}

