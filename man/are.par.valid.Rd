\name{are.par.valid}
\alias{are.par.valid}
\title{Are the Distribution Parameters Consistent with the Distribution}
\description{
This function is a dispatcher on top of the \code{are.parCCC.valid} functions, where \code{CCC} represents the distribution type:  \code{aep4}, \code{cau}, \code{emu}, \code{exp}, \code{gam}, \code{gep}, \code{gev}, \code{glo}, \code{gno}, \code{gov}, \code{gpa}, \code{gum}, \code{kap}, \code{kmu}, \code{kur}, \code{lap}, \code{ln3}, \code{nor}, \code{pe3}, \code{ray}, \code{revgum}, \code{rice}, \code{sla}, \code{st3}, \code{texp}, \code{tri}, \code{wak}, or \code{wei}. For \pkg{lmomco} functionality, \code{\link{are.par.valid}} is called only by \code{\link{vec2par}} in the process of converting a vector into a proper distribution parameter object.
}
\usage{
are.par.valid(para, paracheck=TRUE, ...)
}
\arguments{
  \item{para}{A distribution parameter object having at least attributes \code{type} and \code{para}.}
      \item{paracheck}{A logical controlling whether the parameters are checked for validity and if \code{paracheck=TRUE} then effectively this whole function becomes turned off.}
  \item{...}{Additional arguments for the \code{are.parCCC.valid} call that is made internally.}
}
\value{
  \item{TRUE}{If the parameters are consistent with the distribution specified by the \code{type} attribute of the parameter object.}
  \item{FALSE}{If the parameters are not consistent with the distribution specified by the \code{type} attribute of the parameter object.}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{vec2par}}, \code{\link{dist.list}} }
\examples{
vec  <- c(12,120)           # parameters of exponential distribution
para <- vec2par(vec,'exp')  # build exponential distribution parameter
                            # object
# The following two conditionals are equivalent as are.parexp.valid()
# is called within are.par.valid().
if(are.par.valid(para))    Q <- quaexp(0.5,para)
if(are.parexp.valid(para)) Q <- quaexp(0.5,para)
}
\keyword{utility (distribution)}
\keyword{utility (distribution/parameter validation)}
\keyword{distribution (utility)}

