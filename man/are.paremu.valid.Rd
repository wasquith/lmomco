\name{are.paremu.valid}
\alias{are.paremu.valid}
\title{Are the Distribution Parameters Consistent with the Eta-Mu Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfemu}}, \code{\link{pdfemu}}, \code{\link{quaemu}}, \code{\link{lmomemu}}), and \code{\link{lmomemu}} require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively.  These functions internally use the \code{\link{are.paremu.valid}} function. The documentation for \code{\link{pdfemu}} provides the conditions for valid parameters.
}
\usage{
are.paremu.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{paremu}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{emu} consistent.}
  \item{FALSE}{If the parameters are not \code{emu} consistent.}
}
\note{
This function calls \code{\link{is.emu}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.emu}}, \code{\link{paremu}} }
\examples{
\dontrun{
para <- vec2par(c(0.4, .04), type="emu")
if(are.paremu.valid(para)) Q <- quaemu(0.5,para) # }
}
\keyword{utility (distribution)}
\keyword{Distribution: Eta-Mu}
\keyword{utility (distribution, parameter validation)}

