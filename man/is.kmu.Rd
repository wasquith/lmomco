\name{is.kmu}
\alias{is.kmu}
\title{Is a Distribution Parameter Object Typed as Kappa-Mu}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parkmu}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{kmu} for the Kappa-Mu (\eqn{\kappa:\mu}) distribution.
}
\usage{
is.kmu(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parkmu}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{kmu}.}

  \item{FALSE}{If the \code{type} is not \code{kmu}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parkmu}} }
\examples{
para <- vec2par(c(3.1, 1.4), type='kmu')
if(is.kmu(para)) {
  Q <- quakmu(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Kappa-Mu}
