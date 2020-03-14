\name{is.ray}
\alias{is.ray}
\title{Is a Distribution Parameter Object Typed as Rayleigh}
\description{
The distribution parameter object returned by functions of this module such as by \code{\link{parray}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{ray} for the Rayleigh distribution.
}
\usage{
is.ray(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parray}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{ray}.}

  \item{FALSE}{If the \code{type} is not \code{ray}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parray}} }
\examples{
para <- vec2par(c(.9252, .1636, .7),type='ray')
if(is.ray(para)) {
  Q <- quaray(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Rayleigh}
