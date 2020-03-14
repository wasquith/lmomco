\name{is.rice}
\alias{is.rice}
\title{Is a Distribution Parameter Object Typed as Rice}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parrice}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{rice} for the Rice distribution.
}
\usage{
is.rice(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parrice}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{rice}.}

  \item{FALSE}{If the \code{type} is not \code{rice}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parrice}} }
\examples{
para <- vec2par(c(3, 4),type='rice')
if(is.rice(para)) {
  Q <- quarice(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Rice}
