\name{is.tri}
\alias{is.tri}
\title{Is a Distribution Parameter Object Typed as Asymmetric Triangular}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{partri}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{tri} for the Asymmetric Triangular distribution.
}
\usage{
is.tri(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{partri}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{tri}.}

  \item{FALSE}{If the \code{type} is not \code{tri}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{partri}} }
\examples{
para <- partri(lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52)))
if(is.tri(para) == TRUE) {
  Q <- quatri(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}
