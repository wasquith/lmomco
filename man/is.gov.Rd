\name{is.gov}
\alias{is.gov}
\title{Is a Distribution Parameter Object Typed as Govindarajulu}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{pargov}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{gov} for the Govindarajulu distribution.
}
\usage{
is.gov(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{pargov}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{gov}.}

  \item{FALSE}{If the \code{type} is not \code{gov}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pargov}} }
\examples{
para <- pargov(lmoms(c(123,34,4,654,37,78)))
if(is.gov(para) == TRUE) {
  Q <- quagov(0.5,para)
}
}
\keyword{utility (distribution, type check)}
\keyword{Distribution: Govindarajulu}
