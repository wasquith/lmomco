\name{is.gno}
\alias{is.gno}
\title{Is a Distribution Parameter Object Typed as Generalized Normal}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{pargno}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{gno} for the Generalized Normal distribution.
}
\usage{
is.gno(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{pargno}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{gno}.}

  \item{FALSE}{If the \code{type} is not \code{gno}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pargno}} }
\examples{
para <- pargno(lmoms(c(123,34,4,654,37,78)))
if(is.gno(para) == TRUE) {
  Q <- quagno(0.5,para)
}
}
\keyword{utility (distribution, type check)}
\keyword{Distribution: Generalized Normal}
