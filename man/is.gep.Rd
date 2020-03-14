\name{is.gep}
\alias{is.gep}
\title{Is a Distribution Parameter Object Typed as Generalized Extreme Value}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{pargep}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{gep} for the Generalized Extreme Value distribution.
}
\usage{
is.gep(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{pargep}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{gep}.}

  \item{FALSE}{If the \code{type} is not \code{gep}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pargep}} }
\examples{
#para <- pargep(lmoms(c(123,34,4,654,37,78)))
#if(is.gep(para) == TRUE) {
#  Q <- quagep(0.5,para)
#}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Generalized Exponential Poisson}
