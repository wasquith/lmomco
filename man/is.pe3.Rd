\name{is.pe3}
\alias{is.pe3}
\title{Is a Distribution Parameter Object Typed as Pearson Type III}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parpe3}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{pe3} for the Pearson Type III distribution.
}
\usage{
is.pe3(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parpe3}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{pe3}.}

  \item{FALSE}{If the \code{type} is not \code{pe3}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parpe3}} }
\examples{
para <- parpe3(lmoms(c(123,34,4,654,37,78)))
if(is.pe3(para) == TRUE) {
  Q <- quape3(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Pearson Type III}
