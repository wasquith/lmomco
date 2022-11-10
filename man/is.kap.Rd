\name{is.kap}
\alias{is.kap}
\title{Is a Distribution Parameter Object Typed as Kappa}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parkap}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{kap} for the Kappa distribution.
}
\usage{
is.kap(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parkap}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{kap}.}

  \item{FALSE}{If the \code{type} is not \code{kap}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{parkap}} }
\examples{
para <- parkap(lmoms(c(123,34,4,654,37,78)))
if(is.kap(para) == TRUE) {
  Q <- quakap(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Kappa}
