\encoding{utf8}
\name{is.pdq4}
\alias{is.pdq4}
\title{Is a Distribution Parameter Object Typed as Polynomial Density-Quantile4}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parpdq4}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{pdq4} for the Polynomial Density-Quantile4 distribution.
}
\usage{
is.pdq4(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parpdq4}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{pdq4}.}

  \item{FALSE}{If the \code{type} is not \code{pdq4}.}
}
\author{W.H. Asquith}
\seealso{ \code{\link{parpdq4}} }
\examples{
para <- parpdq4(lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52)))
if(is.pdq4(para) == TRUE) {
  Q <- quapdq4(0.5, para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Polynomial Density-Quantile4}
