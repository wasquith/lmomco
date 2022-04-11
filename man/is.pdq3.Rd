\encoding{utf8}
\name{is.pdq3}
\alias{is.pdq3}
\title{Is a Distribution Parameter Object Typed as Polynomial Density-Quantile3}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parpdq3}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{pdq3} for the Polynomial Density-Quantile3 distribution.
}
\usage{
is.pdq3(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parpdq3}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{pdq3}.}

  \item{FALSE}{If the \code{type} is not \code{pdq3}.}
}
\author{W.H. Asquith}
\seealso{ \code{\link{parpdq3}} }
\examples{
para <- parpdq3(lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52)))
if(is.pdq3(para) == TRUE) {
  Q <- quapdq3(0.5, para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Polynomial Density-Quantile3}
