\name{is.texp}
\alias{is.texp}
\title{Is a Distribution Parameter Object Typed as Truncated Exponential}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{partexp}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{texp} for the Truncated Exponential distribution.
}
\usage{
is.texp(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{partexp}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{texp}.}

  \item{FALSE}{If the \code{type} is not \code{texp}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{partexp}} }
\examples{
yy <- vec2par(c(123, 2.3, TRUE),  type="texp")
zz <- vec2par(c(123, 2.3, FALSE), type="texp")
if(is.texp(yy) & is.texp(zz)) {
   print(lmomtexp(yy)$lambdas)
   print(lmomtexp(zz)$lambdas)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Exponential (trimmed)}
