\name{is.kur}
\alias{is.kur}
\title{Is a Distribution Parameter Object Typed as Kumaraswamy}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parkur}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{kur} for the Kumaraswamy distribution.
}
\usage{
is.kur(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parkur}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{kur}.}

  \item{FALSE}{If the \code{type} is not \code{kur}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{parkur}} }
\examples{
para <- parkur(lmoms(c(0.25, 0.4, 0.6, 0.65, 0.67, 0.9)))
if(is.kur(para) == TRUE) {
  Q <- quakur(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Kumaraswamy}
