\name{is.nor}
\alias{is.nor}
\title{Is a Distribution Parameter Object Typed as Normal}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parnor}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{nor} for the Normal distribution.
}
\usage{
is.nor(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parnor}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{nor}.}

  \item{FALSE}{If the \code{type} is not \code{nor}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{parnor}} }
\examples{
para <- parnor(lmoms(c(123,34,4,654,37,78)))
if(is.nor(para) == TRUE) {
  Q <- quanor(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Normal}
