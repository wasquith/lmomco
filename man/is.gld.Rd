\name{is.gld}
\alias{is.gld}
\title{Is a Distribution Parameter Object Typed as Generalized Lambda}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{pargld}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{gld} for the Generalized Lambda distribution.
}
\usage{
is.gld(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{pargld}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{gld}.}

  \item{FALSE}{If the \code{type} is not \code{gld}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{pargld}} }
\examples{
\dontrun{
para <- vec2par(c(123,120,3,2),type="gld")
if(is.gld(para) == TRUE) {
  Q <- quagld(0.5,para)
}
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Generalized Lambda}
