\name{is.smd}
\alias{is.smd}
\title{Is a Distribution Parameter Object Typed as Singh--Maddala}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parsmd}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{smd} for the Singh--Maddala distribution.
}
\usage{
is.smd(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parsmd}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{smd}.}

  \item{FALSE}{If the \code{type} is not \code{smd}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{parsmd}} }
\examples{
para <- parsmd(lmoms(c(123, 34, 4, 654, 37, 78)))
if(is.smd(para) == TRUE) {
  Q <- quasmd(0.5, para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
