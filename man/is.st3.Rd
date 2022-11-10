\name{is.st3}
\alias{is.st3}
\title{Is a Distribution Parameter Object Typed as 3-Parameter Student t Distribution}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parst3}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{st3} for the 3-parameter Student t distribution.
}
\usage{
is.st3(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parst3}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{st3}.}
  \item{FALSE}{If the \code{type} is not \code{st3}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{parst3}} }
\examples{
para <- vec2par(c(3, 4, 5), type='st3')
if(is.st3(para)) {
  Q <- quast3(0.25,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Student t (3-parameter)}
