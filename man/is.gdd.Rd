\name{is.gdd}
\alias{is.gdd}
\title{Is a Distribution Parameter Object Typed as Gamma Difference}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{pargdd}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{gdd} for the Gamma Difference distribution.
}
\usage{
is.gdd(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{pargdd}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{gdd}.}

  \item{FALSE}{If the \code{type} is not \code{gdd}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{pargdd}} }
\examples{
#
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Gamma Difference}
