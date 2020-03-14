\name{is.cau}
\alias{is.cau}
\title{Is a Distribution Parameter Object Typed as Cauchy}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parcau}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{cau} for the Cauchy distribution.
}
\usage{
is.cau(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parcau}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{cau}.}

  \item{FALSE}{If the \code{type} is not \code{cau}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parcau}} }
\examples{
para <- vec2par(c(12,12),type='cau')
if(is.cau(para) == TRUE) {
  Q <- quacau(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Cauchy}
