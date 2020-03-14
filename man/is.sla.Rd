\name{is.sla}
\alias{is.sla}
\title{Is a Distribution Parameter Object Typed as Slash}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parst3}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{sla} for the Slash distribution.
}
\usage{
is.sla(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parsla}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{sla}.}

  \item{FALSE}{If the \code{type} is not \code{sla}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parsla}} }
\examples{
para <- vec2par(c(12,1.2), type='sla')
if(is.sla(para) == TRUE) {
  Q <- quasla(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Slash}
