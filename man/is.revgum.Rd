\name{is.revgum}
\alias{is.revgum}
\title{Is a Distribution Parameter Object Typed as Reverse Gumbel}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parrevgum}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{revgum} for the Reverse Gumbel distribution.
}
\usage{
is.revgum(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parrevgum}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{revgum}.}

  \item{FALSE}{If the \code{type} is not \code{revgum}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parrevgum}} }
\examples{
para <- vec2par(c(.9252, .1636, .7),type='revgum')
if(is.revgum(para)) {
  Q <- quarevgum(0.5,para)
}
}
\keyword{utility (distribution, type check)}
\keyword{Distribution: Gumbel (reversed)}
