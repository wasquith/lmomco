\name{is.ln3}
\alias{is.ln3}
\title{Is a Distribution Parameter Object Typed as 3-Parameter Log-Normal}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parln3}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{ln3} for the 3-parameter Log-Normal distribution.
}
\usage{
is.ln3(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parln3}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{ln3}.}

  \item{FALSE}{If the \code{type} is not \code{ln3}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parln3}} }
\examples{
para <- vec2par(c(.9252, .1636, .7),type='ln3')
if(is.ln3(para)) {
  Q <- qualn3(0.5,para)
}
}
\keyword{utility (distribution/type check)}
\keyword{Distribution: Log-Normal (3-parameter)}
