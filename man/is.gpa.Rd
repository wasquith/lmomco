\name{is.gpa}
\alias{is.gpa}
\title{Is a Distribution Parameter Object Typed as Generalized Pareto}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{pargpa}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{gpa} for the Generalized Pareto distribution.
}
\usage{
is.gpa(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{pargpa}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{gpa}.}

  \item{FALSE}{If the \code{type} is not \code{gpa}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pargpa}} }
\examples{
para <- pargpa(lmoms(c(123,34,4,654,37,78)))
if(is.gpa(para) == TRUE) {
  Q <- quagpa(0.5,para)
}
}
\keyword{utility (distribution, type check)}
\keyword{Distribution: Generalized Pareto}
