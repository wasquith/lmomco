\name{par2vec}
\alias{par2vec}
\title{Convert a Parameter Object to a Vector of Parameters}
\description{
This function converts a parameter object to a vector of parameters using the \code{$para} component of the parameter list such as returned by \code{\link{vec2par}}.
}
\usage{
par2vec(para, ...)
}
\arguments{
  \item{para}{A parameter object of a distribution.}
  \item{...}{Additional arguments should they even be needed.}
}
\value{
  An \R \code{vector} is returned in moment order.
}
\author{ W.H. Asquith}
\seealso{\code{\link{vec2par}}}
\examples{
para <- vec2par(c(12,123,0.5), type="gev")
par2vec(para)
#   xi alpha kappa
# 12.0 123.0   0.5
}
\keyword{utility (vector)}
