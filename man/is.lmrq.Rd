\name{is.lmrq}
\alias{is.lmrq}
\title{Is a Distribution Parameter Object Typed as Linear Mean Residual Quantile Function}
\description{
The distribution parameter object returned by functions of \pkg{lmomco} such as by \code{\link{parlmrq}} are typed by an attribute \code{type}. This function checks that \code{type} is \code{lmrq} for the Linear Mean Residual Quantile Function distribution.
}
\usage{
is.lmrq(para)
}
\arguments{
  \item{para}{A parameter \code{list} returned from \code{\link{parlmrq}} or \code{\link{vec2par}}.}
}
\value{
  \item{TRUE}{If the \code{type} attribute is \code{lmrq}.}

  \item{FALSE}{If the \code{type} is not \code{lmrq}.}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{parlmrq}} }
\examples{
para <- parlmrq(lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2)))
if(is.lmrq(para) == TRUE) {
  Q <- qualmrq(0.5,para)
}
}
\keyword{utility (distribution, type check)}
\keyword{Distribution: Linear Mean Residual Quantile Function}
