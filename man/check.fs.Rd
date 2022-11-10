\name{check.fs}
\alias{check.fs}
\title{Check Vector of Nonexceedance Probabilities }
\description{
This function checks that a nonexceedance probability (\eqn{F}) is in the \eqn{0 \le F \le 1} range.
It does not check that the distribution specified by parameters for \eqn{F = 0} or \eqn{F = 1} is valid. End point checking is left to additional internal checks within the functions implementing the distribution. The function is intended for internal use to build a flow of logic throughout the distribution functions. Users are not anticipated to need this function themselves. The \code{\link{check.fs}} function is separate because of the heavy use of the logic across a myriad of functions in \pkg{lmomco}.
}
\usage{
check.fs(fs)
}
\arguments{
  \item{fs}{A vector of nonexceedance probablity values.}
}
\value{
  \item{TRUE}{The nonexceedance probabilities are valid.}
  \item{FALSE}{The nonexceedance probabilities are invalid.}
}
\author{W.H. Asquith}
\seealso{\code{\link{quaaep4}},
\code{\link{quaaep4kapmix}},
\code{\link{quacau}},
\code{\link{quaemu}},
\code{\link{quaexp}},
\code{\link{quagam}},
\code{\link{quagep}},
\code{\link{quagev}},
\code{\link{quagld}},
\code{\link{quaglo}},
\code{\link{quagno}},
\code{\link{quagov}},
\code{\link{quagpa}},
\code{\link{quagum}},
\code{\link{quakap}},
\code{\link{quakmu}},
\code{\link{quakur}},
\code{\link{qualap}},
\code{\link{qualmrq}},
\code{\link{qualn3}},
\code{\link{quanor}},
\code{\link{quape3}},
\code{\link{quaray}},
\code{\link{quarevgum}},
\code{\link{quarice}},
\code{\link{quasla}},
\code{\link{quast3}},
\code{\link{quatexp}},
\code{\link{quawak}},
\code{\link{quawei}}
}
\examples{
F <- c(0.5,0.7,0.9,1.1)
if(check.fs(F) == FALSE) cat("Bad nonexceedances 0<F<1\n")
}
\keyword{utility (probability)}
