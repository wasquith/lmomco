\name{pars2x}
\alias{pars2x}
\title{Estimate Quantiles from an Ensemble of Parameters}
\description{
This function acts as a frontend to estimate quantiles from an ensemble of parameters from the methods of L-moments (\code{\link{lmr2par}}), maximum likelihood (MLE, \code{\link{mle2par}}), and maximum product of spacings (MPS, \code{\link{mps2par}}) for nonexceedance probabilities. The mean, standard deviation, and number of unique quantiles for each nonexceedance probability are computed too. The unique quantiles are used because the MLE and MPS methods could fall back to L-moments or other and thus it should be considered that one of the methods might have failed.
}
\usage{
pars2x(f, paras, na.rm=FALSE, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{paras}{An ensemble of parameters from \code{\link{x2pars}}.}
  \item{na.rm}{A logical to pass to the mean and standard deviation computations.}
  \item{...}{The additional arguments, if ever used.}
}
\value{
  A \code{data.frame} having, if at least one of the parameter estimation methods is not \code{NULL}, the following columns in addition to attributes that are demonstrated in the \bold{Examples} section:
  \item{lmr}{Quantiles based on parameters from method of L-moments.}
  \item{mle}{Quantiles based on parameters from MLE.}
  \item{mps}{Quantiles based on parameters from MPS.}
  \item{f}{The nonexceedance probabilities.}
  \item{mean}{The mean of the unique quantiles (usually three) seen for each probability. Results can be affected by \code{na.rm}.}
  \item{sd}{The standard deviation of the unique quantiles (usually three) seen for each probability. Results can be affected by \code{na.rm}.}
  \item{n}{The number of unique quantiles (usually three) seen for each probability and quantiles computed as \code{NA} are not counted.}
}
\author{W.H. Asquith}
\seealso{\code{\link{x2pars}} }
\examples{
\dontrun{
# Simulate from GLO and refit it. Occasionally, the simulated data
# will result in MLE or MPS failing to converge, just a note to users.
# This example also shows the use of the attributes of the Results.
set.seed(3237)
x <- rlmomco(32, vec2par(c(2.5, 0.7, -0.39), type="glo"))
three.para.est <- x2pars(x, type="glo")
FF <- nonexceeds() # a range in nonexceedance probabilities
# In the event of MLE or MPS failure, one will see NA's in the Results.
Results <- pars2x(FF, three.para.est, na.rm=FALSE)
sum <- attr(Results, "all.summary")
plot(pp(x), sort(x), type="n", ylim=range(sum), log="y")
polygon(attr(Results, "f.poly"), attr(Results, "x.poly"), col=8, lty=0)
points(pp(x), sort(x), col=3)
lines(Results$f, Results$lmr,  col=1) # black line
lines(Results$f, Results$mle,  col=2) # red   line
lines(Results$f, Results$mps,  col=4) # blue  line
lines(Results$f, Results$mean, col=6, lty=2, lwd=2) # purple mean # }
}
\keyword{quantile function (from parameter ensemble)}

