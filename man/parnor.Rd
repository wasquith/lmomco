\name{parnor}
\alias{parnor}
\title{Estimate the Parameters of the Normal Distribution}
\description{
This function estimates the parameters of the Normal distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relation between distribution parameters and L-moments is seen under \code{\link{lmomnor}}. 

There are interesting parallels between \eqn{\lambda_2} (L-scale) and \eqn{\sigma} (standard deviation). The \eqn{\sigma} estimated from this function will not necessarily equal the output of the \code{sd} function of \R, and in fact such equality is not expected. This disconnect between the parameters of the Normal distribution and the moments (sample) of the same name can be most confusing to young trainees in statistics. The Pearson Type III is similar. See the extended example for further illustration.
}
\usage{
parnor(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}}
or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{nor}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parnor}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomnor}},
\code{\link{cdfnor}}, \code{\link{pdfnor}}, \code{\link{quanor}}
}
\examples{
lmr <- lmoms(rnorm(20))
parnor(lmr)

# A more extended example to explore the differences between an
# L-moment derived estimate of the standard deviation and R's sd()
true.std <- 15000 # select a large standard deviation
std         <- vector(mode = "numeric") # vector of sd()
std.by.lmom <- vector(mode = "numeric") # vector of L-scale values
sam <- 7   # number of samples to simulate
sim <- 100 # perform simulation sim times
for(i in seq(1,sim)) {
  Q <- rnorm(sam,sd=15000) # draw random normal variates
  std[i] <- sd(Q) # compute standard deviation
  lmr <- lmoms(Q) # compute the L-moments
  std.by.lmom[i] <- lmr$lambdas[2] # save the L-scale value
}
# convert L-scale values to equivalent standard deviations
std.by.lmom      <- sqrt(pi)*std.by.lmom

# compute the two biases and then output
# see how the standard deviation estimated through L-scale
# has a smaller bias than the usual (product moment) standard
# deviation. The unbiasness of L-moments is demonstrated.
std.bias         <- true.std - mean(std)
std.by.lmom.bias <- true.std - mean(std.by.lmom)
cat(c(std.bias,std.by.lmom.bias,"\n"))
}
\keyword{distribution (parameters)}
\keyword{Distribution: Normal}
