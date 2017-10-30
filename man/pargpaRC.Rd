\name{pargpaRC}
\alias{pargpaRC}
\title{Estimate the Parameters of the Generalized Pareto Distribution with Right-Tail Censoring}
\description{
This function estimates the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) of the Generalized Pareto distribution given the \dQuote{B}-type L-moments (through the B-type probability-weighted moments) of the data under right censoring conditions (see \code{\link{pwmRC}}). The relations between distribution parameters and L-moments are seen under \code{\link{lmomgpaRC}}.
}
\usage{
pargpaRC(lmom, zeta=1, xi=NULL, lower=-1, upper=20, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{A B-type L-moment object created by a function such as \code{\link{pwm2lmom}} from B-type probability-weighted moments from \code{\link{pwmRC}}.}
  \item{zeta}{The compliment of the right-tail censoring fraction. The number of samples observed (noncensored) divided by the total number of samples.}
  \item{xi}{The lower limit of the distribution. If \eqn{\xi} is known, then alternative algorithms are used.}
  \item{lower}{The lower value for \eqn{\kappa} for a call to the \code{optimize} function. For the L-moments of the distribution to be valid \eqn{\kappa > -1}.}
  \item{upper}{The upper value for \eqn{\kappa} for a call to the \code{optimize} function. Hopefully, a large enough default is chosen for real-world data sets.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gpa}.}
  \item{para}{The parameters of the distribution.}
  \item{zeta}{The compliment of the right-tail censoring fraction.}
  \item{source}{The source of the parameters: \dQuote{pargpaRC}.}
  \item{optim}{The \code{list} returned by the \R function \code{optimize}.}
}
\details{
The \code{optimize} \R function is used to numerically solve for the shape parameter \eqn{\kappa}. No test or evaluation is made on the quality of the minimization. Users should consult the contents of the \code{optim} portion of the returned list. Finally, this function should return the same parameters if \eqn{\zeta=1} as the \code{\link{pargpa}} function.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data, in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan, chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomgpa}}, \code{\link{lmomgpaRC}}, \code{\link{pargpa}}, \code{\link{cdfgpa}}, \code{\link{pdfgpa}}, \code{\link{quagpa}}
}
\examples{
n         <- 60 # samplesize
para      <- vec2par(c(1500,160,.3),type="gpa") # build a GPA parameter set
fakedata  <- quagpa(runif(n),para) # generate n simulated values
threshold <- 1700 # a threshold to apply the simulated censoring
fakedata  <- sapply(fakedata,function(x) { if(x > threshold)
                                           return(threshold) else return(x) })
lmr       <- lmoms(fakedata) # Ordinary L-moments without considering
                             # that the data is censored
estpara   <- pargpa(lmr) # Estimated parameters of parent

pwm2     <- pwmRC(fakedata,threshold=threshold) # compute censored PWMs
typeBpwm <- pwm2$Bbetas # the B-type PWMs
zeta     <- pwm2$zeta # the censoring fraction

cenpara <- pargpaRC(pwm2lmom(typeBpwm),zeta=zeta) # Estimated parameters
F       <- nonexceeds() # nonexceedance probabilities for plotting purposes

# Visualize some data
plot(F,quagpa(F,para), type='l', lwd=3) # The true distribution
lines(F,quagpa(F,estpara), col=3) # Green estimated in the ordinary fashion
lines(F,quagpa(F,cenpara), col=2) # Red, consider that the data is censored
# now add in what the drawn sample looks like.
PP <- pp(fakedata) # plotting positions of the data
points(PP,sort(fakedata)) # sorting is needed!
# Interpretation. You should see that the red line more closely matches
# the heavy black line. The green line should be deflected to the right
# and pass through the values equal to the threshold, which reflects the
# much smaller L-skew of the ordinary L-moments compared to the type-B
# L-moments.

# Assertion, given some PWMs or L-moments, if zeta=1 then the parameter
# estimates must be identical. The following provides a demonstration.
para1 <- pargpaRC(pwm2lmom(typeBpwm),zeta=1)
para2 <- pargpa(pwm2lmom(typeBpwm))
str(para1); str(para2)

# Assertion as previous assertion, let us trigger different optimizer
# algorithms with a non-NULL xi parameter and see if the two parameter
# lists are the same.
para1 <- pargpaRC(pwm2lmom(typeBpwm), zeta=zeta)
para2 <- pargpaRC(pwm2lmom(typeBpwm), xi=para1$para[1], zeta=zeta)
str(para1); str(para2)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Generalized Pareto}
\keyword{data censoring}
