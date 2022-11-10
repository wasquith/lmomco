\name{parkmu}
\alias{parkmu}
\title{Estimate the Parameters of the Kappa-Mu Distribution}
\description{
This function estimates the parameters (\eqn{\nu} and \eqn{\alpha}) of the Kappa-Mu (\eqn{\kappa:\mu})  distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomkmu}}.

The basic approach for parameter optimization is to extract initial guesses for the parameters from the table \code{KMU_lmompara_bykappa} in the \code{.lmomcohash} environment. The parameters having a minimum Euclidean error as controlled by three arguments are used for initial guesses in a Nelder-Mead simplex multidimensional optimization using the \R function \code{optim} and default arguments.

Limited testing indicates that of the \dQuote{error term controlling options} that the default values as shown in the Usage section seem to provide superior performance in terms of recovering the \emph{a priori known} parameters in experiments. It seems that only Euclidean optimization using L-skew and L-kurtosis is preferable, but experiments show the general algorithm to be slow.
}
\usage{
parkmu(lmom, checklmom=TRUE, checkbounds=TRUE,
         alsofitT3=FALSE, alsofitT3T4=FALSE, alsofitT3T4T5=FALSE,
         justfitT3T4=TRUE, boundary.tolerance=0.001,
         verbose=FALSE, trackoptim=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{pwm2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality).}
   \item{checkbounds}{Should the L-skew and L-kurtosis boundaries of the distribution be checked.}
   \item{alsofitT3}{Logical when true will add the error term \eqn{(\hat\tau_3 - \tau_3)^2} to the sum of square errors for the mean and L-CV.}
   \item{alsofitT3T4}{Logical when true will add the error term \eqn{(\hat\tau_3 - \tau_3)^2 + (\hat\tau_4 - \tau_4)^2} to the sum of square errors for the mean and L-CV.}
   \item{alsofitT3T4T5}{Logical when true will add the error term \eqn{(\hat\tau_3 - \tau_3)^2 + (\hat\tau_4 - \tau_4)^2 + (\hat\tau_5 - \tau_5)^2} to the sum of square errors for the mean and L-CV.}
   \item{justfitT3T4}{Logical when true will only consider the sum of squares errors for L-skew and L-kurtosis as mathematically shown for \code{alsofitT3T4}.}
   \item{boundary.tolerance}{A fudge number to help guide how close to the boundaries an arbitrary list of \eqn{\tau_3} and \eqn{\tau_4} can be to consider them formally in or out of the attainable \eqn{\{\tau_3, \tau_4\}} domain.}
   \item{verbose}{A logical to control a level of diagnostic output.}
   \item{trackoptim}{A logical to control specific messaging through each iteration of the objective function.}
   \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{kmu}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parkmu}.}
}
\references{
Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomkmu}}, \code{\link{cdfkmu}}, \code{\link{pdfkmu}}, \code{\link{quakmu}}
}
\examples{
\dontrun{
   par1 <- vec2par(c(0.7, 0.2), type="kmu")
   lmr1 <- lmomkmu(par1, nmom=4)
   par2.1 <- parkmu(lmr1, alsofitT3=TRUE,   verbose=TRUE, trackoptim=TRUE)
   par2.1$para
   par2.2 <- parkmu(lmr1, alsofitT3T4=TRUE, verbose=TRUE, trackoptim=TRUE)
   par2.2$para
   par2.3 <- parkmu(lmr1, alsofitT3=FALSE,  verbose=TRUE, trackoptim=TRUE)
   par2.3$para
   par2.4 <- parkmu(lmr1, justfitT3T4=TRUE, verbose=TRUE, trackoptim=TRUE)
   par2.4$para
   x <- seq(0,3,by=.01)
   plot(x,  pdfkmu(x, par1), type="l", lwd=6, col=8, ylim=c(0,5))
   lines(x, pdfkmu(x, par2.1), col=2, lwd=2, lty=2)
   lines(x, pdfkmu(x, par2.2), col=4)
   lines(x, pdfkmu(x, par2.3), col=3, lty=3, lwd=2)
   lines(x, pdfkmu(x, par2.4), col=5, lty=2, lwd=2)
}
\dontrun{
   par1 <- vec2par(c(1, 0.65), type="kmu")
   lmr1 <- lmomkmu(par1, nmom=4)
   par2.1 <- parkmu(lmr1, alsofitT3=TRUE,   verbose=TRUE, trackoptim=TRUE)
   par2.1$para # eta=1.0  mu=0.65
   par2.2 <- parkmu(lmr1, alsofitT3T4=TRUE, verbose=TRUE, trackoptim=TRUE)
   par2.2$para # eta=1.0  mu=0.65
   par2.3 <- parkmu(lmr1, alsofitT3=FALSE,  verbose=TRUE, trackoptim=TRUE)
   par2.3$para # eta=8.5779  mu=0.2060
   par2.4 <- parkmu(lmr1, justfitT3T4=TRUE, verbose=TRUE, trackoptim=TRUE)
   par2.4$para # eta=1.0 mu=0.65
   x <- seq(0,3,by=.01)
   plot(x,  pdfkmu(x, par1), type="l", lwd=6, col=8, ylim=c(0,1))
   lines(x, pdfkmu(x, par2.1), col=2, lwd=2, lty=2)
   lines(x, pdfkmu(x, par2.2), col=4)
   lines(x, pdfkmu(x, par2.3), col=3, lty=3, lwd=2)
   lines(x, pdfkmu(x, par2.4), col=5, lty=2, lwd=2)
   lines(x, dlmomco(x, lmom2par(lmr1, type="gam")),  lwd=2, col=2)
   lines(x, dlmomco(x, lmom2par(lmr1, type="ray")),  lwd=2, col=2, lty=2)
   lines(x, dlmomco(x, lmom2par(lmr1, type="rice")), lwd=2, col=4, lty=2)
}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Kappa-Mu}
