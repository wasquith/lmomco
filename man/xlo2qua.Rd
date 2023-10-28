\name{xlo2qua}
\alias{xlo2qua}
\title{Conversion of a Vector through a Left-Hand Threshold to Setup Conditional Probability Computations}
\description{
This function takes a vector of nonexceedance probabilities, a parameter object, and the object of the conditional probabability structure and computes the quantiles. This function only performs very simple vector operations. The real features for conditional probability application are found in the \code{\link{x2xlo}} and \code{\link{f2flo}} functions.
}
\usage{
xlo2qua(f, para=NULL, xlo=NULL, augasNA=FALSE,
           retrans=function(x) x, paracheck=TRUE, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}). Be aware, these are sorted internally.}
  \item{para}{Parameters from \code{\link{parpe3}} or \code{\link{vec2par}}.}
  \item{xlo}{Mandatory result from \code{\link{x2xlo}} containing the content needed for internal call to \code{\link{f2flo}} and then vector augmentation with the threshold within the \code{xlo}.}
  \item{augasNA}{A logical to switch out the threshold of \code{xlo} for \code{NA}.}
  \item{retrans}{A retransformation function for the quantiles after they are computed according to the \code{para}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity.}
  \item{...}{Additional arguments, if needed, dispatched to \code{\link{par2qua}}.}
}
\value{
  A vector of quantiles (sorted) for the nonexceedance probabilities and padding as needed to the threshold within the \code{xlo} object.
}
\author{W.H. Asquith}
\seealso{\code{\link{f2flo}}, \code{\link{flo2f}}, \code{\link{f2f}}, \code{\link{x2xlo}} }
\examples{
# This seed produces a quantile below the threshold for the FF nonexceedances and
# triggers the qua[qua <= xlo$thres] <- xlo$thres inside xlo2qua().

set.seed(2)
FF  <- nonexceeds();  LOT <- 0 # low-outlier threshold

XX  <- 10^rlmomco(20, vec2par(c(3, 0.7, 0.3), type="pe3"))
XX  <- c(rep(LOT, 5), XX)
# Pack the LOT values to the simulation, note that in most practical applications
# involving logarithms, that zeros rather than LOTs would be more apt, but this
# demonstration is useful because of the qua[qua <= xlo$thres] (see sources).
# Now, make the xlo object using the LOT as the threshold---the out of sample flag.

xlo <- x2xlo(XX, leftout=LOT)
pe3 <- parpe3( lmoms( log10(xlo$xin) ) )
# Fit the PE3 to the log10 of those values remaining in the sample.

QQ  <- xlo2qua(FF, para=pe3, xlo=xlo, retrans=function(x) 10^x)
# This line does all the work. Saves about four lines of code and streamlines
# logic when making frequency curves from the parameters and the xlo.

# Demonstrate this frequency curve to the observational sample.
plot(FF, QQ, log="y", type="l", col=grey(0.8))
points(pp(XX), sort(XX), col="red")

# Notice that with logic here and different seeds that XX could originally have
# values less than the threshold, so one would not have the lower tail all
# plotting along the threshold and a user might want to make other decisions.
QZ  <- xlo2qua(FF, para=pe3, xlo=xlo, augasNA=TRUE, retrans=function(x) 10^x)
lines(FF, QZ, col="blue")
# See how the QZ does not plot until about FF=0.2 because of the augmentation
# as NA (augasNA) being set true.
}
\keyword{utility (probability)}
\keyword{conditional probability}
\keyword{probability (conditional)}
