\name{gen.freq.curves}
\alias{gen.freq.curves}
\title{Plot Randomly Generated Frequency Curves from a Parent Distribution}
\description{
This function generates random samples of specified size from a specified parent distribution. Subsequently, the type of parent distribution is fit to the L-moments of the generated sample. The fitted distribution is then plotted. It is the user's responsibility to have an active \code{plot} already drawn; unless the \code{callplot} option is \code{TRUE}. This function is useful to demonstration of sample size on the uncertainty of a fitted distribution---a motivation for this function is as a classroom exercise.
}
\usage{
gen.freq.curves(n, para, F=NULL, nsim=10, callplot=TRUE, aslog=FALSE,
                asprob=FALSE, showsample=FALSE, showparent=FALSE,
                lowerCI=NA, upperCI=NA, FCI=NA, ...)
}
\arguments{
  \item{n}{Sample size to draw from parent as specified by \code{para}.}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{F}{The nonexceedance probabilities for horizontal axis---defaults to \code{\link{nonexceeds}} when the argument is \code{NULL}.}
  \item{nsim}{The number of simulations to perform (frequency curves to draw)---the default is 10.}
  \item{callplot}{Calls \code{plot} to acquire a graphics device---default is \code{TRUE}, but the called \code{plot} is left empty.}
  \item{aslog}{Compute \code{log10} of quantiles---note that \preformatted{NaNs produced in: log(x, base)} will be produced for less than zero values. Otherwise this is a harmless message.}
  \item{asprob}{The \code{qnorm} function is used to convert nonexceedance probabilities, which are produced by \code{\link{nonexceeds}}, to Standard Normal variates. The Normal distribution will be a straight line when this argument is \code{TRUE} and \code{aslog=FALSE}.}
  \item{showsample}{Each simulated sample is drawn through plotting positions (\code{\link{pp}}).}
  \item{showparent}{The curve for the parent distribution is plotted on exit from the function if \code{TRUE}. Further plotting options can not be controlled---unlike the situation with the drawing of the simulated frequency curves.}
  \item{lowerCI}{An optional estimate of the lower confidence limit for the \code{FCI} nonexceedance probability.}
  \item{upperCI}{An optional estimate of the upper confidence limit for the \code{FCI} nonexceedance probability.}
  \item{FCI}{The nonexeedance probability of interest for the confidence limits provided in \code{lowerCI} and \code{upperCI}.}
  \item{...}{Additional parameters are passed to the \code{lines} call within the function---except for the drawing of the parent distribution (see argument \code{showparent}).}
}
\value{
  This function is largely used for its graphical side effects, but if estimates of the lower and upper confidence limits are known (say from \code{\link{genci.simple}}) then this function can be used to evaluate the counts of simulations at nonexceedance probability \code{FCI} outside the limits provided in \code{lowerCI} and \code{upperCI}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmom2par}}, \code{\link{nonexceeds}}, \code{\link{rlmomco}}, \code{\link{lmoms}}
}
\examples{
\dontrun{
# 1-day rainfall Travis county, Texas
para <- vec2par(c(3.00, 1.20, -.0954), type="gev")
F <- .99 # the 100-year event
n <- 46 # sample size derived from 75th percentile of record length distribution
# for Edwards Plateau from Figure 3 of USGS WRIR98-4044 (Asquith, 1998)
# Argument for 75th percentile is that the contours of distribution parameters
# in that report represent a regionalization of the parameters and hence
# record lengths such as the median or smaller for the region seem too small
# for reasonable exploration of confidence limits of precipitation.
nsim <- 5000 # simulation size
seed <- runif(1, min=1, max=10000)
set.seed(seed)
CI <- genci.simple(para, n, F=F, nsim=nsim, edist="nor")
lo.nor <- CI$lower; hi.nor <- CI$upper

set.seed(seed)
CI <- genci.simple(para, n, F=F, nsim=nsim, edist="aep4")
lo.aep4 <- CI$lower; hi.aep4 <- CI$upper
message("NORMAL ERROR DIST: lowerCI = ",lo.nor, " and upperCI = ",hi.nor)
message("  AEP4 ERROR DIST: lowerCI = ",lo.aep4," and upperCI = ",hi.aep4)
qF <- qnorm(F)
# simulated are grey, parent is black
set.seed(seed)
counts.nor  <- gen.freq.curves(n, para, nsim=nsim,
                   asprob=TRUE, showparent=TRUE, col=rgb(0,0,1,0.025),
                   lowerCI=lo.nor, upperCI=hi.nor, FCI=F)
set.seed(seed)
counts.aep4 <- gen.freq.curves(n, para, nsim=nsim,
                   asprob=TRUE, showparent=TRUE, col=rgb(0,0,1,0.025),
                   lowerCI=lo.aep4, upperCI=hi.aep4, FCI=F)
lines( c(qF,qF), c(lo.nor, hi.nor),  lwd=2, col=2)
points(c(qF,qF), c(lo.nor, hi.nor),  pch=1, lwd=2, col=2)
lines( c(qF,qF), c(lo.aep4,hi.aep4), lwd=2, col=2)
points(c(qF,qF), c(lo.aep4,hi.aep4), pch=2, lwd=2, col=2)
percent.nor  <- (counts.nor$count.above.upperCI +
                 counts.nor$count.below.lowerCI) /
                 counts.nor$count.valid.simulations
percent.aep4 <- (counts.aep4$count.above.upperCI +
                 counts.aep4$count.below.lowerCI) /
                 counts.aep4$count.valid.simulations
percent.nor  <- 100 * percent.nor
percent.aep4 <- 100 * percent.aep4
message("NORMAL ERROR DIST: ",percent.nor)
message("  AEP4 ERROR DIST: ",percent.aep4)
# Continuing on, we are strictly focused on F being equal to 0.99
# Also we are no restricted to the example using the GEV distribution
# The vargev() function is from Handbook of Hydrology
"vargev" <-
function(para, n, F=c("F080", "F090", "F095", "F099", "F998", "F999")) {
   F <- as.character(F)
   if(! are.pargev.valid(para)) return()
   F <- match.arg(F)
   A <- para$para[2]
   K <- para$para[3]
   AS <- list(F080=c(-1.813,  3.017, -1.4010, 0.854),
              F090=c(-2.667,  4.491, -2.2070, 1.802),
              F095=c(-3.222,  5.732, -2.3670, 2.512),
              F098=c(-3.756,  7.185, -2.3140, 4.075),
              F099=c(-4.147,  8.216, -0.2033, 4.780),
              F998=c(-5.336, 10.711, -1.1930, 5.300),
              F999=c(-5.943, 11.815, -0.6300, 6.262))
   AS <- as.environment(AS); CO <- get(F, AS)
   varx <- A^2 * exp( CO[1] + CO[2]*exp(-K) + CO[3]*K^2 + CO[4]*K^3 ) / n
   names(varx) <- NULL
   return(varx)
}
sdx <- sqrt(vargev(para, n, F="F099"))
VAL  <- qlmomco(F, para)
lo.vargev <- VAL + qt(0.05, df=n) * sdx # minus covered by return of qt()
hi.vargev <- VAL + qt(0.95, df=n) * sdx

set.seed(seed)
counts.vargev <- gen.freq.curves(n, para, nsim=nsim,
                   xlim=c(0,3), ylim=c(3,15),
                   asprob=TRUE, showparent=TRUE, col=rgb(0,0,1,0.01),
                   lowerCI=lo.vargev, upperCI=hi.vargev, FCI=F)
percent.vargev  <- (counts.vargev$count.above.upperCI +
                    counts.vargev$count.below.lowerCI) /
                    counts.vargev$count.valid.simulations
percent.vargev  <- 100 * percent.vargev
lines(c(qF,qF),  range(c(lo.nor,   hi.nor,
                         lo.aep4,  hi.aep4,
                         lo.vargev,hi.vargev)), col=2)
points(c(qF,qF), c(lo.nor,      hi.nor), pch=1, lwd=2, col=2)
points(c(qF,qF), c(lo.aep4,    hi.aep4), pch=3, lwd=2, col=2)
points(c(qF,qF), c(lo.vargev,hi.vargev), pch=2, lwd=2, col=2)
message("NORMAL ERROR DIST: ",percent.nor)
message("  AEP4 ERROR DIST: ",percent.aep4)
message("VARGEV ERROR DIST: ",percent.vargev)
}
}
\keyword{misc}
\keyword{Graphics}
