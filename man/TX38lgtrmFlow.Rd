\name{TX38lgtrmFlow}
\docType{data}
\alias{TX38lgtrmFlow}
\title{First six L-moments of logarithms of annual mean streamflow and variances for 35 selected long-term U.S. Geological Survey streamflow-gaging stations in Texas}
\description{
L-moments of annual mean streamflow for 35 long-term U.S. Geological Survey (USGS) streamflow-gaging stations (streamgages) with at least 49 years of natural and unregulated record through water year 2012 (Asquith and Barbie, 2014). Logarithmic transformations of annual mean streamflow at each of the 35 streamgages were done. For example, logarithmic transformation of strictly positive hydrologic data is done to avoid conditional probability adjustment for the zero values; values equal to zero must be offset to avoid using a logarithm of zero. A mathematical benefit of using logarithmic transformation is that probability distributions with infinite lower and upper limits become applicable. An arbitrary value of 10 cubic feet per second was added to the streamflows for each of the 35 streamgages prior to logarithmic transformation to accommodate mean annual streamflows equal to zero (no flow). These data should be referred to as the offset-annual mean streamflow. The offsetting along the real-number line permits direct use of logarithmic transformations without the added complexity of conditional probability adjustment for zero values in magnitude and frequency analyses.

The first six sample L-moments of the base-10 logarithms of the offset-annual mean streamflow were computed using the \code{lmoms(..., nmom=6)}. The sampling variances of each corresponding L-moment are used to compute regional or study-area values for the L-moments through weighted-mean computation.  The available years of record for each of 35 stations is so large as to produce severe numerical problems in matrices needed for sampling variances using the recently developed the exact-analytical bootstrap for L-moments method (Wang and Hutson, 2013) (\code{\link{lmoms.bootbarvar}}). In order to compute sampling variances for each of the sample L-moments for each streamgage, replacement-bootstrap simulation using the \code{sample(..., replace=TRUE)} function with 10,000 replications with replacement.
}
\usage{
data(TX38lgtrmFlow)
}
\format{
An \R \code{data.frame} with
\describe{
  \item{STATION}{The USGS streamgage number.}
  \item{YEARS}{The number of years of data record.}
  \item{Mean}{The arthimetic mean (\eqn{\lambda_1}) of \eqn{\log_{10}(x + 10)}, where \eqn{x} is the vector of data.}
  \item{Lscale}{The L-scale (\eqn{\lambda_2}) of the log10-offset data.}
  \item{LCV}{The coefficient of L-variation (\eqn{\tau_2}) of the log10-offset data.}
  \item{Lskew}{The L-skew (\eqn{\tau_3}) of the log10-offset data.}
  \item{Lkurtosis}{The L-kurtosis (\eqn{\tau_4}) of the log10-offset data.}
  \item{Tau5}{The \eqn{\tau_5} of the log10-offset data.}
  \item{Tau6}{The \eqn{\tau_6} of the log10-offset data.}
  \item{VarMean}{The estimated sampling variance for \eqn{\lambda_1} multiplied by 1000.}
  \item{VarLscale}{The estimated sampling variance for \eqn{\lambda_2} multiplied by 1000.}
  \item{VarLCV}{The estimated sampling variance for \eqn{\tau_2} multiplied by 1000.}
  \item{VarLskew}{The estimated sampling variance for \eqn{\tau_3} multiplied by 1000.}
  \item{VarLkurtosis}{The estimated sampling variance for \eqn{\tau_4} multiplied by 1000.}
  \item{VarTau5}{The estimated sampling variance for \eqn{\tau_5} multiplied by 1000.}
  \item{VarTau6}{The estimated sampling variance for \eqn{\tau_6} multiplied by 1000.}
}
}
\note{
The title of this dataset indicates 35 stations, and 35 stations is the length of the data. The name of the dataset \code{TX38lgtrmFlow} and the source of the data (Asquith and Barbie, 2014) reflects 38 stations. It was decided to not show the data for 3 of the stations because a trend was detected but the dataset had already been named. The inconsistency will have to stand.
}
\references{
Asquith, W.H., and Barbie, D.L., 2014, Trend analysis and selected summary statistics of annual mean streamflow for 38 selected long-term U.S. Geological Survey streamflow-gaging stations in Texas, water years 1916--2012: U.S. Geological Survey Scientific Investigations Report 2013--5230, 16 p.

Wang, D., and Hutson, A.D., 2013, Joint confidence region estimation of L-moments with an extension to right censored data: Journal of Applied Statistics, v. 40, no. 2, pp. 368--379.
}
\examples{
data(TX38lgtrmFlow)
summary(TX38lgtrmFlow)
\dontrun{
# Need to load libraries in this order
library(lmomco); library(lmomRFA)
data(TX38lgtrmFlow)
TxDat <- TX38lgtrmFlow
TxDat <- TxDat[,-c(4)]; TxDat <- TxDat[,-c(8:15)]
summary(regtst(TxDat))
TxDat2 <- TxDat[-c(11, 28),] # Remove 08082700 Millers Creek near Munday
                             # Remove 08190500 West Nueces River at Brackettville
# No explanation for why Millers Creek is so radically discordant with the other
# streamgages with the possible exception that its data record does not span the
# drought of the 1950s like many of the other streamgages.
# The West Nueces is a highly different river from even nearby streamgages. It
# is a problem in flood frequency analysis too. So not surprizing to see this
# streamgage come up as discordant.
summary(regtst(TxDat2))
S <- summary(regtst(TxDat2))
# The results suggest that none of the 3-parameter distributions are suitable.
# The bail out solution using the Wakeby distribution is accepted. Our example
# will continue on by consideration of the two 4-parameter distributions
# available. A graphical comparison between three frequency curves will be made.
kap <- S$rpara
rmom <- S$rmom
lmr <- vec2lmom(rmom, lscale=FALSE)
aep <- paraep4(lmr)
F <- as.numeric(unlist(attributes(S$quant)$dimnames[2]))
plot(qnorm(F), S$quant[6,], type="l", lwd=3, lty=2,
     xlab="Nonexceedance probability (as standard normal variate)",
     ylab="Frequency factor (dimensionless)")
lines(qnorm(F), quakap(F, kap), col=4, lwd=2)
lines(qnorm(F), quaaep4(F, aep), col=2)
legend(-1, 0.8, c("Wakeby distribution (5 parameters)",
                  "Kappa distribution (4 parameters)",
                  "Asymmetrical Exponential Power distribution (4 parameters)"),
       bty = "n", cex=0.75, lwd=c(3,2,1), lty=c(2,1,1), col=c(1,4,2)
      )
# Based on general left tail behavior the Wakeby distribution is not acceptable.
# Based on general right tail behavior the AEP is preferred.
#
# It is recognized that the regional analysis provided by regtst() indicates
# substantial heterogeneity by all three definitions of that statistic. Further
# analysis to somehow compensate for climatological and general physiographic
# differences between the watersheds might be able to compensate for the
# heterogeneity. Such an effort is outside scope of this example.
#
# Suppose that the following data set is available for particular stream site from
# a short record streamgage, let us apply the dimensionless frequency curve as
# defined by the asymmetric exponential power distribution. Lettuce also use the
# 50-year drought as an example. This recurrence interval has a nonexceedance
# probability of 0.02. Lastly, there is the potential with this particular process
# to compute a negative annual mean streamflow, when this happens truncate to zero.
data <- c(11.9, 42.8, 36, 20.4, 43.8, 30.7, 91.1, 54.7, 43.7, 17, 28.7, 20.5, 81.2)
xbar <- mean(log10(data + 10)) # shift, log, and mean
# Note the application of the "the index method" within the exponentiation.
tmp.quantile <- 10^(xbar*quaaep4(0.02, aep)) - 10 # detrans, offset
Q50yeardrought <- ifelse(tmp.quantile < 0, 0, tmp.quantile)
# The value is 2.53 cubic feet per second average streamflow.
}
}
\keyword{datasets}
