\name{par2qua2lo}
\alias{par2qua2lo}
\title{ Equivalent Quantile Function of Two Distributions Stemming from Left-Hand Threshold to Setup Conditional Probability Computations}
\description{
\bold{EXPERIMENTAL!} This function computes the nonexceedance probability of a given quantile from a linear weighted combination of two quantile functions---a mixed distribution---when the data have been processed through the \code{\link{x2xlo}} function setting up left-hand thresholding and conditional probability compuation. The \code{par2qua2lo} function is a partial generalization of the \code{\link{par2qua2}} function (see there for the basic mathematics). The \bold{Examples} section has an exhaustive demonstration. The resulting weighted- or mixed-quantile function is not rigorously checked for monotonic increase with \eqn{F}, which is a required property of quantile functions. However, a first-order difference on the mixed quantiles with the probabilities is computed and a warning issued if not monotonic increasing.
}
\usage{
par2qua2lo(f, para1, para2, xlo1, xlo2,
              wfunc=NULL, weight=NULL, addouts=FALSE,
              inf.as.na=TRUE, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para1}{The first distribution parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{para2}{The second distribution parameters from \code{\link{x2xlo}}.}
  \item{xlo1}{The first distribution parameters from \code{\link{x2xlo}}.}
  \item{xlo2}{The second distribution parameters from \code{\link{lmom2par}} or similar.}
  \item{wfunc}{A function taking the argument \code{f} and computing a weight for the \code{para2} curve for which the complement of the computed weight is used for the weight on \code{para1}.}
  \item{weight}{An optional weighting argument to use in lieu of \code{F}. If \code{NULL} then weights are a function of \code{length(xlo1$xin)} and \code{length(xlo2$xin)} for the first and second distribution respectively, if \code{weight} has length 1, then weight on first distribution is the complement of the weight, and the weight on second distribution is \code{weight[1]}, and if \code{weight} had length 2, then \code{weight[1]} is the weight on the first distribution, and \code{weight[2]} is the weight on the second distribution.}
  \item{addouts}{In the computation of weight factors when the \code{xlo1$xin} and \code{xlo2$xin} are used by other argument settings, the \code{addouts} arguments triggers the inclusion of the lengths of the \code{xlo1$xout} and \code{xlo2$xout} (see source code).}
  \item{inf.as.na}{A logical controlling whether quantiles for each distribution that are non-finite are to be converted to \code{NA}s. If they are converter to \code{NA}s, then when the application of the weight or weights are made then that those indices of \code{NA} quantiles become a zero and the weight for the other quantile will become unity. It is suggested to review the source code.}
  \item{...}{Additional arguments to pass if needed.}
}
\value{
  The mixed quantile values for likely a subset of the provided \code{f} from the two distributions depending on the internals of \code{xlo1} and \code{xlo2} require the quantiles to actually start. This requires this function to return an \R \code{data.frame} that was only optional for \code{\link{par2qua2}}:
 \item{f}{Nonexceedance probabilities.}
 \item{quamix}{The mixed quantiles.}
 \item{delta_curve1}{The computation \code{quamix} minus curve for \code{para1}.}
 \item{delta_curve2}{The computation \code{quamix} minus curve for \code{para2}.}

Alternatively, the returned value could be a weighting function for subsequent calls as \code{wfunc} to \code{par2qua2lo} (see \bold{Examples}). This alternative operation is triggered by setting \code{wfunc} to an arbitrary character string, and internally the contents of \code{xlo1} and \code{xlo2}, which themselves have to be called as named arguments, are recombined. This means that the \code{xin} and code{xout} are recombined, into their respective samples. Each data point is then categorized with probability zero for the \code{xlo1} values and probability unity for the \code{xlo2} values. A logistic regression is fit using logit-link function for a binomial family using a generalized linear model. The binomial (0 or 1) is regressed as a function of the plotting positions of a sample composed of \code{xlo1} and \code{xlo2}. The coefficients of the regression are extracted, and a function created to predict the probability of event \dQuote{\code{xlo2}}. The \code{attributes} of the computed value inside the function store the coefficients, the regression model, and potentially useful for graphical review, a \code{data.frame} of the data used for the regression. This sounds more complicated than it really is (see source code and \bold{Examples}).
}
\author{ W.H. Asquith}
\seealso{\code{\link{par2qua}}, \code{\link{par2cdf2}}, \code{\link{par2qua2}}, \code{\link{x2xlo}}
}
\examples{
\dontrun{
XloSNOW <- list( # data from "snow events" from prior call to x2xlo()
   xin=c(4670, 3210, 4400, 4380, 4350, 3380, 2950, 2880, 4100),
   ppin=c(0.9444444, 0.6111111, 0.8888889, 0.8333333, 0.7777778, 0.6666667,
          0.5555556, 0.5000000, 0.7222222),
   xout=c(1750, 1610, 1750, 1460, 1950, 1000, 1110, 2600),
   ppout=c(0.27777778, 0.22222222, 0.33333333, 0.16666667, 0.38888889,
           0.05555556, 0.11111111, 0.44444444),
   pp=0.4444444, thres=2600, nin=9, nout=8, n=17, source="x2xlo")
# RAIN data from prior call to x2xlo() are
XloRAIN <- list( # data from "rain events" from prior call to x2xlo()
   xin=c(5240, 6800, 5990, 4600, 5200, 6000, 4500, 4450, 4480, 4600,
         3290, 6700, 10600, 7230, 9200, 6540, 13500, 4250, 5070,
         6640, 6510, 3610, 6370, 5530, 4600, 6570, 6030, 7890, 8410),
   ppin=c(0.41935484, 0.77419355, 0.48387097, 0.25806452, 0.38709677, 0.51612903,
          0.22580645, 0.16129032, 0.19354839, 0.29032258, 0.06451613, 0.74193548,
          0.93548387, 0.80645161, 0.90322581, 0.64516129, 0.96774194, 0.12903226,
          0.35483871, 0.70967742, 0.61290323, 0.09677419, 0.58064516, 0.45161290,
          0.32258065, 0.67741935, 0.54838710, 0.83870968, 0.87096774),
   xout=c(1600), ppout=c(0.03225806),
   pp=0.03225806, thres=2599, nin=29, nout=1, n=30, source="x2xlo")

QSNOW <- c(XloSNOW$xin,  XloSNOW$xout ) # collect all of the snow
QRAIN <- c(XloRAIN$xin,  XloRAIN$xout ) # collect all of the rain
PSNOW <- c(XloSNOW$ppin, XloSNOW$ppout) # probabilities collected
PRAIN <- c(XloRAIN$ppin, XloRAIN$ppout) # probabilities collected

# Logistic regression to blend the proportion of snow versus rain events as
# ***also*** a function of nonexceedance probability
wfunc <- par2qua2lo(xlo1=XloSNOW, xlo2=XloRAIN, wfunc="wfunc") # weight function

# Plotting the data and the logistic regression. This shows how to gain access
# to the attributes, in order to get the data, so that we can visualize the
# probability mixing between the two samples. If the two samples are not a
# function of probability, then each systematically would have a regression-
# predicted weight of 50/50. For the RAIN and SNOW, the SNOW is likely to
# produce the smaller events and RAIN the larger.
 opts <- par(las=1) # Note the 0.5 in the next line is arbitrary, we simply
 bin <- attr(wfunc(0.5), "data") # have to use wfunc() to get its attributes.
 FF <- seq(0,1,by=0.01); HH <- wfunc(FF); n <- length(FF)
 plot(bin$f, bin$prob, tcl=0.5, col=2*bin$prob+2,
      xlab="NONEXCEEDANCE PROBABILITY", ylab="RAIN-CAUSED EVENT RELATIVE TO SNOW")
 lines(c(-0.04,1.04), rep(0.5,2), col=8, lwd=0.8) # origin line at 50/50 chance
 text(0, 0.5, "50/50 chance line", pos=4, cex=0.8)
 segments(FF[1:(n-1)], HH[1:(n-1)], x1=FF[2:n], y1=HH[2:n], lwd=1+4*abs(FF-0.5),
          col=rgb(1-FF,0,FF)) # line grades from one color to other
 text(1, 0.1, "Events caused by snow", col=2, cex=0.8, pos=2)
 text(0, 0.9, "Events caused by rain", col=4, cex=0.8, pos=4)
 par(opts)

# Suppose that the Pearson type III is thought applicable to the SNOW
# and the AEP4 for the RAIN, now estimate respective parameters.
parSNOW <- lmr2par(log10(XloSNOW$xin), type="nor" )
parRAIN <- lmr2par(log10(XloRAIN$xin), type="wak")
# Two distributions are chosen to show the user than we are not constrained to one.

Qall   <- c(QSNOW, QRAIN)                # combine into a "whole" sample
XloALL <- x2xlo(Qall, leftout=2600, a=0) # apply the low-outlier threshold
parALL <- lmr2par(log10(XloALL$xin), type="nor") # estimate Wakeby
# Wakey has five parameters and is very flexible.

FF <- nonexceeds() # useful nonexceedance probabilities
col <- c(rep(0,length(QSNOW)), rep(2,length(QRAIN))) # for coloring
plot(0, 0, col=2+col, ylim=c(1000,20000), xlim=qnorm(range(FF)), log="y",
           xlab="STANDARD NORMAL VARIATE", ylab="QUANTILE", type="n")
lines(par()$usr[1:2], rep(2600, 2), col=6, lty=2, lwd=0.5) # draw threshold
points(qnorm(pp(Qall, sort=FALSE)), Qall, col=2+col, lwd=0.98) # all record
points(qnorm(PSNOW), QSNOW, pch=16, col=2) # snow events
points(qnorm(PRAIN), QRAIN, pch=16, col=4) # rain events
lines(     qnorm(f2f(  FF, xlo=XloSNOW)), # show fitted curve for snow events
      10^par2qua(f2flo(FF, xlo=XloSNOW ), parSNOW), col=2)
lines(     qnorm(f2f(  FF, xlo=XloRAIN)), # show fitted curve for rain events
      10^par2qua(f2flo(FF, xlo=XloRAIN ), parRAIN), col=4)
lines(     qnorm(f2f(  FF, xlo=XloALL )), # show fitted curve for all events combined
      10^par2qua(f2flo(FF, xlo=XloALL  ), parALL ), col=1, lty=3)
PQ <- par2qua2lo(      FF, parSNOW, parRAIN, XloSNOW, XloRAIN, wfunc=wfunc)
lines(qnorm(PQ$f), 10^PQ$quamix, lwd=2)                  # draw the mixture
legend(-3,20000, c("Rain curve", "Snow curve", "All combined (all open circles)",
                    "MIXED CURVE by par2qua2lo()"),
                  bty="n", lwd=c(1,1,1,2), lty=c(1,1,3,1), col=c(4,2,1,1))
text(-3, 15000, "A low-outlier threshold of 2,600 is used throughout.", col=6, pos=4)
text(-3,  2600, "2,600", cex=0.8, col=6, pos=4)
mtext("Mixed population frequency computation of snow and rainfall streamflow")#}

\dontrun{
nsim <- 50000; FF <- runif(nsim); WF <- wfunc(FF)
rB <- rbinom(nsim, 1, WF)
RF <- FF[rB == 1]; SF <- FF[rB == 0]
RAIN <- 10^qlmomco(f2flo(runif(length(RF)), xlo=XloRAIN), parRAIN)
SNOW <- 10^qlmomco(f2flo(runif(length(SF)), xlo=XloRAIN), parSNOW)
RAIN[RAIN < XloRAIN$thres] <- XloRAIN$thres
SNOW[SNOW < XloSNOW$thres] <- XloSNOW$thres
RAIN <- c(RAIN,rep(XloRAIN$thres, length(RF)-length(RAIN)))
SNOW <- c(SNOW,rep(XloSNOW$thres, length(SF)-length(SNOW)))
ALL <- c(RAIN,SNOW)
lines(qnorm(pp(ALL)), sort(ALL), cex=0.6, lwd=0.8, col=3)

RF <- FF[rB == 1]; SF <- FF[rB == 0]
RAIN <- 10^qlmomco(RF, parRAIN)
SNOW <- 10^qlmomco(SF, parSNOW)
RAIN[RAIN < XloRAIN$thres] <- XloRAIN$thres
SNOW[SNOW < XloSNOW$thres] <- XloSNOW$thres
RAIN <- c(RAIN,rep(XloRAIN$thres, length(RF)-length(RAIN)))
SNOW <- c(SNOW,rep(XloSNOW$thres, length(SF)-length(SNOW)))
ALL <- c(RAIN,SNOW)
lines(qnorm(pp(ALL)), sort(ALL), cex=0.6, lwd=0.8, col=3)

RF <- FF[rB == 1]; SF <- FF[rB == 0]
RAIN <- 10^qlmomco(f2flo(RF, xlo=XloRAIN), parRAIN)
SNOW <- 10^qlmomco(f2flo(SF, xlo=XloRAIN), parSNOW)
RAIN[RAIN < XloRAIN$thres] <- XloRAIN$thres
SNOW[SNOW < XloSNOW$thres] <- XloSNOW$thres
RAIN <- c(RAIN,rep(XloRAIN$thres, length(RF)-length(RAIN)))
SNOW <- c(SNOW,rep(XloSNOW$thres, length(SF)-length(SNOW)))
ALL <- c(RAIN,SNOW)
lines(qnorm(pp(ALL)), sort(ALL), cex=0.6, lwd=0.8, col=3) #}
}
\keyword{quantile mixture function}
\keyword{distribution (mixture)}
\keyword{mixed distribution}
\keyword{conditional probability}
\keyword{probability (conditional)}
\keyword{logistic regression}
\keyword{generalized linear model}
\keyword{binomial family}
\keyword{logit link}
