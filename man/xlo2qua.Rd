\name{xlo2qua}
\alias{xlo2qua}
\title{Conversion of a Vector through a Left-Hand Threshold to Setup Conditional Probability Computations}
\description{
This function takes a vector of nonexceedance probabilities, a parameter object, and the object of the conditional probabability structure and computes the quantiles. This function only performs very simple vector operations. The real features for conditional probability application are found in the \code{\link{x2xlo}} and \code{\link{f2flo}} functions.
}
\usage{
xlo2qua(f, para=NULL, xlo=NULL, augasNA=FALSE, sort=FALSE, fillthres=TRUE,
           retrans=function(x) x, paracheck=TRUE, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}). Be aware, these are sorted internally.}
  \item{para}{Parameters from \code{\link{parpe3}} or \code{\link{vec2par}}.}
  \item{xlo}{Mandatory result from \code{\link{x2xlo}} containing the content needed for internal call to \code{\link{f2flo}} and then vector augmentation with the threshold within the \code{xlo}. If this is left as \code{NULL}, then the function simply calls the quantile function for the parameters in \code{para}.}
  \item{augasNA}{A logical to switch out the threshold of \code{xlo} for \code{NA}.}
  \item{sort}{A logical whose default adheres to long-term assembly of \pkg{lmomco} behavior with working with conditional trunction. Setting this to true, triggers hand assembly of the the unsorted returned quantiles with support for \code{NA} and more flexibility than \code{\link{x2xlo}} as originally designed. If sort is true, then the \code{f} is permitted to contain \code{NA} values.}
  \item{fillthres}{A logical to trigger \code{qua[qua <= xlo$thres] <- xlo$thres} or replacement of computed values less than the threshold with the threshold. The argument \code{augasNA} is consulted after \code{fillthres}.}
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

\dontrun{
# Needs library(copBasic); library(MGBT) # too
Asite <- "08148500"; Bsite <- "08150000"; dtype <- "gev"
AB    <- MGBT::jointPeaks(Asite, Bsite) # tables of the peaks and pairwise peaks
A     <- AB$Asite_no[AB$Asite_no$appearsSystematic == TRUE, ] # only record when
B     <- AB$Bsite_no[AB$Bsite_no$appearsSystematic == TRUE, ] # monitoring occurring
QA    <- A$peak_va; Alot <- 0 # cfs (just protection from zeros, more sophisticated)
QB    <- B$peak_va; Blot <- 0 # cfs (work might be needed for better thresholds)
Alo   <- x2xlo(QA, leftout=Alot) # A xlo object
Blo   <- x2xlo(QB, leftout=Blot) # B xlo object
Apara <- lmr2par(log10(Alo$xin), type=dtype) # note log10
Bpara <- lmr2par(log10(Blo$xin), type=dtype) # note log10
Aupr  <- 10^supdist(Apara)$support[2]
Bupr  <- 10^supdist(Bpara)$support[2]
UVsS  <- AB$AB[, c("U", "V")] # isolate paired empirical probabilities
rhoS  <- copBasic::rhoCOP(as.sample=TRUE,     para=UVsS) # Spearman rho
infS  <- copBasic::LzCOPpermsym(cop=EMPIRcop, para=UVsS, as.vec=TRUE)
# a vector of permutation (variable exchangability) distances

tparf <- function(par) { c(log(par[1] -1), log(par[2]),  # transform for optimization
                   qnorm(punif(par[3],  min=-1, max=1))) }
rparf <- function(par) { c(exp(par[1])+1,  exp(par[2]),  # re-transformation to copula
                   qunif(pnorm(par[3]), min=-1, max=1)) }

ofunc <- function(par) { # objective function
  mypara <- rparf(par)   # re-transform to copula space
  mypara <- list(cop=GHcop, para=mypara[1:2], breve=mypara[3]) # asymmetry by breveCOP()
  rhoT   <- copBasic::rhoCOP(cop=breveCOP, para=mypara) # Spearman rho
  infT   <- copBasic::LzCOPpermsym(cop=breveCOP, para=mypara, as.vec=TRUE)
  err    <- mean( (infT - infS)^2 ) + (rhoT - rhoS)^2 # sum of square-like errors
  return(err)
}
init.par <- tparf(c(2, 1, 0)); rt <- NULL # init parameters and root
try( rt <- optim(init.par, ofunc) )
cpara <- rparf(rt$par) # re-transformation
cpara <- list(cop=GHcop, para=cpara[1:2], breve=cpara[3]) # copula parameters for
# an double-parameter Gumbel copula with permutation asymmetry via the breve.

ns <- 1000 # years of bivariate simulation
UVsim <- copBasic::rCOP(ns, cop=breveCOP, para=cpara, resamv01=TRUE) # simulation
AS <- xlo2qua(UVsim[,1], para=Apara, xlo=Alo, sort=FALSE,  # **** see xlo2qua() use
                         retrans=function(x) 10^x, paracheck=FALSE)
BS <- xlo2qua(UVsim[,2], para=Bpara, xlo=Blo, sort=FALSE,  # **** see xlo2qua() use
                         retrans=function(x) 10^x, paracheck=FALSE)

FF  <- seq(0.001, 0.999, by=0.001); qFF <- qnorm(FF) # probabilities for marginal curve
AF <- xlo2qua(FF, para=Apara, xlo=Alo, sort=FALSE,         # **** see xlo2qua() use
                  retrans=function(x) 10^(x), paracheck=FALSE)
BF <- xlo2qua(FF, para=Bpara, xlo=Blo, sort=FALSE,         # **** see xlo2qua() use
                  retrans=function(x) 10^(x), paracheck=FALSE)
# There might be a small region in the lower-left corner that is not attainable by the
# use of the thresholding. Let us add the complexity to the example by working out
# about the minimum points on the curves w/o more sophisticated computation.
mx <- min(c(AS, AF), na.rm=TRUE); my <- min(c(BS, BF), na.rm=TRUE)
# The use of the mx and my help us with a polygon to come, but also help us to set
# some axis limits that are especially suitable to see the entire situation of the
# simulation canvasing [0,1]^2 but the quantiles through the univariate margins might
# have truncation because of handling of the lower-tail by the threshold.

# finally plot the bivariate relation
plot(AB$AB$Apeak_va, AB$AB$Bpeak_va, log="xy", type="n",
     xlim=range(c(mx, QA, AS, ifelse(is.finite(Aupr), Aupr, NA)), na.rm=TRUE),
     ylim=range(c(my, QB, BS, ifelse(is.finite(Bupr), Bupr, NA)), na.rm=TRUE),
     xlab=paste0("Paired water-year peak streamflow for streamgage ", Asite),
     ylab=paste0("Paired water-year peak streamflow for streamgage ", Bsite))
cr <- 10^par()$usr[c(1, 3)]             # finish forming the region in the lower-left
px <- c(cr[1], mx, mx, cr[1], cr[1])    # corner that is truncated away; we do this
py <- c(cr[2], cr[2], my, my, cr[2])    # this because log10() used and in practical
polygon(px, py, col="wheat", border=NA) # applications at best zeros might be data
abline(v=mx, lty=2, lwd=0.8); abline(h=my, lty=2, lwd=0.8) # further demarcation
if( is.finite(Aupr) ) abline(v=Aupr, lty=2, lwd=1.5, col="purple") # upper limit
if( is.finite(Bupr) ) abline(h=Bupr, lty=2, lwd=1.5, col="purple") # upper limit
points(AS, BS, pch=21, col="red", bg="white") # now plot the simulations
points(AB$AB$Apeak_va, AB$AB$Bpeak_va, cex=AB$AB$cex, # now plot the observed data that
       col="black", bg=grey(AB$AB$cex/2), pch=21) # defined the parameter estimation of
legend("bottomright",                             # the copula then draw a legend.
     c("Paired streamflow (fill lightens/size increases as days apart increases)",
       paste0(ns, " years simulated by copula and GEV margins")), bty="o", cex=0.8,
       pch=c(21,21), col=c("black","red"), pt.cex=c(1.3,1), pt.bg=c(grey(0.7),"white"))

ST <- round(1/(1-kfuncCOP(0.99, cop=breveCOP, para=cpara)), digits=0)
message("Super-critical return period for ",
               "primary return period of 100 years is ", ST, " years.")

#  move on to showing the univariate margins by parametric fit with left-truncation
plot(qnorm(pp(QA)), sort(QA), log="y", pch=21, bg="white", main=Asite,
     ylim=range(c(QA, AF, Aupr), na.rm=TRUE),
     xlab="Standard normal variate", ylab="Peak streamflow, in cfs")
abline(h=Aupr, lty=2, lwd=1.5, col="purple")
lines(qFF, AF, lwd=3, col="seagreen")
legend("bottomright",
     c(paste0("Marginal distribution by ", toupper(dtype)),
       "Upper bounds of fitted distribution",
       "Systematic peaks by Weibull plotting position"), bty="o", seg.len=3,
       pch=c(NA,NA,21), col=c("seagreen","purple","black"), bg="white", cex=0.8,
       lty=c(1, 2, NA), lwd=c(3, 1.5, NA), pt.bg=c(NA, NA, "white"))

plot(qnorm(pp(QB)), sort(QB), log="y", pch=21, bg="white", main=Bsite,
     ylim=range(c(QB, BF, Bupr), na.rm=TRUE),
     xlab="Standard normal variate", ylab="Peak streamflow, in cfs")
abline(h=Bupr, lty=2, lwd=1.5, col="purple")
lines(qFF, BF, lwd=3, col="seagreen")
legend("bottomright",
     c(paste0("Marginal distribution by ", toupper(dtype)),
       "Upper bounds of fitted distribution",
       "Systematic peaks by Weibull plotting position"), bty="o", seg.len=3,
       pch=c(NA,NA,21), col=c("seagreen","purple","black"), bg="white", cex=0.8,
       lty=c(1, 2, NA), lwd=c(3, 1.5, NA), pt.bg=c(NA, NA, "white")) # }
}
\keyword{utility (probability)}
\keyword{conditional probability}
\keyword{probability (conditional)}
