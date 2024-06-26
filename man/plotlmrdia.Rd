\name{plotlmrdia}
\alias{plotlmrdia}
\title{Plot L-moment Ratio Diagram (Tau3 and Tau4)}
\description{
Plot the Tau3-Tau4 L-moment ratio diagram of L-skew and L-kurtosis from a Tau3-Tau4 L-moment ratio diagram object returned by \code{\link{lmrdia}}. This diagram is useful for selecting a distribution to model the data. The application of L-moment diagrams is well documented in the literature. This function is intended to function as a demonstration of L-moment ratio diagram plotting with enough user settings for many practical applications.
}
\usage{
plotlmrdia(lmr=NULL, nopoints=FALSE, nolines=FALSE, nolimits=FALSE,
           noaep4=FALSE, nogev=FALSE, noglo=FALSE,  nogno=FALSE, nogov=FALSE,
           nogpa=FALSE,  nope3=FALSE, nopdq3=FALSE, nowei=TRUE,
           nocau=TRUE,   noexp=FALSE, nonor=FALSE,  nogum=FALSE,
           noray=FALSE, nosla=TRUE, nouni=FALSE,
           xlab="L-skew (Tau3), dimensionless",
           ylab="L-kurtosis (Tau4), dimensionless", add=FALSE, empty=FALSE,
           autolegend=FALSE, xleg=NULL, yleg=NULL, legendcex=0.9,
           ncol=1, text.width=NULL, lwd.cex=1, expand.names=FALSE, ...)
}
\arguments{
  \item{lmr}{L-moment diagram object from \code{\link{lmrdia}}, if \code{NULL}, then \code{empty} is internally set to \code{TRUE}.}
  \item{nopoints}{If \code{TRUE} then point distributions are not drawn.}
  \item{nolines}{If \code{TRUE} then line distributions are not drawn.}
  \item{nolimits}{If \code{TRUE} then theoretical limits of L-moments are not drawn.}
  \item{noaep4}{If \code{TRUE} then the lower bounds line of Asymmetric Exponential Power distribution is not drawn.}
  \item{nogev}{If \code{TRUE} then line of Generalized Extreme Value distribution is not drawn.}
  \item{noglo}{If \code{TRUE} then line of Generalized Logistic distribution is not drawn.}
  \item{nogno}{If \code{TRUE} then line of Generalized Normal (Log-Normal3) distribution is not drawn.}
  \item{nogov}{If \code{TRUE} then line of Govindarajulu distribution is not drawn.}
  \item{nogpa}{If \code{TRUE} then line of Generalized Pareto distribution is not drawn.}
  \item{nope3}{If \code{TRUE} then line of Pearson Type III distribution is not drawn.}
  \item{nopdq3}{If \code{TRUE} then line of Polynomial Density-Quantile3 distribution is not drawn.}
 \item{nowei}{If \code{TRUE} then line of the Weibull distribution is not drawn. The Weibull is a reverse of the Generalized Extreme Value. Traditionally in the literature, the Tau3-Tau4 L-moment ratio diagram have usually included the Weibull distribution and therefore the default setting of this argument is to not plot the Weibull.}
  \item{nocau}{If \code{TRUE} then point (TL-moment [trim=1]) of the Cauchy distribution is not drawn.}
  \item{noexp}{If \code{TRUE} then point of Exponential distribution is not drawn.}
  \item{nonor}{If \code{TRUE} then point of Normal distribution is not drawn.}
  \item{nogum}{If \code{TRUE} then point of Gumbel distribution is not drawn.}
  \item{noray}{If \code{TRUE} then point of Rayleigh distribution is not drawn.}
  \item{nouni}{If \code{TRUE} then point of Uniform distribution is not drawn.}
  \item{nosla}{If \code{TRUE} then point (TL-moment [trim=1]) of the Slash distribution is not drawn.}
  \item{xlab}{Horizonal axis label passed to \code{xlab} of the \code{plot} function.}
  \item{ylab}{Vertical axis label passed to \code{ylab} of the \code{plot} function.}
  \item{add}{A logical to toggle a call to \code{plot} to start a new plot, otherwise, just the trajectories are otherwise plotted.}
  \item{empty}{A logical to return before any trajectories are plotted but after the condition of the \code{add} has been evaluated.}
  \item{autolegend}{Generate the legend by built-in algorithm.}
  \item{xleg}{X-coordinate of the legend. This argument is checked for being a character versus a numeric. If it is a character, then \code{yleg} is not needed and \code{xleg} and take on \dQuote{location may also be specified by setting x to a single keyword} as per the functionality of \code{graphics::legend()} itself.}
  \item{yleg}{Y-coordinate of the legend.}
  \item{legendcex}{The \code{cex} to pass to \code{graphics::legend()}.}
    \item{ncol}{The number of columns in which to set the legend items (default is 1, which differs from \code{legend()} default of 1).}
  \item{text.width}{Argument of the same name for \code{legend}. Setting to 0.1 for \code{ncol} set to 2 seems to work pretty well when two columns are desired.}
  \item{lwd.cex}{Expansion factor on the line widths.}
  \item{expand.names}{Expand the distribution names in the legend.}
  \item{...}{Additional arguments passed into \code{plot()} and \code{legend()} functions..}
}
\note{
This function provides hardwired calls to \code{lines} and \code{points} to produce the diagram. The plot symbology for the shown distributions is summarized here. The Asymmetric Exponential Power and Kappa (four parameter) and Wakeby (five parameter) distributions are not well represented on the diagram as each constitute an area (Kappa) or hyperplane (Wakeby) and not a line (3-parameter distributions) or a point (2-parameter distributions). However, the Kappa demarks the area bounded by the Generalized Logistic (\code{glo}) on the top and the
theoretical L-moment limits on the bottom. The Asymmetric Exponential Power demarks its own unique lower boundary and extends up in the \eqn{\tau_4} direction to \eqn{\tau_4 = 1}. However, parameter estimation with L-moments has lost considerable accuracy for \eqn{\tau_4} that large (see Asquith, 2014).

\tabular{cl}{
  \bold{GRAPHIC TYPE} \tab \bold{GRAPHIC NATURE} \cr
  L-moment Limits \tab line width 2 and color a medium-dark grey \cr
  Asymmetric Exponential Power (4-p) \tab line width 1, line type 4 (dot), and color red \cr
  Generalized Extreme Value (GEV) \tab line width 1, line type 1 (solid), and color darkred \cr
  Generalized Logistic \tab line width 1 and color green \cr
  Generalized Normal \tab line width 1, line type 2 (dash), and color blue \cr
  Govindarajulu \tab line width 1, line type 2 (dash), and color 6 (magenta) \cr
  Generalized Pareto \tab line width 1, line type 1 (solid), and color blue \cr
  Pearson Type III \tab line width 1, line type 1 (solid), and color 6 (purple) \cr
  Polynomial Density-Quantile3 \tab line width 1.3, line type 2 (dash), and color darkgreen \cr
  Weibull (reversed GEV) \tab line width 1, line type 1 (solid), and color darkorange \cr
  Exponential \tab symbol 16 (filled circle) and color red \cr
  Normal \tab symbol 15 (filled square) and color red \cr
  Gumbel \tab symbol 17 (filled triangle) and color red) \cr
  Rayleigh \tab symbol 18 (filled diamond) and color red \cr
  Uniform \tab symbol 12 (square and a plus sign) and color red \cr
  Cauchy \tab symbol 13 (circle with over lapping \eqn{\times}) and color turquoise4 \cr
  Slash  \tab symbol 10 (cicle containing \eqn{+}) and color turquoise4 \cr
}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.

Hosking, J.R.M., 1986, The theory of probability weighted moments: Research Report RC12210, IBM Research Division, Yorkton Heights, N.Y.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis--An approach based on L-moments: Cambridge University Press.

Vogel, R.M., and Fennessey, N.M., 1993, L moment diagrams should replace product moment diagrams: Water Resources Research, v. 29, no. 6, pp. 1745--1752.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmrdia}}, \code{\link{plotlmrdia46}}, \code{\link{plotradarlmr}}}
\examples{
plotlmrdia(lmrdia()) # simplest of all uses

\dontrun{
# A more complex example follows: for a given mean, L-scale, L-skew, and L-kurtosis,
# use sample size of 30, use 500 simulations, set L-moments, fit the Kappa distribution
T3 <- 0.34; T4 <- 0.21; n <- 30; nsim <- 500
lmr <- vec2lmom(c(10000, 7500, T3, T4)); kap <- parkap(lmr)

# create vectors for storing simulated L-skew (t3) and L-kurtosis (t4)
t3 <- t4 <- vector(mode="numeric")

# perform nsim simulations by randomly drawing from the Kappa distribution
# and compute the L-moments in sim.lmr and store the t3 and t4 of each sample
for(i in 1:nsim) {
  sim.lmr <- lmoms(rlmomco(n, kap))
  t3[i] <- sim.lmr$ratios[3]; t4[i] <- sim.lmr$ratios[4]
}

# plot the diagram and "zoom" by manually setting the axis limits
plotlmrdia(xlim=c(-0.1, 0.5), ylim=c(-0.1, 0.4), las=1, empty=TRUE)

# Follow up by plotting the {t3, t4} values and the mean of the values
points(t3, t4, pch=21, bg="white", lwd=0.8) # plot each simulation

# plot crossing dashed lines at true values of L-skew and L-kurtosis
abline(v=T3, col="salmon4", lty=2, lwd=3) # Theoretical values for the
abline(h=T4, col="salmon4", lty=2, lwd=3) # distribution as fit

points(mean(t3), mean(t4), pch=16, cex=3) # mean of simulations and
# should plot reasonably close to the salmon4-colored crossing lines

# plot the trajectories of the distributions
plotlmrdia(lmrdia(), add=TRUE, nopoints=TRUE, inset=0.01,
           autolegend=TRUE, xleg="topleft", lwd.cex=1.5) # }
}
\keyword{L-moment ratio diagram}
\keyword{Graphics}
