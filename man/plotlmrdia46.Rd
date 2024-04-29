\name{plotlmrdia46}
\alias{plotlmrdia46}
\title{Plot L-moment Ratio Diagram (Tau4 and Tau6)}
\description{
Plot the Tau4-Tau6 L-moment ratio diagram showing trajectories of \eqn{\tau_4} and \eqn{\tau_6} for strictly symmetrical distributions from a Tau4-Tau6 L-moment ratio diagram object returned by \code{\link{lmrdia46}}. This diagram is useful for selecting among symmetrical distributions to model the data. This function is intended to function as a demonstration of Tau4-Tau6 L-moment ratio diagram plotting with enough user settings for many practical applications.
}
\usage{
plotlmrdia46(lmr=NULL, nopoints=FALSE, nolines=FALSE,
             noaep4=FALSE, nogld_byt5opt=TRUE,
             nopdq4=FALSE, nost3=FALSE, nosymstable=FALSE, notukey=FALSE,
             nocau=TRUE,   nonor=FALSE, nosla=TRUE, trucate.tau4.to.gtzero=TRUE,
             xlab="L-kurtosis (Tau4), dimensionless",
             ylab="Sixth L-moment ratio (Tau6), dimensionless",
             add=FALSE, empty=FALSE,
             autolegend=FALSE, xleg=NULL, yleg=NULL, legendcex=0.9,
             ncol=1, text.width=NULL, lwd.cex=1, expand.names=FALSE, ...)
}
\arguments{
  \item{lmr}{L-moment diagram object from \code{\link{lmrdia46}}, if \code{NULL}, then \code{empty} is internally set to \code{TRUE}.}
  \item{nopoints}{If \code{TRUE} then point distributions are not drawn.}
  \item{nolines}{If \code{TRUE} then line distributions are not drawn.}
  \item{noaep4}{If \code{TRUE} then the Symmetric Exponential Power distribution is not drawn.}
  \item{nogld_byt5opt}{If \code{TRUE} then line of Generalized Lambda distribution through it solution optimization on \eqn{\tau_5 = 0} is not drawn.}
  \item{nopdq4}{If \code{TRUE} then line of Polynomial Density-Quantile4 distribution is not drawn.}
  \item{nost3}{If \code{TRUE} then line of Student 3t distribution is not drawn.}
  \item{nosymstable}{If \code{TRUE} then line of Symmetric Stable distribution is not drawn.}
  \item{notukey}{If \code{TRUE} then line of Tukey Lambda distribution is not drawn.}
  \item{nocau}{If \code{TRUE} then point of Cauchy distribution (trim=1 L-moments) is not drawn.}
  \item{nonor}{If \code{TRUE} then point of Normal distribution is not drawn.}
  \item{nosla}{If \code{TRUE} then point of Slash distribution (trim=1 L-moments) is not drawn.}
  \item{trucate.tau4.to.gtzero}{Truncate the distributions that can extend to negative \eqn{\tau_4} to zero. This is a reasonable default and prevents line drawing to the left into a clipping region for easier handling of post processing of a graphic in vector editing software.}
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
  \item{...}{Additional arguments passed into the \code{plot()} and \code{legend()} functions.}
}
\note{
This function provides hardwired calls to \code{lines} and \code{points} to produce the diagram. The plot symbology for the shown distributions is summarized here.

\tabular{cl}{
  \bold{GRAPHIC TYPE} \tab \bold{GRAPHIC NATURE} \cr
  Symmetric Exponential Power \tab line width 1, line type 4 (dot), and color red \cr
  Generalized Lambda \tab line width 1, line type 1 (solid), and color purple \cr
  Polynomial Density-Quantile4 \tab line width 1, line type 1 (solid), and color darkgreen \cr
  Student t \tab line width 1, line type 1 (solid), and color blue\cr
  Symmetric Stable \tab line width 2, line type 1 (solid), and color a medium-dark grey\cr
  Tukey Lambda (1-p) \tab line width 1, line type 2 (dash), and color purple\cr
  Normal \tab symbol 15 (filled square) and color red \cr
  Cauchy \tab symbol 13 (circle with over lapping \eqn{\times}) and color turquoise4 \cr
  Slash  \tab symbol 10 (cicle containing \eqn{+}) and color turquoise4 \cr

}
}
\references{
Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmrdia46}}, \code{\link{plotlmrdia}}}
\examples{
plotlmrdia46(lmrdia46(), nogld_byt5opt=FALSE, autolegend=TRUE,
             xleg="topleft")


\dontrun{
# A more complex example follows: for a given mean, L-scale, L-skew = 0 (symmetry), and
# L-kurtosis, use sample size of 30, use 500 simulations, set L-moments,
# fit the Asymmetric Exponential Power4 distribution, which is symmetrical when the
# L-skew is zero and thus the distribution is the Exponential Power.
T3  <- 0; T4 <- 0.21; n <- 30; nsim <- 500
lmr <- vec2lmom(c(10000, 7500, T3, T4, 0)); aep4 <- paraep4(lmr)
T6  <- theoLmoms(aep4, nmom=6)$ratios[6]

# create vectors for storing simulated L-kurtosis (t4) and Tau6 (t6)
t4 <- t6 <- vector(mode="numeric")

# perform nsim simulations by randomly drawing from the AEP4 distribution
# and compute the L-moments in sim.lmr and store the t4 and t6 of each sample
for(i in 1:nsim) {
  sim.lmr <- lmoms(rlmomco(n, aep4), nmom=6)
  t4[i] <- sim.lmr$ratios[4]; t6[i] <- sim.lmr$ratios[6]
}

# plot the diagram and "zoom" by manually setting the axis limits
plotlmrdia46(xlim=c(-0.05, 0.5), ylim=c(-0.1, 0.35), las=1, empty=TRUE)

# follow up by plotting the {t3, t4} values and the mean of the values
points(t4, t6, cex=0.8, pch=21, bg="white", lwd=0.8) # plot each simulation

# plot crossing dashed lines at true values of L-skew and L-kurtosis
abline(v=T4, col="salmon4", lty=2, lwd=3) # Theoretical values for the
abline(h=T6, col="salmon4", lty=2, lwd=3) # distribution as fit

points(mean(t4), mean(t6), pch=16, cex=3) # mean of simulations and
# should plot reasonably close to the salmon4-colored crossing lines

# plot the trajectories of the distributions
plotlmrdia46(lmrdia46(), add=TRUE, nopoints=TRUE, inset=0.01,
             autolegend=TRUE, xleg="topleft", lwd.cex=1.5) # }
}
\keyword{L-moment ratio diagram}
\keyword{L-moment ratio diagram (Tau4 and Tau6)}
\keyword{Distribution: Tukey Lambda}
\keyword{Distribution: Stable}
\keyword{Graphics}
