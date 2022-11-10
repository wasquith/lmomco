\name{lmrdiscord}
\alias{lmrdiscord}
\title{Compute Discordance on L-CV, L-skew, and L-kurtosis }
\description{
This function computes the Hosking and Wallis discordancy of the first three L-moment ratios (L-CV, L-skew, and L-kurtosis) according to their implementation in Hosking and Wallis (1997) and earlier. Discordancy triplets of these L-moment ratios is heuristically measured by effectively locating the triplet from the mean center of the 3-dimensional cloud of values. The \pkg{lmomRFA} provides for discordancy embedded in the \dQuote{L-moment method} of regional frequency analysis. The author of \pkg{lmomco} chooses to have a separate \dQuote{high level} implementation for emergent ideas of his in evaluating unusual sample distributions outside of the \code{regdata} object class envisioned by Hosking in the \pkg{lmomRFA} package.

Let \eqn{\bm{\mu_i}} be a row vector of the values of \eqn{\tau^{[i]}_2, \tau^{[i]}_3, \tau^{[i]}_4} and these are the L-moment ratios for the \eqn{i}th group or site out of \eqn{n} sites. Let \eqn{\bm{\overline\mu}} be a row vector of mean values of all the \eqn{n} sites. Defining a sum of squares and cross products  \eqn{3\times 3} matrix as
\deqn{\bm{S} = \sum_i^n (\bm{\mu} - \bm{\overline\mu})(\bm{\mu} - \bm{\overline\mu})^{T}}
compute the discorancy of the \eqn{i}th site as
\deqn{
D_i = \frac{n}{3} (\bm{\mu} - \bm{\overline\mu})^T \bm{S}^{-1} (\bm{\mu} - \bm{\overline\mu}\mbox{.})
}

The L-moments of a sample for a location are judged to be discordance if \eqn{D_i} exceeds a critical value. The critical value is a function of sample size. Hosking and Wallis (1997, p. 47) provide a table for general application. By about \eqn{n=14}, the critical value is taken as \eqn{D_c = 3}, although the \eqn{D_{max}} increases with sample size. Specifically, the \eqn{D_i} has an upper limit of
\deqn{D_i \le (n-1)/3\mbox{.}}
However, Hosking and Wallis (1997, p. 47) recommend \dQuote{that any site with \eqn{D_i > 3} be regarded as discordant.} A statistical test of \eqn{D_i} can be constructed. Hosking and Wallis (1997, p. 47) report that the \eqn{D_{critical}} is
\deqn{D_{critical, n, \alpha} = \frac{(n - 1)Z}{n - 4 + 3Z}\mbox{,}}
where
\deqn{Z = F(\alpha/n, 3, n - 4)\mbox{,}}
upper-tail quantile of the F distribution with degrees of freedom 3 and \eqn{n - 4}. A table of critical values is preloaded into the \code{\link{lmrdiscord}} function as this mimics the table of Hosking and Wallis (1997, table 3.1) as a means for cross verification. This table corresponds to an \eqn{\alpha = 0.1} significance.
}
\usage{
lmrdiscord(site=NULL, t2=NULL, t3=NULL, t4=NULL,
           Dcrit=NULL, digits=4, lmrdigits=4, sort=TRUE,
           alpha1=0.10, alpha2=0.01, ...)
}
\arguments{
  \item{site}{An optional group or site identification; it will be sequenced from 1 to \eqn{n} if \code{NULL}.}
  \item{t2}{L-CV values; emphasis that L-scale is not used.}
  \item{t3}{L-skew values.}
  \item{t4}{L-kurtosis values.}
  \item{Dcrit}{An optional (user specified) critical value for discordance. This value will override the Hosking and Wallis (1997, table 3.1) critical values.}
  \item{digits}{The number of digits in rounding operations.}
  \item{lmrdigits}{The numer of digits in rounding operation for the echo of the L-moment ratios.}
  \item{sort}{A logical on the sort status of the returned data frame.}
  \item{alpha1}{A significance level that is greater (less significant, although in statistics we need to avoid assigning less or more in this context) than \code{alpha2}.}
  \item{alpha2}{A significance level that is less (more significant, although in statistics we need to avoid assigning less or more in this context) than \code{alpha1}.}
  \item{...}{Other arguments that might be used. The author added these because it was found that the function was often called by higher level functions that aggregated much of the discordance computations.}
}
\value{
  An \R \code{data.frame} is returned.
  \item{site}{The group or site identification as used by the function.}
  \item{t2}{L-CV values.}
  \item{t3}{L-skew values.}
  \item{t4}{L-kurtosis.}
  \item{Dmax}{The maximum discordancy \eqn{D_{max} = (n-1)/3}.}
  \item{Dalpha1}{The critical value of \eqn{D} for \eqn{\alpha_1 = 0.10} (default) significance as set by \code{alpha1} argument.}
  \item{Dalpha2}{The critical value of \eqn{D} for \eqn{\alpha_2 = 0.01} (default) significance as set by \code{alpha1} argument.}
  \item{Dcrit}{The critical value of discordancy (user or tabled).}
  \item{D}{The discordancy of the L-moment ratios used to trigger the logical in \code{isD}.}
  \item{isD}{Are the L-moment ratios discordant (if starred).}
  \item{signif}{A hyphen, star, or double star based on the \code{Dalpha1} and \code{Dalpha2} values.}
}
\source{
Consultation of the \code{lmomRFA.f} and \code{regtst()} function of the \pkg{lmomRFA} \R package by J.R.M. Hosking. Thanks Jon and Jim Wallis for such a long advocation of the discordancy issue that began at least as early as the 1993 Water Resources Research Paper (-wha).
}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith }
\seealso{\code{\link{lmoms}}
}
\examples{
\dontrun{
# This is the canonical test of lmrdiscord().
library(lmomRFA) # Import lmomRFA, needs lmom package too
data(Cascades)   # Extract Hosking's data use in his examples
data <- as.regdata(Cascades) # A "regional" data structure
Dhosking <- sort(regtst(data)$D, decreasing=TRUE) # Discordancy

Dlmomco <- lmrdiscord(site=data$name, t2=data$t, t3=data$t_3, t4=data$t_4)

Dasquith <- Dlmomco$D
# Now show the site id, and the two discordancy computations
print(data.frame(NAME=data$name, Dhosking=Dhosking,
                                 Dasquith=Dasquith))
# The Dhosking and Dasquith columns had better match!

set.seed(3) # This seed produces a "*" and "**", but users
# are strongly encouraged to repeat the folowing code block
# over and over with an unspecified seed and look at the table.
n <- 30 # simulation sample size
par1 <- lmom2par(vec2lmom(c(1, .23, .2, .1)), type="kap")
par2 <- lmom2par(vec2lmom(c(1, .5, -.1)),      type="gev")
name <- t2 <- t3 <- t4 <- vector(mode="numeric")
for(i in 1:20) {
  X <- rlmomco(n, par1); lmr <- lmoms(X)
  t2[i] <- lmr$ratios[2]
  t3[i] <- lmr$ratios[3]
  t4[i] <- lmr$ratios[4]
  name[i] <- "kappa"
}
j <- length(t2)
for(i in 1:3) {
  X <- rlmomco(n, par2); lmr <- lmoms(X)
  t2[j + i] <- lmr$ratios[2]
  t3[j + i] <- lmr$ratios[3]
  t4[j + i] <- lmr$ratios[4]
  name[j + i] <- "gev"
}
D <- lmrdiscord(site=name, t2=t2, t3=t3, t4=t4)
print(D)

plotlmrdia(lmrdia(), xlim=c(-.2,.6), ylim=c(-.1, .4),
           autolegend=TRUE, xleg=0.1, yleg=.4)
points(D$t3,D$t4)
text(D$t3,D$t4,D$site, cex=0.75, pos=3)
text(D$t3,D$t4,D$D, cex=0.75, pos=1) #}
}
\keyword{L-moment (hypothesis)}

