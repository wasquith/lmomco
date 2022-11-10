\name{qua2ci.simple}
\alias{qua2ci}
\alias{qua2ci.simple}
\title{Estimate a Confidence Interval for a Single Quantile of a Parent Distribution by a Simple Algorithm}
\description{
This function estimates the lower and upper limits of a specified confidence interval for an aribitrary quantile value of a specified parent distribution [quantile function \eqn{Q(F,\theta)} with parameters \eqn{\theta}] using Monte Carlo simulation. The quantile value, actually the nonexceedance probability (\eqn{F} for \eqn{0 \le F \le 1}) of the value, is specified by the user. The user also provides the parameters of the parent distribution (see \code{\link{lmom2par}}). This function does consider an estimate of the variance-covariance structure of the sample data (for that see \code{\link{qua2ci.cov}}).  The \code{qua2ci.simple} is the original implementation and dates close to the initial releases of \pkg{lmomco} and was originally named \code{qua2ci}. That name is now deprecated but retained as an alias, which will be removed at some later release.

For \code{nsim} simulation runs (ideally a large number), samples of size \eqn{n} are drawn from \eqn{Q(F,\theta)}. The L-moments of each simulated sample are computed using \code{\link{lmoms}} and a distribution of the same type is fit. The \eqn{F}-quantile of the just-fitted distribution is computed and placed into a vector. The process of simulating the sample, computing the L-moments, computing the parameters, and solving for the \eqn{F}-quantile is repeated for the specified number of simulation runs.

To estimate the confidence interval, the L-moments of the vector simulated quantiles are computed. Subsequently, the parameters of a user-specified distribution \dQuote{error} distribution (\code{edist}) are computed. The two quantiles of this error distribution for the specified confidence interval are computed. These two quantiles represent the estimated lower and upper limits for the confidence interval of the parent distribution for samples of size \eqn{n}. The error distribution defaults to the Generalized Normal (see \code{\link{pargno}}) because this distribution has the Normal as a special case but extends the fit to the 3rd L-moment (\eqn{\tau_3}) for exotic situations in which some asymmetry in the quantile distribution might exist.

Finally, it is often useful to have vectors of lower and upper limits for confidence intervals for a vector of \eqn{F} values. The function \code{\link{genci.simple}} does just that and uses \code{\link{qua2ci.simple}} as the computational underpinning.
}
\usage{
qua2ci.simple(f,para,n, level=0.90, edist="gno", nsim=1000, showpar=FALSE,
                        empdist=TRUE, verbose=FALSE, maxlogdiff=6, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}) of the quantile for which the confidence interval is needed. This function is not vectorized and therefore only the first value will be used. This is in contrast to the vectorization of \eqn{F} in the conceptually similar function \code{\link{qua2ci.cov}}.}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}---these parameters represent the \dQuote{true} parent.}
  \item{n}{The sample size for each Monte Carlo simulation will use.}
  \item{level}{The confidence interval (\eqn{0 \le } \code{level} \eqn{ < 1}). The interval is specified as the size of the interval. The default is 0.90 or the 90th percentile. The function will return the 5th [\eqn{(1-0.90)/2}] and 95th [\eqn{(1-(1-0.90)/2)}] percentile cumulative probability of the simulated quantile distribution as specified by the nonexceedance probability argument. The arguments \code{level} and \code{f} therefore are separate features.}
  \item{edist}{The model for the error distribution. Although the Normal (the default) commonly is  assumed in error analyses, it need not be, as support for other distributions supported by \pkg{lmomco} is available. The default is the Generalized Normal so the not only is the Normal possible but asymmetry is also accomodated (\code{\link{lmomgno}}).  For example, if the L-skew  (\eqn{\tau_4}) or L-kurtosis (\eqn{\tau_4}) values depart considerably from those of the Normal (\eqn{\tau_3 = 0} and \eqn{\tau_4 = 0.122602}), then the Generalized Normal or some alternative distribution would likely provide more reliable confidence interval estimation.}
  \item{nsim}{The number of simulations (replications) for the sample size \code{n} to perform. Large numbers produce more refined confidence limit estimates at the cost of CPU time. The default is anticipated to be large enough for evaluative-useage without too much computational delay. Larger simulation numbers are recommended.}
  \item{showpar}{The parameters of the \code{edist} for each simulation are printed.}
  \item{empdist}{If \code{TRUE}, then an \R \code{environment} is appended onto the element \code{empdist} in the returned list, otherwise \code{empdist} is \code{NA}.}
  \item{verbose}{The verbosity of the operation of the function.}
  \item{maxlogdiff}{The maximum permitted difference in log10 space between a simulated quantile and the true value. It is possible that a well fit simulated sample to the parent distribution type provides crazy quantile estimates in the far reaches of either tail.  The default value of 6 was chosen based on experience with the Kappa distribution fit to a typical heavy-right tail flood magnitude data set. The concern motivating this feature is that as the number of parameters increases, it seems progressively there is more chance for a distribution tail to swing wildy into regions for which an analyst would not be comfortable with given discipline-specific knowledge. The choice of 6-log cycles is \emph{ad hoc} at best, and users are encouraged to do their own exploration. If \code{verbose=TRUE} then a message will be printed when the \code{maxlogdiff} condition is tripped.}
  \item{...}{Additional arguments to pass such as to \code{\link{lmom2par}}.}
}
\value{
  An \R \code{list} is returned. The \code{lwr} and \code{upr} match the nomenclature of \code{\link{qua2ci.cov}} but because \code{qua2ci.simple} is provided the parent, the \code{true} value is returned, whereas \code{\link{qua2ci.cov}} returns the \code{fit}.

  \item{lwr}{The lower value of the confidence interval having nonexceedance probability equal to \eqn{(1-}\code{level}\eqn{)/2}.}
  \item{true}{The value returned by \code{par2qua(f,para)}.}
  \item{upr}{The upper  value of the confidence interval having nonexceedance probability equal to \eqn{1-(1-}\code{level}\eqn{)/2}.}
  \item{elmoms}{The L-moments from \code{\link{lmoms}} of the distribution of simulated of quantiles.}
  \item{epara}{The parameters of the error distribution fit using the \code{elmoms}.}
  \item{empdist}{An \R \code{environment} (see below).}
  \item{ifail}{A diagnostic value. A value of zero means that successful exit was made.}
  \item{ifailtext}{A descriptive message related to the \code{ifail} value.}
  \item{nsim}{An echoing of the \code{nsim} argument for the function.}
  \item{sim.attempts}{The number of executions of the \code{while} loop (see Note below).}

The \code{empdist} element in the returned \code{list} is an \R \code{environment} that contains:

\item{simquas}{A \code{nsim}-long vector of the simulated quantiles for \code{f}.}
\item{empir.dist.lwr}{The \emph{lower} limit derived from the \R \code{quantile} function for \code{type=6}, which uses \eqn{i/(n+1)}.}
\item{empir.dist.upr}{The \emph{upper} limit derived from the \R \code{quantile} function for \code{type=6}, which uses \eqn{i/(n+1)}.}
\item{bern.smooth.lwr}{The \emph{lower} limit estimated by the Bernstein smoother in \code{\link{dat2bernqua}} for \cr \code{poly.type = "Bernstein"} and \code{bound.type = "none"}.}
\item{bern.smooth.upr}{The \emph{upper} limit estimated by the Bernstein smoother in \code{\link{dat2bernqua}} for \cr \code{poly.type = "Bernstein"} and \code{bound.type = "none"}.}
\item{epmoms}{The product moments of the simulated quantiles from \code{\link{pmoms}}.}
}
\note{
This function relies on a \code{while} loop that runs until \code{nsim} have successfully completed. Some reasons for an early \code{next} in the loop include invalid L-moments by \code{\link{are.lmom.valid}} of the simluated data or invalid fitted parameters by \code{\link{are.par.valid}} to simulated L-moments.  See the source code for more details.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmoms}}, \code{\link{pmoms}}, \code{\link{par2qua}}, \code{\link{genci.simple}}, \code{\link{qua2ci.cov}}
}
\examples{
\dontrun{
# It is well known that standard deviation (sigma) of the
# sample mean is equal to sigma/sample_size. Let is look at the
# quantile distribution of the median (f=0.5)
mean   <- 0; sigma <- 100
parent <- vec2par(c(mean,sigma), type='nor')
CI     <- qua2ci.simple(0.5, parent, n=10, nsim=20)
# Theoretrical sample mean sigma = 100/10 = 10
# L-moment theory: L-scale * sqrt(pi) = sigma
# Thus, it follows that the quantity
CI$elmoms$lambdas[2]/sqrt(pi)
# approaches 10 as nsim --> Inf.}

# Another example.
D   <- c(123, 34, 4, 654, 37, 78, 93, 95, 120) # fake sample
lmr <- lmoms(D)    # compute the L-moments of the data
WEI <- parwei(lmr) # estimate Weibull distribution parameters
CI  <- qua2ci.simple(0.75,WEI,20, nsim=20, level=0.95)
# CI contains the estimate 95-percent confidence interval for the
# 75th-percentile of the parent Weibull distribution for 20 sample size 20.
\dontrun{
pdf("Substantial_qua2ci_example.pdf")
level <- 0.90; cilo <- (1-level)/2; cihi <- 1 - cilo
para <- lmom2par(vec2lmom(c(180,50,0.75)), type="gev")
A <- qua2ci.simple(0.98, para, 30, edist="gno", level=level, nsim=3000)
Apara <- A$epara; Aenv <- A$empdist
Bpara <- lmom2par(A$elmoms, type="aep4")

lo <- log10(A$lwr); hi <- log10(A$upr)
xs <- 10^(seq(lo-0.2, hi+0.2, by=0.005))
lo <- A$lwr; hi <- A$upr; xm <- A$true; sbar <- mean(Aenv$simquas)
dd <- density(Aenv$simquas, adjust=0.5)
pk <- max(dd$y, dlmomco(xs, Apara), dlmomco(xs, Bpara))
dx <- dd$x[dd$x >= Aenv$empir.dist.lower & dd$x <= Aenv$empir.dist.upper]
dy <- dd$y[dd$x >= Aenv$empir.dist.lower & dd$x <= Aenv$empir.dist.upper]
dx <- c(dx[1], dx, dx[length(dx)]); dy <- c(0, dy, 0)

plot(c(0), c(0), type="n", xlim=range(xs), ylim=c(0,pk),
                 xlab="X VALUE", ylab="PROBABILITY DENSITY")
polygon(dx, dy, col=8)
lines(xs, dlmomco(xs, Apara)); lines(xs, dlmomco(xs, Bpara), col=2, lwd=2)
lines(dd, lty=2, lwd=2, col=8)
lines(xs, dlmomco(xs, para), col=6); lines(c(xm,xm), c(0,pk), lty=4, lwd=3)
lines(c(lo,lo,NA,hi,hi), c(0,pk,NA,0,pk), lty=2)

xlo <- qlmomco(cilo, Apara); xhi <- qlmomco(cihi, Apara)
points(c(xlo, xhi), c(dlmomco(xlo, Apara), dlmomco(xhi, Apara)), pch=16)
xlo <- qlmomco(cilo, Bpara); xhi <- qlmomco(cihi, Bpara)
points(c(xlo, xhi), c(dlmomco(xlo, Bpara), dlmomco(xhi, Bpara)), pch=16, col=2)
lines(rep(Aenv$empir.dist.lwr, 2), c(0,pk), lty=3, lwd=2, col=3)
lines(rep(Aenv$empir.dist.upr, 2), c(0,pk), lty=3, lwd=2, col=3)
lines(rep(Aenv$bern.smooth.lwr,2), c(0,pk), lty=3, lwd=2, col=4)
lines(rep(Aenv$bern.smooth.upr,2), c(0,pk), lty=3, lwd=2, col=4)
cat(c(  "F(true) = ",             round(plmomco(xm,   Apara), digits=2),
      "; F(mean(sim), edist) = ", round(plmomco(sbar, Apara), digits=2), "\n"), sep="")
dev.off()}
\dontrun{
ty <- "nor" # try running with "glo" (to get the L-skew "fit", see below)
para <- lmom2par(vec2lmom(c(-180,70,-.5)), type=ty)
f <- 0.99; n <- 41; ns <- 1000; Qtrue <- qlmomco(f, para)
Qsim1 <- replicate(ns, qlmomco(f, lmom2par(lmoms(rlmomco(n, para)), type=ty)))
Qsim2 <- qua2ci.simple(f, para, n, nsim=ns, edist="gno")
Qbar1 <- mean(Qsim1); Qbar2 <- mean(Qsim2$empdist$simquas)
epara <- Qsim2$epara; FT <- plmomco(Qtrue, epara)
F1 <- plmomco(Qbar1, epara); F2 <- plmomco(Qbar2, epara)
cat(c(  "F(true) = ",      round(FT, digits=2),
      "; F(via sim.) = ",  round(F1, digits=2),
      "; F(via edist) = ", round(F2, digits=2), "\n"), sep="")
# The given L-moments are highly skewed, but a Normal distribution is fit so
# L-skew is ignored. The game is deep tail (f=0.99) estimation. The true value of the
# quantile has a percentile on the error distribution 0.48 that is almost exactly 0.5
# (median = mean = symmetrical error distribution).  A test run shows nice behavior:
# F(true) =  0.48; F(via sim.) =  0.49; F(via edist) =  0.5
# But another run with ty <- "glo" (see how 0.36 << [0.52, 0.54]) has
# F(true) =  0.36; F(via sim.) =  0.54; F(via edist) =  0.52
# So as the asymmetry becomes extreme, the error distribution becomes asymmetrical too.}
}
\keyword{quantile function (confidence limit)}
