\encoding{utf8}
\name{disfitgovloc}
\alias{disfitgovloc}
\title{Fit a Govindarajulu Distribution to Bounds and Location }
\description{
Fits a \emph{Govindarajulu} distribution to specified lower and upper bounds and a given location measure (either mean and median). Fitting occurs through \eqn{3}-dimensional minimization using the \code{optim} function. Objective function forms are either root mean-square error (RMSE) or mean absolute deviation (MAD), and the objective functions are expected to result in slightly different estimates of distribution parameters. The RMSE form (\eqn{\sigma_{\mathrm{RMSE}}}) is defined as

\deqn{\sigma_{\mathrm{RMSE}} = \biggl[ \frac{1}{3}\,\sum_{i=1}^3 \bigl[x_i - \hat{x}_i\bigr]^2\biggr]^{1/2}\mbox{,}}

where \eqn{x_i} is a vector of the targeted lower bounds (\code{lwr} argument), location measure (\code{loc} argument), and upper bounds (\code{upr} argument), and \eqn{\hat{x}_i} is a similar vector of Govindarajulu properties for \dQuote{current} iteration of the optimization. Similarly, the MAD form (\eqn{\sigma_{\mathrm{MAD}}}) is defined as

\deqn{\sigma_{\mathrm{MAD}} = \frac{1}{3}\,\sum_{i=1}^3 \mid x_i - \hat{x}_i \mid  \mbox{.}}

The premise of this function is that situations might exist in practical applications wherein the user has an understanding or commitment to certain bounding conditions of a distribution. The user also has knowledge of a particular location measure (the mean or median) of a distribution. The bounded nature of the Govindarajulu might be particularly of interest because the quantile function (\code{\link{quagov}}) is explicit. The curvatures that the distribution can attain also provide it more flexibility to fitting to a given location measure than say the \emph{Triangular} distribution (\code{\link{quatri}}).
}
\usage{
disfitgovloc(x=NULL, loc=NULL, lwr=0, upr=NA, init.para=NULL,
             loctype=c("mean", "median"), objfun=c("rmse", "mad"),
             ptransf=function(p) return(log(p)),
             pretransf=function(p) return(exp(p)),
             silent=TRUE, verbose=FALSE, ...)
}
\arguments{
  \item{x}{Optional vector to help guide the initial parameter estimates for the optimization, if given and if \code{loc=NULL}, then \code{loc} by \code{loctype} will be computed from the \code{x}.}
  \item{loc}{Optional value for the location statistic, which if not given will be computed from mean or median of the \code{x}. The \code{loc} however can also be given if an \code{x} is given and at which point the user's setting prevails.}
  \item{lwr}{Lower bounds for the distribution with default supposing that most often positive domain bounds might be of interest.}
  \item{upr}{Upper bounds for the distribution, which must be specified.}
  \item{init.para}{Optional initial values for the parameters used for starting values for the \code{optim} function. If this argument is not set nor is \code{x}, then an unrigorous attempt is made to guess at the initial parameters using heuristics and the triangular quantile function (because the triangle is trivial and also bounded) (see sources).}
  \item{loctype}{The type of location measure constraint.}
  \item{objfun}{The form of the objective function as previously described.}
  \item{ptransf}{The parameter transformation function that is useful to guide the optimization run. The distribution requires its second and third parameters to be nonzero without constraint on the first parameter; however, the default treats the first parameter as also nonzero. This is potentially suboptimal for some situations (see \bold{Examples}).}
  \item{pretransf}{The parameter retransformation function that is useful to guide the optimization run. The distribution requires its second and third parameters to be nonzero without constraint on the first parameter; however, the default treats the first parameter as also nonzero. This is potentially suboptimal for some situations (see \bold{Examples}).}
  \item{silent}{A logical to silence the \code{try()} function wrapping the \code{optim()} function.}
  \item{verbose}{A logical to trigger verbose output within the objective function.}
  \item{...}{Additional arguments to pass to the \code{optim} function.}
}
\details{
Support of the Govindarajulu for the optimized parameter set is computed by internally and reported as part of the returned values. This enhances the documentation a bit more---the computed parameters might not always have full convergence and result in slightly difference bounds than targeted. Finally, this function was developed using some heredity to \code{\link{disfitqua}}.
}
\value{
  An \R \code{list} is returned.  This list should contain at least the following items.
  \item{type}{The type of distribution in three character (minimum) format.}
  \item{para}{The parameters of the Govindarajulu distribution.}
  \item{source}{Attribute specifying source of the parameters.}
  \item{supdist}{A list of confirming the distribution support from \code{quagov(c(0,1), gov)} where \code{gov} are the final computed parameters before return.}
  \item{init.para}{A vector of the initial parameters actually passed to the \code{optim} function to serve only as a reminder.}
  \item{optim}{The returned \code{list} of the \code{optim()} function.}
  \item{message}{Helpful messages on the computations.}
}
\author{W.H. Asquith}
\seealso{\code{\link{disfitqua}}, \code{\link{quagov}}}
\examples{
# EXAMPLE 1 --- Example of strictly positive domain.
disfitgovloc(loc=125, lwr=99, upr=175, loctype="mean")$para
#        xi     alpha      beta
# 99.000000 76.000000  3.846154
# These parameters have a lmomgov()$lambdas[1] mean of 124.9999999.

# EXAMPLE 2 --- Operations spanning zero and revision to the default parameter
# transform functions. Testing indicates that these, ideally align to need of
# the Govindarajulu, such do not work for all strictly positive domain, which
# led to a decision to have the defaults different than this example.
disfitgovloc(loc=100, lwr=-99, upr=175, loctype="median",
               ptransf=function(p) c(p[1], log(p[2:3])),
             pretransf=function(p) c(p[1], exp(p[2:3])))$para
#         xi        alpha         beta
# -99.000002   274.000004   1.08815151

\dontrun{
  # EXTENDED EXAMPLE 3
  r <- function(r) round(r, 1)
  X <- c(8751, 14507, 4061, 22056, 6330, 3130, 5180, 6700, 22409, 3380, 17902,
         8956,  4523, 1604,  4460, 4239, 3010, 9155, 5107, 4821,  5221, 20700)
  mu  <-   mean(X); med <- median(X)
  for(objfun in c("rmse", "mad")) {
    gov <- disfitgovloc(x=X,  loc=mu,  upr=41000, objfun=objfun, loctype="mean"    )
    message(objfun, ": seek   mean=", r(mu),
                    ", GOV   mean=",  r(lmomgov(gov)$lambdas[1]))
    gov <- disfitgovloc(x=X, loc=med,  upr=41000, objfun=objfun, loctype="median"  )
    message(objfun, ": seek median=", r(med),
                    ", GOV median=",  r(quagov(0.5, gov)))
  }
  for(objfun in c("rmse", "mad")) {
    gov <- disfitgovloc(x=NULL,  loc=mu,  upr=41000, objfun=objfun, loctype="mean"  )
    message(objfun, ": seek   mean=", r(mu),
                    ", GOV   mean=",  r(lmomgov(gov)$lambdas[1]) )
    gov <- disfitgovloc(x=NULL, loc=med,  upr=41000, objfun=objfun, loctype="median")
    message(objfun, ": seek median=", r(med),
                    ", GOV median=",  r(quagov(0.5, gov)))
  } # end of loop
  # *** That last message() : mad: seek median=5200.5, GOV median=5226.2
  print(gov$para) # 64.521326, 40935.479117, 4.740232 # last parameters in prior loop
  ngv <- vec2par( c(64.521326, 40935.479117, 4.740232), type="gov") # for reuse
  # We see (at least in testing) that the last message in the sequence shows that
  # the median is not recovered via the guessed at initial parameters, let us turn
  # the gov parameters back into disfitgovloc() as the initial parameters.
  mgv <- disfitgovloc(init.para=ngv, loc=med, upr=41000, objfun=objfun,loctype="median")
  message(objfun, ": seek median=", r(med),
                   ", GOV median=", r(quagov(0.5, mgv)))
  # *** BETTER FIT mad: seek median=5200.5, GOV median=5200.5
  print(mgv$para) # 1.227568, 40998.903644, 4.729768 # last parameters
  # So, conveniently in this example, we can see that there are cases wherein an
  # apparent convergence can be made even better. But, need to be aware that
  # feed fack a very good solution can in turn cause optim() itself to NULL out. }

\dontrun{
  # EXTENDED EXAMPLE 4 --- Continuing from the previous example
  FF    <- seq(0.001, 0.999, by=0.001)
  maxes <- as.integer(10^(seq(4, 5, by=0.02))); n <- length(maxes)
  for(max in maxes) {
    govA <- disfitgovloc(x=X,  loc=mu,     upr=max, loctype="mean"  , lwr=0)
    govB <- disfitgovloc(x=X,  loc=median, upr=max, loctype="median", lwr=0)
    plot( FF, quagov(FF, govA), col="red",  lwd=2, type="l", ylim=c(0, maxes[n]),
         xlab="Nonexceedance probability", ylab="Quantile of Govindarajulu",
         main=paste0("Maximum = ", max))
    lines(FF, quagov(FF, govB), col="blue", lwd=2); quagov(0.5, govB)
    legend("topleft", c("Govindarajulu constrained given mean (dashed red)",
                        "Govindarajulu constrained given median (dashed blue)",
                        "disfitgovloc() computed mean (red dot)",
                        "disfitgovloc() computed median (blue dot)"),
                    lwd=c( 2,  2, NA, NA), col=c("red", "blue"), inset=0.02,
                    pch=c(NA, NA, 16, 16), pt.cex=1.5, cex=0.9)
    abline(h=mu,  lty=2, col="red" ); abline(h=med, lty=2, col="blue")
    tmu <- lmomgov(govA)$lambdas[1]
    points(cdfgov(tmu, govA), tmu, cex=1.5, pch=16, col="red" )
    points(0.5, quagov(0.5, govB), cex=1.5, pch=16, col="blue")
  } # end of loop }
}
\keyword{distribution (parameters)}
\keyword{Distribution: Govindarajulu}
