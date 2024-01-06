\encoding{utf8}
\name{disfitqua}
\alias{disfitqua}
\title{Fit a Distribution using Minimization of Available Quantiles }
\description{
This function fits a distribution to available quantiles (or irregular quantiles) through \eqn{n}-dimensional minimization using the \code{optim} function. Objective function forms are either root mean-square error (RMSE) or mean absolute deviation (MAD), and the objective functions are expected to result in slightly different estimates of distribution parameters. The RMSE form (\eqn{\sigma_{\mathrm{RMSE}}}) is defined as
\deqn{\sigma_{\mathrm{RMSE}} = \biggl[ \frac{1}{m}\,\sum_{i=1}^m \bigl[x_o(f_i) - \hat{x}(f_i)\bigr]^2\biggr]^{1/2}\mbox{,}}
where \eqn{m} is the length of the vector of \eqn{o}bserved quantiles \eqn{x_o(f_i)} for nonexceedance probability \eqn{f_i} for \eqn{i \in 1, 2, \cdots, m}, and \eqn{\hat{x}(f_i)} for \eqn{i \in 1, 2, \cdots, m} are quantile estimates based on the \dQuote{current} iteration of the parameters for the selected distribution having \eqn{n} parameters for \eqn{n \le m}. Similarly, the MAD form (\eqn{\sigma_{\mathrm{MAD}}}) is defined as
\deqn{\sigma_{\mathrm{MAD}} = \frac{1}{m}\,\sum_{i=1}^m \mid x_o(f_i) - \hat{x}(f_i) \mid \mbox{.}}

The \code{disfitqua} function is not intended to be an implementation of the \emph{method of percentiles} but rather is intended for circumstances in which the available quantiles are restricted to either the left or right tails of the distribution. It is evident that a form of the method of percentiles however could be pursued by \code{disfitqua} when the length of \eqn{x(f)} is equal to the number of distribution parameters (\eqn{n = m}).  The situation of \eqn{n < m} however is thought to be the most common application.

The right-tail restriction is the general case in flood-peak hydrology in which the median and select quantiles greater than the median can be available from empirical studies (e.g. Asquith and Roussel, 2009) or rainfall-runoff models. The available quantiles suit engineering needs and thus left-tail quantiles simply are not available. This circumstance might appear quite unusual to users from most statistical disciplines but quantile estimates can exist from regional study of observed data. The \bold{Examples} section provides further motivation and discussion.
}
\usage{
disfitqua(x, f, objfun=c("rmse", "mad"),
                init.lmr=NULL, init.para=NULL, type=NA,
                ptransf=  function(t) return(t),
                pretransf=function(t) return(t), verbose=FALSE, ... )
}
\arguments{
  \item{x}{The quantiles \eqn{x_o(f)} for the nonexceedance probabilities in \code{f}.}
  \item{f}{The nonexceedance probabilities \eqn{f} of the quantiles \eqn{x_o(f)} in \code{x}.}
  \item{objfun}{The form of the objective function as previously described.}
  \item{init.lmr}{Optional initial values for the L-moments from which the initial starting parameters for the optimization will be determined. The optimizations by this function are not performed on the L-moments during the optimization. The form of \code{init.lmr} is that of an L-moment object from the \pkg{lmomco} package (e.g. \code{\link{lmoms}}).}
  \item{init.para}{Optional initial values for the parameters used for starting values for the \code{optim} function. If this argument is not set nor is \code{init.lmr}, then unrigorous estimates of the mean \eqn{\lambda_1} and L-scale \eqn{\lambda_2} are made from the available quantiles, higher L-moment ratios \eqn{\tau_r} for \eqn{r \ge 3} are set to zero, and the L-moments converted to the initial parameters.}
  \item{type}{The distribution type specified by the abbreviations listed under \code{\link{dist.list}}.}
  \item{ptransf}{An optional parameter transformation function (see \bold{Examples}) that is useful to guide the optimization run. For example, suppose the first parameter of a three parameter distribution resides in the positive domain, then \cr \code{ptransf(t) = }
\code{function(t) c(log(t[1]), t[2], t[3])}.}
  \item{pretransf}{An optional parameter retransformation function (see \bold{Examples}) that is useful to guide the optimization run. For example, suppose the first parameter of a three parameter distribution resides in the positive domain, then \cr \code{pretransf(t) = } \code{function(t) c(exp(t[1]), t[2], t[3])}.}
  \item{verbose}{A logical switch on the verbosity of output.}
  \item{...}{Additional arguments to pass to the \code{optim} function.}
}
\value{
  An \R \code{list} is returned, and this \code{list} contains at least the following items:

  \item{type}{The type of distribution in character format (see \code{\link{dist.list}}).}
  \item{para}{The parameters of the distribution.}
  \item{source}{Attribute specifying source of the parameters---\dQuote{disfitqua}.}
  \item{init.para}{A vector of the initial parameters actually passed to the \code{optim} function to serve only as a reminder.}
  \item{disfitqua}{The returned \code{list} from the \code{optim} function. This \code{list} contains a repeat of the parameters, the value of the objective function (\eqn{\sigma_{\mathrm{RMSE}}} or \eqn{\sigma_{\mathrm{MAD}}}), the interation count, and convergence status.}
}
\note{
The \code{disfitqua} function is likely more difficult to apply for \eqn{n > 3} (high parameter) distributions because of the inherent complexity of the mathematics of such distributions and their applicable parameter (and thus valid L-moment ranges). The complex interplay between parameters and L-moments can make identification of suitable initial parameters \code{init.para} or initial L-moments \code{init.lmr} more difficult than is the case for \eqn{n \le 3} distributions. The default initial parameters are computed from an assumed condition that the L-moments ratios \eqn{\tau_r = 0} for \eqn{r \ge 3}. This is not ideal, however, and the \bold{Examples} show how to move into high parameter distributions using the results from a previous fit.
}
\references{
Asquith, W.H., and Roussel, M.C., 2009, Regression equations for estimation of annual peak-streamflow frequency for undeveloped watersheds in Texas using an L-moment-based, PRESS-minimized, residual-adjusted approach: U.S. Geological Survey Scientific Investigations Report 2009--5087, 48 p., \doi{10.3133/sir20095087}.
}
\author{W.H. Asquith}
\seealso{\code{\link{dist.list}}, \code{\link{lmoms}}, \code{\link{lmom2vec}}, \code{\link{par2lmom}}, \code{\link{par2qua}}, \code{\link{vec2lmom}}, \code{\link{vec2par}} }
\examples{
# Suppose the following quantiles are estimated using eight equations provided by
# Asquith and Roussel (2009) for some watershed in Texas:
Q <- c(1480, 3230, 4670, 6750, 8700, 11000, 13600, 17500)
# These are real estimates from a suite of watershed properties; the watershed
# itself and location are not germane to demonstrate this function.
LQ <- log10(Q) # transform to logarithms of cubic feet per second
# Convert the averge annual return periods for the quantiles into probability
P <- T2prob(c(2, 5, 10, 25, 50, 100, 200, 500)); qP <- qnorm(P) # std norm variates
# The log-Pearson type III (LPIII) is immensely popular for flood-risk computations.
# Let us compute LPIII parameters to the available quantiles and probabilities for
# the watershed. The log-Pearson type III is "pe3" in the lmomco with logarithms.
par1 <- disfitqua(LQ, P, type="pe3", objfun="rmse") # root mean square error
par2 <- disfitqua(LQ, P, type="pe3", objfun="mad" ) # mean absolute deviation
# Now express the fitted distributions in forms of an LPIII.
LQfit1 <- qlmomco(P, par1); LQfit2 <- qlmomco(P, par2)

plot( qP, LQ, pch=5, xlab="STANDARD NORMAL VARIATES",
                     ylab="FLOOD QUANTILES, CUBIC FEET PER SECOND")
lines(qP, LQfit1, col=2); lines(qP, LQfit2, col=4) # red and blue lines

\dontrun{
# Now demonstrate how a Wakeby distribution can be fit. This is an example of how a
# three parameter distribution might be fit, and then the general L-moments secured for
# an alternative fit using a far more complicated distribution. The Wakeby for the
# above situation does not fit out of the box.
lmr1 <- theoLmoms(par1) # We need five L-moments but lmompe3() only gives four,
# therefore must compute the L-moment by numerical integration provided by theoLmoms().
par3 <- disfitqua(LQ, P, type="wak", objfun="rmse", init.lmr=lmr1)
lines(qP, par2qua(P, par3), col=6, lty=2) # dashed line, par2qua alternative to qlmomco

# Finally, the initial L-moment equivalents and then the L-moments of the fitted
# distribution can be computed and compared.
par2lmom(vec2par(par3$init.para, type="wak"))$ratios # initial L-moments
par2lmom(vec2par(par3$para,      type="wak"))$ratios # final   L-moments}
}
\keyword{distribution (parameters)}
\keyword{quantile function}
\keyword{method of percentiles}
\concept{irregular quantiles}
\concept{available quantiles}
\keyword{Ensemble list of major parameter estimation methods}

