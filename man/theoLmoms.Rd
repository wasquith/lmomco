\name{theoLmoms}
\alias{theoLmoms}
\title{The Theoretical L-moments and L-moment Ratios using Integration of the Quantile Function}
\description{
Compute the theoretrical L-moments for a vector. A theoretrical L-moment in integral form is
\deqn{ \lambda_r = \frac{1}{r}
                               \sum^{r-1}_{k=0}{(-1)^k {r-1 \choose k}
        \frac{r!\:I_r}{(r-k-1)!\,k!}
       } \mbox{,}}
in which
\deqn{ I_r = \int^1_0 x(F) \times F^{r-k-1}(1-F)^{k}\,\mathrm{d}F \mbox{,}}
where \eqn{x(F)} is the quantile function of the random variable \eqn{X} for nonexceedance probability \eqn{F}, and \eqn{r} represents the order of the L-moments. This function actually dispatches to \code{\link{theoTLmoms}} with \code{trim=0} argument.
}
\usage{
theoLmoms(para, nmom=5, minF=0, maxF=1, quafunc=NULL,
                nsim=50000, fold=5,
                silent=TRUE, verbose=FALSE, ...)
}
\arguments{
  \item{para}{A distribution parameter object such as from \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{minF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-zero (but very small)  if the integral is divergent.}
  \item{maxF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-unity (but still very close [perhaps \code{1 - minF}]) if the integral is divergent.}
  \item{quafunc}{An optional and arbitrary quantile function that simply needs to except a nonexceedance probability and the parameter object in \code{para}. This is a feature that permits computation of the L-moments of a quantile function that does not have to be implemented in the greater overhead hassles of the \pkg{lmomco} style. This feature might be useful for estimation of quantile function mixtures or those distributions not otherwise implemented in this package.}
  \item{nsim}{Simulation size for Monte Carlo integration is such is internally deemed necessary (see \code{silent} argument).}
   \item{fold}{The number of fractions or number of folds of \code{nsim}, which in other words, means that \code{nsim} is divided by \code{folds} and a loop creating \code{folds} integrations of \code{nsim/folds} is used from which the mean and mean absolute error of the integrand are computed. This is to try to recover similar output as \code{integrate()}.}
   \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}. If set true and the integral is probability divergent, Monte Carlo integration is triggered using \code{nsim} and \code{folds}. The user would have to set \code{verbose=TRUE} to then acquire the returned table in \code{integration_table} of the integration passes including those are or are not Monte Carlo.}
   \item{verbose}{Toggle verbose output. Because the \R function \code{integrate} is used to perform the numerical integration, it might be useful to see selected messages regarding the numerical integration.}
   \item{...}{Additional arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the TL-moments. First element is \eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is \eqn{\tau_2}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal zero (the ordinary L-moments) because this function dispatches to \code{\link{theoTLmoms}}.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{theoLmoms}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{W.H. Asquith}
\seealso{\code{\link{theoTLmoms}}}
\examples{
para <- vec2par(c(0,1), type='nor') # standard normal
TL00 <- theoLmoms(para) # compute ordinary L-moments
}
\keyword{L-moment (theoretical)}

