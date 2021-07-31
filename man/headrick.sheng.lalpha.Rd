\name{headrick.sheng.lalpha}
\alias{headrick.sheng.lalpha}
\alias{lalpha}
\title{ The Sample Headrick and Sheng L-alpha }
\description{
Compute the sample Headrick and Sheng \dQuote{L-alpha} (Headrick and Sheng, 2013) by
\deqn{\alpha_L = \frac{d}{d-1}
   \biggl(1 - \frac{\sum_j \lambda^{(j)}_2}{\sum_j \lambda^{(j)}_2 + \sum\sum_{j\ne j'} \lambda_2^{(jj')}} \biggr)\mbox{,}}
where \eqn{j = 1,\ldots,d} for dimensions \eqn{d}, the \eqn{\sum_j \lambda^{(j)}_2} is the summation of all the 2nd order (univariate) L-moments (L-scales, \eqn{\lambda^{(j)}_2}), and the double summation is the summation of all the 2nd order L-comoments (\eqn{\lambda_2^{(jj')}}). In other words, the double summation is the sum of all entries in both the lower and upper triangles (not the primary diagonal) of the L-comoment matrix (the L-scale and L-coscale [L-covariance] matrix).
}
\usage{
headrick.sheng.lalpha(x, ...)

lalpha(x, ...)
}
\arguments{
  \item{x}{An \R \code{data.frame} of the random observations for the \eqn{d} random variables \eqn{X}, which must be suitable for internal dispatch to the \code{\link{Lcomoment.matrix}} function for the \code{k=2} order L-comoment. Alternatively, \code{x} can be a precomputed 2nd order L-comoment matrix (L-scale and L-coscale matrix).}
  \item{...}{Additional arguments to pass.}
}
\details{
Headrick and Sheng (2013) propose \eqn{\alpha_L} to be an alternative estimator of reliability based on L-comoments. They describe its context as follows: \dQuote{Consider [a statistic] alpha (\eqn{\alpha}) in terms of a model that decomposes an observed score into the sum of two independent components: a true unobservable score \eqn{t_i} and a random error component \eqn{\epsilon_{ij}}.} And the authors continue \dQuote{The model can be summarized as
\eqn{X_{ij} = t_i + \epsilon_{ij}\mbox{,}} where \eqn{X_{ij}} is the observed score associated with the \eqn{i}th examinee on the \eqn{j}th test item, and where \eqn{i = 1,...,n} [for sample size \eqn{n}]; \eqn{j = 1,\ldots,d}; and the error terms (\eqn{\epsilon_{ij}}) are independent with a mean of zero.} The authors go on to observe that \dQuote{inspection of [this model] indicates that this particular model restricts the true score \eqn{t_i} to be the same across all \eqn{d} test items.}

Headrick and Sheng (2013) show empirical results for a simulation study, which indicate that \eqn{\alpha_L} can be \dQuote{substantially superior} to [a different formulation of \eqn{\alpha} (Cronbach's Alpha) based on product moments (the variance-covariance matrix)] in \dQuote{terms of relative bias and relative standard error when distributions are heavy-tailed and sample sizes are small.}

The authors remind the reader that the second L-moments associated with \eqn{X_j} and \eqn{X_{j'}} can alternatively be expressed as
\eqn{\lambda_2(X_j) = 2\mathrm{Cov}(X_j,F(X_j))} and \eqn{\lambda_2(X_{j'}) = 2\mathrm{Cov}(X_{j'},F(X_{j'}))}. And that the second L-comoments of \eqn{X_j} toward (with respect to) \eqn{X_{j'}} and \eqn{X_{j'}} toward (with respect to) \eqn{X_j} are \eqn{\lambda_2^{(jj')} = 2\mathrm{Cov}(X_j,F(X_{j'}))} and \eqn{\lambda_2^{(j'j)} = 2\mathrm{Cov}(X_{j'},F(X_j))}. The respective cumulative distribution functions are denoted \eqn{F(x_j)}. Evidently the authors present the L-moments and L-comoments this way because their first example (thanks for detailed numerics!) already contain nonexceedance probabilities. Thus the function \code{headrick.sheng.lalpha} is prepared for two different contents of the \code{x} argument. One for a situation in which only the value for the random variables are available, and one for a situation in which the nonexceedances are already available. The numerically the two \eqn{\alpha_L} will not be identical as the example shows.
}
\value{
  An \R \code{list} is returned.

  \item{alpha}{The \eqn{\alpha_L} statistic.}
  \item{title}{The formal name \dQuote{Headrick and Sheng L-alpha}.}
  \item{source}{An attribute identifying the computational source of the Headrick and Sheng L-alpha: \dQuote{headrick.sheng.lalpha}.}
}
\note{
Headrick and Sheng (2013) use \eqn{k} to represent \eqn{d} as used here. The change is made because \code{k} is an L-comoment order argument already in use by \code{\link{Lcomoment.matrix}}.
}
\references{
Headrick, T.C., and Sheng, Y., 2013, An alternative to Cronbach's Alpha---An L-moment-based measure of internal-consistency reliability: \emph{in} Millsap, R.E., van der Ark, L.A., Bolt, D.M., Woods, C.M. (eds) New Developments in Quantitative Psychology, Springer Proceedings in Mathematics and Statistics, v. 66, \doi{10.1007/978-1-4614-9348-8_2}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{Lcomoment.matrix}}}
\examples{
# Table 1 in Headrick and Sheng (2013)
TV1 <- # Observations in cols 1:3, estimated nonexceedance probabilities in cols 4:6
c(2, 4, 3, 0.15, 0.45, 0.15,       5, 7, 7, 0.75, 0.95, 1.00,
  3, 5, 5, 0.35, 0.65, 0.40,       6, 6, 6, 0.90, 0.80, 0.75,
  7, 7, 6, 1.00, 0.95, 0.75,       5, 2, 6, 0.75, 0.10, 0.75,
  2, 3, 3, 0.15, 0.25, 0.15,       4, 3, 6, 0.55, 0.25, 0.75,
  3, 5, 5, 0.35, 0.65, 0.40,       4, 4, 5, 0.55, 0.45, 0.40)
T1 <- matrix(ncol=6, nrow=10)
for(r in seq(1,length(TV1), by=6)) T1[(r/6)+1, ] <- TV1[r:(r+5)]
colnames(T1) <- c("X1", "X2", "X3", "FX1", "FX2", "FX3"); T1 <- as.data.frame(T1)

lco2 <- matrix(nrow=3, ncol=3)
lco2[1,1] <- lmoms(T1$X1)$lambdas[2]
lco2[2,2] <- lmoms(T1$X2)$lambdas[2]
lco2[3,3] <- lmoms(T1$X3)$lambdas[2]
lco2[1,2] <- 2*cov(T1$X1, T1$FX2); lco2[1,3] <- 2*cov(T1$X1, T1$FX3)
lco2[2,1] <- 2*cov(T1$X2, T1$FX1); lco2[2,3] <- 2*cov(T1$X2, T1$FX3)
lco2[3,1] <- 2*cov(T1$X3, T1$FX1); lco2[3,2] <- 2*cov(T1$X3, T1$FX2)
headrick.sheng.lalpha(lco2)$alpha     # Headrick and Sheng (2013): alpha = 0.807
# 0.8074766
headrick.sheng.lalpha(T1[,1:3])$alpha #              FXs not used: alpha = 0.781
# 0.7805825
}
\keyword{L-alpha (sample)}
\keyword{Headrick and Sheng L-alpha}
\concept{internal consistency}
