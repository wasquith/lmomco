\encoding{utf8}
\name{headrick.sheng.lalpha}
\alias{headrick.sheng.lalpha}
\alias{lalpha}
\title{Sample Headrick and Sheng L-alpha}
\description{
Compute the sample \dQuote{Headrick and Sheng L-alpha} (\eqn{\alpha_L}) (Headrick and Sheng, 2013) by
\deqn{\alpha_L = \frac{d}{d-1}
   \biggl(1 - \frac{\sum_j \lambda^{(j)}_2}{\sum_j \lambda^{(j)}_2 + \sum\sum_{j\ne j'} \lambda_2^{(jj')}} \biggr)\mathrm{,}}
where \eqn{j = 1,\ldots,d} for dimensions \eqn{d}, the \eqn{\sum_j \lambda^{(j)}_2} is the summation of all the 2nd order (univariate) L-moments (L-scales, \eqn{\lambda^{(j)}_2}), and the double summation is the summation of all the 2nd-order L-comoments (\eqn{\lambda_2^{(jj')}}). In other words, the double summation is the sum total of all entries in both the lower and upper triangles (not the primary diagonal) of the L-comoment matrix (the L-scale and L-coscale [L-covariance] matrix) (\code{\link{Lcomoment.matrix}}).

The \eqn{\alpha_L} is closely related in structural computation as the well-known \dQuote{Cronbach alpha} (\eqn{\alpha_C}). These are coefficients of reliability, which commonly ranges from 0 to 1, that provide what some methodologists portray as an overall assessment of a measure's reliability. If all of the scale items are entirely independent from one another, meaning that they are not correlated or share no covariance, then \eqn{\alpha_C} is 0, and, if all of the items have high covariances, then \eqn{\alpha_C} will approach 1 as the number of items in the scale approaches infinity. The higher the \eqn{\alpha_C} coefficient, the more the items have shared covariance and probably measure the same underlying concept. Theoretically, there is no lower bounds for \eqn{\alpha_{C,L}}, which can add complicating nuances in bootstrap or simulation study of both \eqn{\alpha_C} and \eqn{\alpha_L}. Negative values are considered a sign of something potentially wrong about the measure related to items not being positively correlated with each other, or a scoring system for a question item reversed. (This paragraph in part paraphrases \url{https://data.library.virginia.edu/using-and-interpreting-cronbachs-alpha/}\cr (accessed May 21, 2023) and other general sources.)
}
\usage{
headrick.sheng.lalpha(x, bycovFF=FALSE, a=0.5, digits=8, ...)

lalpha(x, bycovFF=FALSE, a=0.5, digits=8, ...)
}
\arguments{
  \item{x}{An \R \code{data.frame} of the random observations for the \eqn{d} random variables \eqn{X}, which must be suitable for internal dispatch to the \code{\link{Lcomoment.matrix}} function for computation of the 2nd-order L-comoment matrix. Alternatively, \code{x} can be a precomputed 2nd-order L-comoment matrix (L-scale and L-coscale matrix) as shown by the following usage: \code{lalpha(Lcomoment.matrix(x, k=2)$matrix).}}
  \item{bycovFF}{A logical triggering the covariance pathway for the computation and bypassing the call to the L-comoments. The additional arguments can be used to control the \code{\link{pp}} function that is called internally to estimate nonexceedance probabilities and the \dQuote{covariance pathway} (see \bold{Details}). If \code{bycovFF} is \code{FALSE}, then the direct to L-comoment computation is used.}
  \item{a}{The plotting position argument \code{a} to the \code{\link{pp}} function that is hardwired here to Hazen in contrast to the default \code{a=0} of \code{\link{pp}} (Weibull) for reasoning shown in this documentation.}
  \item{digits}{Number of digits for rounding on the returned value(s).}
  \item{...}{Additional arguments to pass.}
}
\details{
Headrick and Sheng (2013) propose \eqn{\alpha_L} to be an alternative estimator of reliability based on L-comoments. Those authors describe its context as follows: \dQuote{Consider [a statistic] alpha (\eqn{\alpha}) in terms of a model that decomposes an observed score into the sum of two independent components: a true unobservable score \eqn{t_i} and a random error component \eqn{\epsilon_{ij}}.}

Those authors continue \dQuote{The model can be summarized as
\eqn{X_{ij} = t_i + \epsilon_{ij}\mathrm{,}} where \eqn{X_{ij}} is the observed score associated with the \eqn{i}th examinee on the \eqn{j}th test item, and where \eqn{i = 1,...,n} [for sample size \eqn{n}]; \eqn{j = 1,\ldots,d}; and the error terms (\eqn{\epsilon_{ij}}) are independent with a mean of zero.} Those authors comment that \dQuote{inspection of [this model] indicates that this particular model restricts the true score \eqn{t_i} to be the same across all \eqn{d} test items.}

Those authors show empirical results for a simulation study, which indicate that \eqn{\alpha_L} can be \dQuote{substantially superior} to [a different formulation of \eqn{\alpha_C} (Cronbach's alpha) based on product moments (the variance-covariance matrix)] in \dQuote{terms of relative bias and relative standard error when distributions are heavy-tailed and sample sizes are small.}

Those authors show (Headrick and Sheng, 2013, eqs. 4 and 5) the reader that the second L-comoments of \eqn{X_j} and \eqn{X_{j'}} can alternatively be expressed as
\eqn{\lambda_2(X_j) = 2\mathrm{Cov}(X_j, F(X_j))} and \eqn{\lambda_2(X_{j'}) = 2\mathrm{Cov}(X_{j'}, F(X_{j'}))}. The second L-comoments of \eqn{X_j} toward (with respect to) \eqn{X_{j'}} and \eqn{X_{j'}} toward (with respect to) \eqn{X_j} are \eqn{\lambda_2^{(jj')} = 2\mathrm{Cov}(X_j, F(X_{j'}))} and \eqn{\lambda_2^{(j'j)} = 2\mathrm{Cov}(X_{j'}, F(X_j))}. The respective cumulative distribution functions are denoted \eqn{F(x_j)} (nonexceedance probabilities). Evidently, those authors present the L-moments and L-comoments this way because their first example (thanks for detailed numerics!) already contain nonexceedance probabilities.

This apparent numerical difference between the version using estimates of nonexceedance probabilities for the data (the \dQuote{covariance pathway}) compared to a \dQuote{direct to L-comoment} pathway might be more than academic concern.

The \bold{Examples} provide comparison and brief discussion of potential issues involved in the direct L-comoments and the covariance pathway. The discussion leads to interest in the effects of ties and their handling and the question of \eqn{F(x_j)} estimation by plotting position (\code{\link{pp}}). The \bold{Note} section of this documentation provides expanded information and insights to \eqn{\alpha_L} computation.
}
\value{
  An \R \code{list} is returned.

  \item{alpha}{The \eqn{\alpha_L} statistic.}
  \item{pitems}{The number of items (column count) in the \code{x}.}
  \item{n}{The sample size (row count), if applicable, to the contents of \code{x}.}
  \item{text}{Any pertinent messages about the computations.}
  \item{source}{An attribute identifying the computational source of the Headrick and Sheng L-alpha: \dQuote{headrick.sheng.lalpha} or \dQuote{lalpha.star()}.}
}
\note{
Headrick and Sheng (2013) use \eqn{k} to represent \eqn{d} as used here. The change is made because \code{k} is an L-comoment order argument already in use by \code{\link{Lcomoment.matrix}}.


\bold{Demonstration of Nuances of L-alpha}---Consider Headrick and Sheng (2013, tables 1 and 2) and the effect of those authors' covariance pathway to \eqn{\alpha_L}:
\preformatted{
  X1 <- c(2, 5, 3, 6, 7, 5, 2, 4, 3, 4) # Table 1 in Headrick and Sheng (2013)
  X2 <- c(4, 7, 5, 6, 7, 2, 3, 3, 5, 4)
  X3 <- c(3, 7, 5, 6, 6, 6, 3, 6, 5, 5)
  X  <- data.frame(X1=X1, X2=X2, X3=X3)
  lcm2 <- Lcomoment.matrix(X, k=2)
  print(lcm2$matrix, 3)
  #       [,1]  [,2]  [,3]
  # [1,] 0.989 0.567 0.722
  # [2,] 0.444 1.022 0.222
  # [3,] 0.644 0.378 0.733
}

Now, compare the above matrix to Headrick and Sheng (2013, table 2) where it is immediately seen that the matrices are not the same before the summations are applied to compute \eqn{\alpha_L}.
\preformatted{
  #       [,1]  [,2]  [,3]
  # [1,] 0.989 0.500 0.789
  # [2,] 0.500 1.022 0.411
  # [3,] 0.667 0.333 0.733
}

Now, consider how the nonexceedances in Headrick and Sheng (2013, table 1) might have been computed w/o following their citation to original sources. It can be shown with reference to the first example above that these nonexceedance probabilities match.
\preformatted{
  FX1 <- rank(X$X1, ties.method="average") / length(X$X1)
  FX2 <- rank(X$X2, ties.method="average") / length(X$X2)
  FX3 <- rank(X$X3, ties.method="average") / length(X$X3)
}
Notice in Headrick and Sheng (2013, table 1) that there is no zero probability, but there is a unity and some of the probabilities are tied. Ties have numerical ramifications. Let us now look at other L-alphas using the nonexceedance pathway and use different definitions of nonexceedance estimation and inspect the results:
\preformatted{
  # lmomco documentation says pp() uses ties.method="first"
  lalpha(X, bycovFF=TRUE, a=0     )$alpha
  # [1] 0.7448583  # unbiased probs all distributions
  lalpha(X, bycovFF=TRUE, a=0.3173)$alpha
  # [1] 0.7671384  # Median probs for all distributions
  lalpha(X, bycovFF=TRUE, a=0.35  )$alpha
  # [1] 0.7695105  # Often used with probs-weighted moments
  lalpha(X, bycovFF=TRUE, a=0.375 )$alpha
  # [1] 0.771334   # Blom, nearly unbiased quantiles for normal
  lalpha(X, bycovFF=TRUE, a=0.40  )$alpha
  # [1] 0.7731661  # Cunnane, appox quantile unbiased
  lalpha(X, bycovFF=TRUE, a=0.44  )$alpha
  # [1] 0.7761157  # Gringorten, optimized for Gumbel
  lalpha(X, bycovFF=TRUE, a=0.5   )$alpha
  # [1] 0.7805825  # Hazen, traditional choice
                   # This the plotting position (i-0.5) / n
}

This is not a particularly pleasing situation because the choice of the plotting position affects the \eqn{\alpha_L}. The Hazen definition \code{lalpha(X[,1:3], bycovFF=FALSE)} using direct to L-comoments matches the last computation shown (\eqn{\alpha_L = 0.7805825}). A question, thus, is does this matching occur because of the nature of the ties and structure of the L-comoment algorithm itself? A note to this question involves a recognition that the answer is yes because L-comoments use a \code{sort()} operation and does not use \code{rank()} because the weights for the linear combinations are used and the covariance pathway \code{2*cov(x$X3, x$FX2)}, for instance.

Recognizing that the direct to L-comoments alpha equals the covariance pathway with Hazen plotting positions, let us look at L-comoments:
\preformatted{
  lmomco::Lcomoment.Lk12 ------> snippet
       X12 <- X1[sort(X2, decreasing = FALSE, index.return = TRUE)$ix]
       n <- length(X1)
       SUM <- sum(sapply(1:n, function(r) { Lcomoment.Wk(k, r, n) * X12[r] }))
       return(SUM/n)
}
Notice that a \code{ties.method} is not present but kind of implicit as ties first by the index return of the \code{sort()} and notice the return of a \code{SUM/n} though this is an L-comoment and not an nonexceedance probability.

Let us run through the tie options using a plotting position definition (\eqn{i / n}) matching the computations of Headrick and Sheng (2013) (\code{"average"}, \code{A=0}, \code{B=0} for \code{\link{pp}}) and the first computation \eqn{\alpha_L = 0.807} matches that reported by Headrick and Sheng (2013, p. 4):
\preformatted{
  for(tie in c("average", "first", "last", "min", "max")) { # random left out
    Lalpha <- lalpha(X, bycovFF=TRUE,
                        a=NULL, A=0, B=0, ties.method=tie)$alpha
    message("Ties method ", stringr::str_pad(tie, 7, pad=" "),
            " provides L-alpha = ", Lalpha)
  }
  # Ties method average provides L-alpha = 0.80747664
  # Ties method   first provides L-alpha = 0.78058252
  # Ties method    last provides L-alpha = 0.83243243
  # Ties method     min provides L-alpha = 0.81363468
  # Ties method     max provides L-alpha = 0.80120709
}

Let us run through the tie options again using a different plotting position estimator (\eqn{(n-0.5) / (n+0.5)}):
\preformatted{
  for(tie in c("average", "first", "last", "min", "max")) { # random left out
    Lalpha <- lalpha(X, bycovFF=TRUE,
                        a=NULL, A=-0.5, B=0.5, ties.method=tie)$alpha
    message("Ties method ", stringr::str_pad(tie, 7, pad=" "),
            " provides L-alpha = ", Lalpha)
  }
  # Ties method average provides L-alpha = 0.78925733
  # Ties method   first provides L-alpha = 0.76230208
  # Ties method    last provides L-alpha = 0.81431215
  # Ties method     min provides L-alpha = 0.79543602
  # Ties method     max provides L-alpha = 0.78296931
}

We see obviously that the decision on how to treat ties has some influence on the computation involving the covariance pathway. This is not an entirely satisfactory situation, but perhaps the distinction does not matter much? The direct L-comoment pathway seems to avoid this issue because \code{sort()} is stable and like \code{ties.method="first"}. Experiments suggest that \code{a=0.5} (Hazen plotting positions) produces the same results as direct L-comoment (see the next section). However, as the following code set shows:
\preformatted{
  for(tie in c("average", "first", "last", "min", "max")) { # random left out
    Lalpha1 <- lalpha(X, bycovFF=TRUE, a=0.5, ties.method=tie)$alpha
    Lalpha2 <- lalpha(X, bycovFF=TRUE, a=NULL, A=-0.5, B=0, ties.method=tie)$alpha
    Lalpha3 <- lalpha(X, bycovFF=TRUE, a=NULL, A=-1  , B=0, ties.method=tie)$alpha
    Lalpha4 <- lalpha(X, bycovFF=TRUE, a=NULL, A=   0, B=0, ties.method=tie)$alpha
    print(c(Lalpha1, Lalpha2, Lalpha3, Lalpha4))
  }
}
The \eqn{\alpha_L} for a given tie setting are all equal as long as the demoninator of the plotting position (\eqn{(i + A) / (n + B)}) has \code{B=0}. The \code{a=0.5} produces Hazen, the \code{a=NULL, A=-0.5} produces Hazen, though \code{a=NULL, A=-1} (lower limit of \code{A}) and \code{a=NULL, A=0} (upper limit of \code{A} given \code{B}) also produces the same. This gives us as-implemented-proof that the sensitivity to the \eqn{\alpha_L} computation is in the sorting and the denominator of the plotting position formula. The prudent default settings for when the \code{bycovFF} argument is true seems to have the \eqn{a=-0.5} as nonexceedance probabilities are computed by the well-known Hazen description and with the tie method as first, the computations match direct to L-comoments.

\bold{Demonstration of Computational Times}---A considerable amount of documentation and examples are provided here about the two pathways that \eqn{\alpha_L} can be computed: (1) direct by L-comoments or (2) covariance pathway requiring precomputed estimates of the nonexceedance probabilities using a \code{ties.method="first"} (default \code{\link{pp}}). The following example shows numerical congruence between the two pathways if the so-called Hazen plotting positions (\code{a=0.5}, see \code{\link{pp}}) are requested with the implicit default of \code{ties.method="first"}. However, the computational time of the direct method is quite a bit longer because of latencies in the weight factor computations involved in the L-comoments and nested \code{for} loops.
\preformatted{
  set.seed(1)
  R <- 1:10; nsam <- 1E5 # random and uncorrelated scores in this measure
  Z <- data.frame( I1=sample(R, nsam, replace=TRUE),
                   I2=sample(R, nsam, replace=TRUE),
                   I3=sample(R, nsam, replace=TRUE),
                   I4=sample(R, nsam, replace=TRUE) )
  system.time(AnF <- headrick.sheng.lalpha(Z, bycovFF=FALSE)$alpha)
  system.time(AwF <- headrick.sheng.lalpha(Z, bycovFF=TRUE )$alpha) # Hazen
  #    user  system elapsed
  #  30.382   0.095  30.501    AnF ---> 0.01370302
  #    user  system elapsed
  #   5.054   0.030   5.092    AwF ---> 0.01370302
}
}
\references{
Headrick, T.C., and Sheng, Y., 2013, An alternative to Cronbach's alpha---An L-moment-based measure of internal-consistency reliability: \emph{in} Millsap, R.E., van der Ark, L.A., Bolt, D.M., Woods, C.M. (eds) New Developments in Quantitative Psychology, Springer Proceedings in Mathematics and Statistics, v. 66, \doi{10.1007/978-1-4614-9348-8_2}.
}
\author{W.H. Asquith}
\seealso{\code{\link{Lcomoment.matrix}}, \code{\link{pp}}}
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
headrick.sheng.lalpha(Lcomoment.matrix(T1[,1:3], k=2)$matrix)$alpha
# 0.7805825
headrick.sheng.lalpha(T1[,1:3])$alpha #              FXs not used: alpha = 0.781
# 0.7805825
headrick.sheng.lalpha(T1[,1:3], bycovFF=TRUE)$alpha  # a=0.5, Hazen by default
# 0.7805825
headrick.sheng.lalpha(T1[,1:3], bycovFF=TRUE, a=0.5)$alpha
# 0.7805825
}
\keyword{L-alpha (coefficient of reliability)}
\keyword{Headrick and Sheng L-alpha}
\keyword{L-alpha (Headrick and Sheng)}
\keyword{Cronbach alpha}
\keyword{alpha (Cronbach)}
\keyword{L-comoment}
\keyword{coefficient of reliability}
\concept{internal consistency}
