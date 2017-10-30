\name{lmoms.bootbarvar}
\alias{lmoms.bootbarvar}
\title{ The Exact Bootstrap Mean and Variance of L-moments }
\description{
This function computes the exact bootstrap mean and variance of L-moments using the exact analytical expressions for the bootstrap mean and variance of any L-estimator described Hutson and Ernst (2000). The approach by those authors is to use the bootstrap distribution of the single order statistic in conjunction with the joint distribution of two order  statistics. The key component is the bootstrap mean vector as well as the variance-covariance matrix of all the order statistics and then performing specific linear combinations of a basic L-estimator combined with the proportion weights used in the computation of L-moments (\code{\link{Lcomoment.Wk}}, see those examples and division by \eqn{n}). Reasonably complex algorithms are used; however, what makes those authors' contribution so interesting is that neither simulation, resampling, or numerical methods are needed as long as the sample size is not too large.

This function provides a uniquely independent method to compute the L-moments of a sample from the vector of exact bootstrap order statistics. It is anticipated that several of the intermediate computations of this function would be of interest in further computations or graphical visualization. Therefore, this function returns many more numerical values than other L-moment functions of \pkg{lmomco}. The variance-covariance matrix for large samples requires considerable CPU time; as the matrix is filled, status output is generated.

 The example section of this function contains the verification of the implementation as well as provides to additional computations of variance through resampling with replacement and simulation from the parent distribution that generated the sample vector shown in the example.
}
\usage{
lmoms.bootbarvar(x, nmom=6, covarinverse=TRUE, verbose=TRUE,
                    force.exact=FALSE, nohatSIGMA=FALSE, nsim=500, bign=40, ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{nmom}{The number of moments to compute. Default is 6 and can not be less than 3.}
  \item{covarinverse}{Logical on computation of the matrix inversions: \cr \code{inverse.varcovar.tau23}, \cr \code{inverse.varcovar.tau34}, and \cr \code{inverse.varcovar.tau46}.}
  \item{verbose}{A logical switch on the verbosity of the construction of the variance-covariance matrix of the order statisitics. This operation is the most time consuming of those inside the function and is provided at default of \code{verbose=TRUE} to make a general user comfortable.}
  \item{force.exact}{A logical switch to attempt a \emph{forced exact bootstrap computation} (empirical bootstrap controlled by \code{nsim} thus is \emph{not} used) even if the sample size is too large as controlled by \code{bign}. See messages during the execution for guidance.}
  \item{nohatSIGMA}{A logical to bypass most of the interesting matrix functions and results. If \code{TRUE}, then only \code{lambdas}, \code{ratios}, and \code{bootstrap.orderstatistics} are populated. This feature is useful if a user is only interested in get the bootstrap estimates of the order statistics.}
  \item{nsim}{Simulation size in case simulations and not the exact bootstrap are used.}
  \item{bign}{A sample size threshold that triggers simulation using \code{nsim} replications for estimation by empirical bootstrap. Some of the \dQuote{exact} operations are extremely expensive and numerical problems in the matrices are known for non-normal data.}
  \item{...}{Additional arguments but not implemented.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the exact bootstrap L-moments. First element is
\eqn{\hat{\lambda}_1}, second element is \eqn{\hat{\lambda}_2}, and so on. This vector is from equation 1.3 and 2.4 of Hutson and Ernst (2000).}
  \item{ratios}{Vector of the exact bootstrap L-moment ratios. Second element is
\eqn{\hat{\tau}}, third element is \eqn{\hat{\tau}_3} and so on.}
  \item{lambdavars}{The exact bootstrap variances of the L-moments from equation 1.4 of Hutson and Ernst (2000) via \code{crossprod} matrix operations.}
  \item{ratiovars}{The exact bootstrap variances of the L-moment ratios with \code{NA} inserted for \eqn{r=1,2} because \eqn{r=1} is the mean and \eqn{r=2} for L-CV is unknown to this author.}
  \item{varcovar.lambdas}{The variance-covariance matrix of the L-moments from which the diagonal are the values \code{lambdavars}.}
  \item{varcovar.lambdas.and.ratios}{The variance-covariance matrix of the first two L-moments and for the L-moment ratios (if \code{nmom}\eqn{>=3}) from which select diagonal are the values \code{ratiovars}.}
  \item{bootstrap.orderstatistics}{The exact bootstrap estimate of the order statistics from equation 2.2 of Hutson and Ernst (2000).}
  \item{varcovar.orderstatistics}{The variance-covariance matrix of the order statistics from equations 3.1 and 3.2 of Hutson and Ernst (2000). The diagonal of this matrix represents the variances of each order statistic.}
  \item{inverse.varcovar.tau23}{The inversion of the variance-covariance matrix of \eqn{\tau_2} and \eqn{\tau_3} by Cholesky decomposition. This matrix may be used to estimate a joint confidence region of (\eqn{\tau_2, \tau_3}) based on asymptotic normality of L-moments.}
  \item{inverse.varcovar.tau34}{The inversion of the variance-covariance matrix of \eqn{\tau_3} and \eqn{\tau_4} by Cholesky decomposition. This matrix may be used to estimate a joint confidence region of (\eqn{\tau_3, \tau_4}) based on asymptotic normality of L-moments; these two L-moment ratios likely represent the most common ratios used in general L-moment ratio diagrams.}
  \item{inverse.varcovar.tau46}{The inversion of the variance-covariance matrix of \eqn{\tau_4} and \eqn{\tau_6} by Cholesky decomposition. This matrix may be used to estimate a joint confidence region of (\eqn{\tau_4, \tau_6}) based on asymptotic normality of L-moments; these two L-moment ratios represent those ratios used in L-moment ratio diagrams of symmetrical distributions.}
  \item{source}{An attribute identifying the computational source of the results: \cr \dQuote{lmoms.bootbarvar}.}
}
\note{
This function internally defines several functions that provide a direct nomenclature connection to Hutson and Ernst (2000). Interested users are invited to adapt these functions as they might see fit. A reminder is made to sort the data vector as needed; the vector is only sorted once within the \code{lmoms.bootbarvar} function.

The \eqn{100(1-\alpha)} percent confidence region of the vector \eqn{{\bm \eta} = (\tau_3, \tau_4)} (for example) based on the sample L-skew and L-kurtosis of the vector \eqn{\hat{\bm \eta} = (\hat\tau_3, \hat\tau_4)} is expressed as

\deqn{({\bm \eta} - \hat{\bm \eta})'\hat{\bm P}^{-1}_{(3,4)}({\bm \eta} - \hat{\bm \eta}) \le \chi^2_{2,\alpha}}

where \eqn{\hat{\bm P}_{(3,4)}} is the variance-covariance matrix of these L-moment ratios subselected from the resulting matrix titled \code{varcovar.lambdas.and.ratios} but extracted and inverted in the resulting matrix titled \code{inverse.varcovar.tau34}, which is \eqn{\hat{\bm P}^{-1}_{(3,4)}}. The value  \eqn{\chi^2_{2,\alpha}} is the upper quantile of the Chi-squared distribution. The inequality represents a standard equal probable ellipse from a Bivariate Normal distribution.
}
\references{
Hutson, A.D., and Ernst, M.D., 2000, The exact bootstrap mean and variance
of an L-estimator: Journal Royal Statistical Society B, v. 62, part 1, pp. 89--94.

Wang, D., and Hutson, A.D., 2013, Joint confidence region estimation of L-moments with an extension to right censored data: Journal of Applied Statistics, v. 40, no. 2, pp. 368--379.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lmoms}} }
\examples{
\dontrun{
   para <- vec2par(c(0,1), type="gum") # Parameters of Gumbel
   n <- 10; nmom <- 6; nsim <- 2000
   # X <- rlmomco(n, para) # This is commented out because
   # the sample below is from the Gumbel distribution as in para.
   # However, the seed for the random number generator was not recorded.
   X <- c( -1.4572506, -0.7864515, -0.5226538,  0.1756959,  0.2424514,
            0.5302202,  0.5741403,  0.7708819,  1.9804254,  2.1535666)
   EXACT.BOOTLMR <- lmoms.bootbarvar(X, nmom=nmom)
   LA <- EXACT.BOOTLMR$lambdavars
   LB <- LC <- rep(NA, length(LA))
   set.seed(n)
   for(i in 1:length(LB)) {
     LB[i] <- var(replicate(nsim,
                  lmoms(sample(X, n, replace=TRUE), nmom=nmom)$lambdas[i]))
   }
   set.seed(n)
   for(i in 1:length(LC)) {
     LC[i] <- var(replicate(nsim,
                  lmoms(rlmomco(n, para), nmom=nmom)$lambdas[i]))
   }
   print(LA) # The exact bootstrap variances of the L-moments.
   print(LB) # Bootstrap variances of the L-moments by actual resampling.
   print(LC) # Simulation of the variances from the parent distribution.

   # The variances for this example are as follows:
   #> print(LA)
   #[1] 0.115295563 0.018541395 0.007922893 0.010726508 0.016459913 0.029079202
   #> print(LB)
   #[1] 0.117719198 0.018945827 0.007414461 0.010218291 0.016290100 0.028338396
   #> print(LC)
   #[1] 0.17348653 0.04113861 0.02156847 0.01443939 0.01723750 0.02512031
   # The variances, when using simulation of parent distribution,
   # appear to be generally larger than those based only on resampling
   # of the available sample of only 10 values.

   # Interested users may inspect the exact bootstrap estimates of the
   # order statistics and the variance-covariance matrix.
   # print(EXACT.BOOTLMR$bootstrap.orderstatistics)
   # print(EXACT.BOOTLMR$varcovar.orderstatistics)

   # The output for these two print functions is not shown, but what follows
   # are the numerical confirmations from A.D. Hutson (personnal commun., 2012)
   # using his personnal algorithms (outside of R).
   # Date: Jul 2012, From: ahutson, To: Asquith
   # expected values the same
   # -1.174615143125091, -0.7537760316881618, -0.3595651823632459,
   # -0.028951905838698,  0.2360931764028858,  0.4614289985084462,
   #  0.713957210869635,  1.0724040932920058,  1.5368435379648948,
   #  1.957207045977329
   # and the first two values on the first row of the matrix are
   # 0.1755400544274771,  0.1306634198810892
}
\dontrun{
# Wang and Hutson (2013): Attempt to reproduce first entry of
# row 9 (n=35) in Table 1 of the reference, which is 0.878.
Xsq  <- qchisq(1-0.05, 2); n <- 35; nmom <- 4; nsim <- 1000
para <- vec2par(c(0,1), type="gum") # Parameters of Gumbel
eta  <- as.vector(lmorph(par2lmom(para))$ratios[3:4])
h <- 0
for(i in 1:nsim) {
   X <- rlmomco(n,para); message(i)
   EB <- lmoms.bootbarvar(X, nmom=nmom, verbose=FALSE)
   lmr    <- lmoms(X); etahat <- as.vector(lmr$ratios[c(3,4)])
   Pinv   <- EB$inverse.varcovar.tau34
   deta   <- (eta - etahat)
   LHS <- t(deta) %*% Pinv %*% deta # Definition 3.3
   if(LHS > Xsq) { # Comparison to Chi-squared distribution
      h <- h + 1 # increment because outside ellipse
      message("Outside: ",i, " ", h, " ", round(h/i, digits=3))
   }
}
message("Empirical Coverage Probability with Alpha=0.05 is ",
        round(1 - h/nsim, digits=3), " and count is", h)
# I have run this loop and recorded an h=123 for the above settings. I compute a
# coverage probability of 0.877, which agrees with Wang and Hutson (2013) within 0.001.
# Hence "very down the line" computations of lmoms.bootbarvar appear to be verified.
}
}
\keyword{L-moment (sample)}
\keyword{L-moment (variation)}
