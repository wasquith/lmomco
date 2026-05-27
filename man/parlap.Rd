\name{parlap}
\alias{parlap}
\title{Estimate the Parameters of the Laplace Distribution}
\description{
This function computes the parameters of the Laplace distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. Relations between distribution parameters and sample L-moments are simple, but there are two methods to discuss.

First, basic theoretical formula of L-moment are simple for the Laplace and only use \eqn{\lambda_1} and \eqn{\lambda_2}, but these are used if \code{minvar=FALSE}:
\deqn{\xi    = \lambda_1\mbox{\ and\ } \alpha = 4\lambda_2/3\mbox{.}}
Hosking (1986, p. 56) reports these estimators as inefficient, and more efficient estimators can be obtained.

Second, Hosking (1986, p. 57) provides minimum-variance unbiased linear estimators with the best combination being shown below, and these are suggested for estimation the default implemented in \pkg{lmomco}. This minimum variance method jointly uses \eqn{\lambda_1, \lambda_2, \lambda_3}, and \eqn{\lambda_4}, and these are minimum variance unbiased estimators by Hosking (1986). The mathematical expressions are
\deqn{\xi = \lambda_1 - 50/31\times\lambda_3 \mbox{\ and}}
\deqn{\alpha = 1.4741\lambda_2 - 0.5960\lambda_4 = 1309/888\times\lambda_2 - 149/250\times\lambda_4\mbox{.}}
Astute readers will notice a fascinating aspect with these estimators. Consider that the distribution is symmetrical and therefore \eqn{\lambda_3 = 0 \equiv \lambda_3/\lambda_2 = \tau_3 = 0}, yet for a minimum-variance estimator, the \eqn{\lambda_3} resides within the estimation scheme.
}
\usage{
parlap(lmom, checklmom=TRUE, minvar=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally, this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{minvar}{A logical to use the simpler parameter estimation scheme shown above and this \code{minvar} will need to be set true if one where to check against output of \code{\link{theoLmoms}}.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.
  \item{type}{The type of distribution: \code{lap}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parlap}.}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: IBM Research Report RC12210, T.J. Watson Research Center, Yorktown Heights, New York.
}
\note{
The decision to use only one of the two systems of equations for Laplace fitting is largely arbitrary, but it seems most fitting to use four L-moments instead of two. The \bold{Examples} section provides evidence that the the \code{minvar=TRUE} is an appropriate default because we get smaller variance on the parameters at a potential cost in minor change in bias.

We can cross check the basic definition of the L-moments of the distribution (refer to \code{minvar=FALSE}) by the recipe in this section and noting (\code{\link{lmomlap}}) that \eqn{\lambda_1 = \xi},
\eqn{\lambda_2 = 3\alpha/4},
\eqn{\tau_3 = 0},
\eqn{\tau_4 = 17/72},
\eqn{\tau_5 = 0}, and
\eqn{\tau_6 = 31/360}. The results show reversibility between parameters and L-moments, which does not occur by the relations for \code{minvar=TRUE}.
\preformatted{
  lap  <- vec2par(c(1, 1), type="lap") # without loss of generalization
  lmomlap(lap) # L1=U, L2=3*A/4, T3=T5=0, T4=17/72=0.2361111, T6=31/360=0.0861111
  # $lambdas [1] 1.00000000 0.75000000 0.00000000 0.17708333 0.00000000 0.0645833
  # $ratios  [1] 0.00000000 0.75000000 0.00000000 0.23611111 0.00000000 0.0861111
  theoLmoms(lap, nmom=6) # by numerical integration of base L-moment definitions
  # $lambdas [1] 1.000e+00 0.7499999  -6.4*e-17  0.1770829  -1.0*e-15  0.06459555
  # $ratios  [1]        NA 0.7499999  -8.6*e-17  0.2361106  -1.4*e-15  0.08612742
}
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomlap}}, \code{\link{cdflap}}, \code{\link{pdflap}}, \code{\link{qualap}}
}
\examples{
lmr <- lmoms(rnorm(20))
parlap(lmr)

\dontrun{
  # Hosking (1986) shows that minimum variance estimators are readily
  # used via L-moments on the Laplace distribution.
  U <- 1; A <- 1; lap <- vec2par(c(1, 1), type="lap")
  nsim <- 50000; nsam <- c(5, 30, 50)
  for(n in nsam) {
    u1 <- u2 <- a1 <- a2 <- vector(mode="numeric", length=nsim)
    for(i in seq_len(nsim)) {
      lmr   <- lmoms(rlmomco(n, para=lap), no.stop=TRUE)
      par1  <- parlap(lmr, minvar=TRUE )$para # the minvar below
      par2  <- parlap(lmr, minvar=FALSE)$para # the direct below
      u1[i] <- par1[1]; a1[i] <- par1[2]
      u2[i] <- par2[1]; a2[i] <- par2[2]
    }
    message("n=", n)
    message("minvar: Bias[u1]=", mean(U - u1), " and Var[u1]=", var(u1))
    message("direct: Bias[u2]=", mean(U - u2), " and Var[u2]=", var(u2))
    message("minvar: Bias[a1]=", mean(A - a1), " and Var[a1]=", var(a1))
    message("direct: Bias[a2]=", mean(A - a2), " and Var[a2]=", var(a2))
  }
  # n=5
  # minvar: Bias[u1]=+0.00386375015371170 and Var[u1]=0.357661731795881
  # direct: Bias[u2]=+0.00560049550273044 and Var[u2]=0.400062259366783
  # minvar: Bias[a1]=-0.00037914342008174 and Var[a1]=0.243574424687960
  # direct: Bias[a2]=-0.00032851490958726 and Var[a2]=0.234478116729675
  # n=30
  # minvar: Bias[u1]=-0.00137763914162150 and Var[u1]=0.045735506455866
  # direct: Bias[u2]=-0.00069032577185824 and Var[u2]=0.066210402905459
  # minvar: Bias[a1]=0.000413407891428103 and Var[a1]=0.033904579056222
  # direct: Bias[a2]=0.000201995978024753 and Var[a2]=0.034924875477396
  # n=50
  # minvar: Bias[u1]=-0.00020392075808535 and Var[u1]=0.027220946727397
  # direct: Bias[u2]=-0.00041597093227831 and Var[u2]=0.040261174739697
  # minvar: Bias[a1]=-0.00053784244653237 and Var[a1]=0.020364246247174
  # direct: Bias[a2]=-0.00049308659560763 and Var[a2]=0.020926375838028 #}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Laplace}
