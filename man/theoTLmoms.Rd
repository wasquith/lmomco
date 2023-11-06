\name{theoTLmoms}
\alias{theoTLmoms}
\title{The Theoretical Trimmed L-moments and TL-moment Ratios using Integration of the Quantile Function}
\description{
Compute the theoretrical trimmed L-moments (TL-moments) for a vector. The level of symmetrical or asymmetrical trimming is
specified. A theoretrical TL-moment in integral form is

\deqn{ \lambda^{(t_1,t_2)}_r = \underbrace{\frac{1}{r}}_{\stackrel{\mbox{average}}{\mbox{of terms}}}
                               \sum^{r-1}_{k=0} \overbrace{(-1)^k}^{\mbox{differences}}
          \underbrace{ r-1 \choose k }_{\mbox{combinations}}
        \frac{\overbrace{(r+t_1+t_2)!}^{\mbox{sample size}}\: I^{(t_1,t_2)}_r}
     {\underbrace{(r+t_1-k-1)!}_{\mbox{left tail}}
      \underbrace{(t_2+k)!}_{\mbox{right tail}}} \mbox{, in which }}

\deqn{ I^{(t_1,t_2)}_r = \int^1_0
                              \underbrace{x(F)}_{\stackrel{\mbox{quantile}}{\mbox{function}}} \times
                              \overbrace{F^{r+t_1-k-1}}^{\mbox{left tail}}
      \overbrace{(1-F)^{t_2+k}}^{\mbox{right tail}} \,\mathrm{d}F \mbox{,}}

where \eqn{x(F)} is the quantile function of the random variable \eqn{X} for nonexceedance probability \eqn{F}, \eqn{t_1} represents the trimming level of the \eqn{t_1}-smallest, \eqn{t_2} represents the trimming level of the \eqn{t_2}-largest values, \eqn{r} represents the order of the L-moments. This function loops across the above equation for each \code{nmom} set in the argument list. The function \eqn{x(F)} is computed through the \code{\link{par2qua}} function. The distribution type is determined using the \code{type} attribute of the \code{para} argument---the parameter object.

As of version 1.5.2 of \pkg{lmomco}, there exists enhanced error trapping on integration failures in \cr
\code{\link{theoTLmoms}}. The function now abandons operations should any of the integrations for the \eqn{r}th L-moment fail for reasons such as divergent integral or round off problems. The function returns NAs for all L-moments in \code{lambdas} and \code{ratios}.
}
\usage{
theoTLmoms(para, nmom=5, trim=NULL, leftrim=NULL, rightrim=NULL,
                 minF=0, maxF=1, quafunc=NULL,
                 nsim=50000, fold=5,
                 silent=TRUE, verbose=FALSE, ...)
}
\arguments{
  \item{para}{A distribution parameter object of this package such as by \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{trim}{Level of symmetrical trimming to use in the computations.
Although \code{NULL} in the argument list, the default is 0---the usual L-moment is returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample.}
  \item{rightrim}{Level of trimming of the right-tail of the sample.}
  \item{minF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-zero (but small) if you have a divergent integral.}
  \item{maxF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-unity (but close) if you have a divergent integral.}
  \item{quafunc}{An optional and arbitrary quantile function that simply needs to except a nonexceedance probability and the parameter object in \code{para}. This is a feature that permits computation of the L-moments of a quantile function that does not have to be implemented in the greater overhead hassles of the \pkg{lmomco} style. This feature might be useful for estimation of quantile function mixtures or those distributions not otherwise implemented in this package.}
   \item{nsim}{Simulation size for Monte Carlo integration is such is internally deemed necessary (see \code{silent} argument).}
   \item{fold}{The number of fractions or number of folds of \code{nsim}, which in other words, means that \code{nsim} is divided by \code{folds} and a loop creating \code{folds} integrations of \code{nsim/folds} is used from which the mean and mean absolute error of the integrand are computed. This is to try to recover similar output as \code{integrate()}.}
   \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}. If set true and the integral is probability divergent, Monte Carlo integration is triggered using \code{nsim} and \code{folds}. The user would have to set \code{verbose=TRUE} to then acquire the returned table in \code{integrations} of the integration passes including those are or are not Monte Carlo.}
   \item{verbose}{Toggle verbose output. Because the \R function \code{integrate} is used to perform the numerical integration, it might be useful to see selected messages regarding the numerical integration.}
   \item{...}{Additional arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the TL-moments. First element is
\eqn{\lambda^{(t_1,t_2)}_1}, second element is \eqn{\lambda^{(t_1,t_2)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau^{(t_1,t_2)}}, third element is \eqn{\tau^{(t_1,t_2)}_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal \code{NULL} if asymmetrical trimming was used.}
  \item{leftrim}{Level of left-tail trimming used in the computation.}
  \item{rightrim}{Level of right-tail trimming used in the computation.}
  \item{nsim}{Echo of the \code{nsim} argument if and only if at least one Monte Carlo integration was required, otherwise this is set to \dQuote{not needed} on the return.}
  \item{folds}{Echo of the \code{folds} argument if and only if at least one Monte Carlo integration was required, otherwise this is set to \dQuote{not needed} on the return.}
  \item{monte_carlo}{A logical vector of whether one or more Monte Carlo integrations was needed for the \code{r}-th index of the vector during the integrations for the \eqn{r}-th L-moment.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{theoTLmoms} or switched to \dQuote{theoLmoms} if this function was dispatched from \code{\link{theoLmoms}}.}
  \item{integrations}{If \code{verbose=TRUE}, then the results of the integrations are a data frame stored here. Otherwise, \code{integrations} is not present in the list.}
}
\note{
An extended example of a \emph{unique application} of the TL-moments is useful to demonstrate capabilities of the \pkg{lmomco} package API.  Consider the following example in which the analyst has 21 years of data for a given spatial location. Based on regional analysis, the highest value (the \code{outlier} = 21.12) is known to be exotically high but also documentable as not representing say a transcription error in the source database. The regional analysis also shows that the Generalized Extreme Value (GEV) distribution is appropriate.

The analyst is using a complex L-moment computational framework (say a software package called \bold{BigStudy.R}) in which only the input data are under the control of the analyst or it is too risky to modify \bold{BigStudy.R}. Yet, it is desired to somehow acquire robust estimation. The \code{outlier} value can be accommodated by estimating a pseudo-value and then simply make a substitution in the input data file for \bold{BigStudy.R}.

The following code initiates pseudo-value estimation by storing the original 20 years of data in variable \code{data.org} and then extending these data with the \code{outlier}. The usual sample L-moments are computed in \code{first.lmr} and will only be used for qualitative comparison. A 3-dimensional optimizer will be used for the GEV so the starting point is stored in \code{first.par}.
\preformatted{
  data.org  <- c(5.19, 2.58, 7.59, 3.22, 7.50, 4.05, 2.54, 9.00, 3.93, 5.15,
                 6.80, 2.10, 8.44, 6.11, 3.30, 5.75, 3.52, 3.48, 6.32, 4.07)
  outlier   <- 21.12;            the.data  <- c(data.org, outlier)
  first.lmr <- lmoms(the.data);  first.par <- pargev(first.lmr)
}

Robustness is acquired by computing the sample TL-moments such that the \code{outlier} is quantitatively removed by single trimming from the right side as the follow code shows:
\preformatted{
  trimmed.lmr <- TLmoms(the.data, rightrim=1, leftrim=0)
}

The objective now is to fit a GEV to the sample TL-moments in \code{trimmed.lmr}. However, the right-trimmed only (\eqn{t_1 = 0} and \eqn{t_2 = 1}) version of the TL-moments is being used and analytical solutions to the GEV for \eqn{t = (0,1)} are lacking or perhaps they are too much trouble to derive. The \code{theoTLmoms} function provides the avenue for progress because of its numerical integration basis for acquistion of the TL-moments. An objective function for the \eqn{t_2 = 1} TL-moments of the GEV is defined and based on the sum of square errors of the first three TL-moments:
\preformatted{
  "afunc" <- function(par, tarlmr=NULL, p=3) {
              the.par  <- vec2par(par, type="gev", paracheck=FALSE)
              fit.tlmr <- theoTLmoms(the.par, rightrim=1, leftrim=0)
              return(sum((tarlmr$lambdas[1:p] - fit.tlmr$lambdas[1:p])^2))
  }
}
and then optimize on this function and make a qualitative comparison between the original sample L-moments (untrimmed) to the equivalent L-moments (untrimmed) of the GEV having TL-moments equaling those in \code{trimmed.lmr}:
\preformatted{
  rt <- optim(first.par$para, afunc, tarlmr=trimmed.lmr)
  last.lmr <- lmomgev(vec2par(rt$par, type="gev"))

  message("# Original sample    L-moment lambdas: ",
           paste(round(first.lmr$lambdas[1:3], digits=4), collapse=" "))
  message("# Targeting back-fit L-moment lambdas: ",
           paste(round(last.lmr$lambdas[ 1:3], digits=4), collapse=" "))
  # Original sample    L-moment lambdas: 5.7981 1.8565 0.7287
  # Targeting back-fit L-moment lambdas: 5.5916 1.6501 0.5223
}
The primary result on comparison of the \eqn{\lambda_r} shows that the L-scale drops substantially as does L-skew: (\eqn{\tau_3 = 0.7287 / 1.8565 = 0.3925 \rightarrow \lambda_3^{(t_2{=}1)} = 0.5223 / 1.6501 = 0.3165}).

Now that the target L-moments (not TL-moments) are known (\code{last.lmr}), it is possible to optimize again on the value for the \code{outlier} that would provide the \code{last.lmr} within the greater computational framework in use by the analyst.
\preformatted{
  "bfunc" <- function(x, tarlmr=NULL, p=3) {
              sam.lmr <- lmoms(c(data.org, x))
              return(sum((tarlmr$lambdas[1:p] - sam.lmr$lambdas[1:p])^2))
  }
  suppressWarnings(outlier.rt <- optim(outlier, bfunc, tarlmr=last.lmr))
  # silence warning about 1D optimization with optim(), well behaved here

  pseudo.outlier <- round(outlier.rt$par, digits=2)
  final.lmr <- lmoms(c(data.org, pseudo.outlier))

  message("# Resulting new L-moment lambdas: ",
          paste(round(final.lmr$lambdas[1:3], digits=4), collapse=" "))
  # Resulting new L-moment lambdas: 5.5914 1.6499 0.5221

  message("# Pseudo-value for highest value: ", round(outlier.rt$par, digits=2))
  # Pseudo-value for highest value: 16.78
}
Where the second optimization shows that if the largest value for the 21 years of data is given a value of \eqn{16.78} instead of its original value of \eqn{21.12} that the sample L-moments (untrimmed) will be consistent as if the TL-moments \eqn{t = (0,1)} has been somehow used without resorting to a risky re-coding of the greater computational framework.
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{W.H. Asquith}
\seealso{\code{\link{theoLmoms}}, \code{\link{TLmoms}}, \code{\link{tlmr2par}}}
\examples{
para <- vec2par(c(0, 1), type='nor') # standard normal
TL00 <- theoTLmoms(para) # compute ordinary L-moments
TL30 <- theoTLmoms(para, leftrim=3, rightrim=0) # trim 3 smallest samples

# Let us look at the difference from simulation to theoretrical using
# L-kurtosis and asymmetrical trimming for generalized Lambda dist.
n     <- 100 # really a much larger sample should be used---for speed
P     <- vec2par(c(10000, 10000, 6, 0.4),type='gld')
Lkurt <- TLmoms(quagld(runif(n),P), rightrim=3, leftrim=0)$ratios[4]
theoLkurt <- theoTLmoms(P, rightrim=3, leftrim=0)$ratios[4]
Lkurt - theoLkurt # as the number for runif goes up, this
                  # difference goes to zero

# Example using the Generalized Pareto Distribution
# to verify computations from theoretical and sample stand point.
n      <- 100 # really a much larger sample should be used---for speed
P      <- vec2par(c(12, 34, 4),type='gpa')
theoTL <- theoTLmoms(P, rightrim=2, leftrim=4)
samTL  <- TLmoms(quagpa(runif(n),P), rightrim=2, leftrim=4)
del    <- samTL$ratios[3] - theoTL$ratios[3] # if n is large difference
                                             # is small
str(del)

\dontrun{
  "cusquaf" <- function(f, para, ...) { # Gumbel-Normal product
     g <- vec2par(c(para[1:2]), type="gum")
     n <- vec2par(c(para[3:4]), type="nor")
     return(par2qua(f,g)*par2qua(f,n))
  }
  para <- c(5.6, .45, 3, .3)
  theoTLmoms(para, quafunc=cusquaf) # L-skew = 0.13038711}

\dontrun{
  # This example has a divergent integral triggered on the last of the inner
  # loop of the 4th L-moment call. Monte Carlo (MC) integration is thus triggered.
  # The verbose=TRUE saves numerical or MC integration result table to the return.
  para   <-  vec2par(c(2.00,  2.00, -0.20, -0.55), type="kap")
  lmrbck <- lmomkap(   para, nmom=5)
  # print(lmrbck$lambdas) 3.1189568 1.9562688 0.4700229 0.4078741 0.1974055
  lmrthe <- theoTLmoms2(para, nmom=5, verbose=TRUE)              # seed dependent
  # print(lmrthe$lambdas) 3.1189569 1.9562686 0.4700227 0.4068539 0.1974049
  parkap(lmrbck)$para # 2.00       2.00     -0.20      -0.55
  parkap(lmrthe)$para # 2.018883  1.986761  -0.202422  -0.570451 # seed dependent}
}
\keyword{L-moment (theoretical)}
\keyword{L-moment (trimmed theoretical)}
