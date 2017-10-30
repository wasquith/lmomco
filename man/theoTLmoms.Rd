\name{theoTLmoms}
\alias{theoTLmoms}
\title{ The Theoretical Trimmed L-moments and TL-moment Ratios using Integration of the Quantile Function }
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
theoTLmoms(para, nmom=5, trim=NULL, leftrim=NULL,
           rightrim=NULL, verbose=FALSE, minF=0, maxF=1, quafunc=NULL)
}
\arguments{
  \item{para}{A distribution parameter object of this package such as by \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{trim}{Level of symmetrical trimming to use in the computations.
Although \code{NULL} in the argument list, the default is 0---the usual L-moment is returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample.}
  \item{rightrim}{Level of trimming of the right-tail of the sample.}
  \item{verbose}{Toggle verbose output. Because the \R function \code{integrate} is used to perform the numerical integration, it might be useful to see selected messages regarding the numerical integration.}
  \item{minF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-zero (but small) if you have a divergent integral.}
  \item{maxF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-unity (but close) if you have a divergent integral.}
  \item{quafunc}{An optional and arbitrary quantile function that simply needs to except a nonexceedance probability and the parameter object in \code{para}. This is a feature that permits computation of the L-moments of a quantile function that does not have to be implemented in the greater overhead hassles of the \pkg{lmomco} style. This feature might be useful for estimation of quantile function mixtures or those distributions not otherwise implemented in this package.}
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
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{theoTLmoms}.}
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
  "gev.afunc" <- function(par, tlmr=NULL) {
              the.par  <- vec2par(par, type="gev", paracheck=FALSE)
              fit.tlmr <- theoTLmoms(the.par, rightrim=1, leftrim=0)
              err1 <- (tlmr$lambdas[1] - fit.tlmr$lambdas[1])^2
              err2 <- (tlmr$lambdas[2] - fit.tlmr$lambdas[2])^2
              err3 <- (tlmr$lambdas[3] - fit.tlmr$lambdas[3])^2
              return(err1 + err2 + err3) # Sum of square errors
  }
}
and then optimize on this function and make a qualitative comparison between the original sample L-moments (untrimmed) to the equivalent L-moments (untrimmed) of the GEV having TL-moments equaling those in \code{trimmed.lmr}:
\preformatted{
  gev.rt   <- optim(first.par$para, gev.afunc, tlmr=trimmed.lmr)
  last.lmr <- lmomgev(vec2par(gev.rt$par, type="gev"))
  message("# Sample L-moment ratios: ",
           paste(round(first.lmr$ratios, digits=4), collapse=" "))
  message("# Target L-moment ratios: ",
           paste(round(last.lmr$ratios, digits=4), collapse=" "))
  # Sample L-moment ratios: NA 0.3202 0.3925 0.3113 0.2852
  # Target L-moment ratios: NA 0.2951 0.3165 0.2251 0.1304
}
The primary result on comparison of the \eqn{\tau_r} shows that the L-skew drops substantially (\eqn{\tau_3 = 0.393 \rightarrow \tau_3^{(t_2{=}1)} = 0.317}). The \eqn{\tau_4} and \eqn{\tau_5} are shown as well but since the GEV is not fit beyond the 3rd L-moment, these are not further considered.

Now that the \dQuote{target L-moments} are known (\code{last.lmr}), it is possible to optimize again on the value for the \code{outlier} that would provide the \code{last.lmr} within the greater computational framework in use by the analyst.
\preformatted{
  "lmr.afunc" <- function(x, target.lmr=NULL) {
              sam.lmr <- lmoms(c(data.org, x))
              return(sam.lmr$lambdas[1] - target.lmr$lambdas[1])
  }
  outlier.rt <- uniroot(lmr.afunc, interval=c(0, outlier), target.lmr=last.lmr)
  message("# Pseudo-value for highest value: ", round(outlier.rt$root, digits=2))
  # Pseudo-value for highest value: 16.78
}
Where the 2nd optimization shows that if the largest value for the 21 years of data is given a value of 16.78 in lieu of 21.12 that the sample L-moments (untrimmed) will be consistent as if the TL-moments \eqn{t = (0,1)} has been somehow used without resorting to a risky re-coding of the greater computational framework. Finally, the analyst can verify the pseudo-value by:
\preformatted{
   pseudo.outlier <- 16.78; print(lmoms(data.org, pseudo.outlier))
}
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{ W.H. Asquith}
\seealso{\code{\link{theoLmoms}}, \code{\link{TLmoms}}}
\examples{
para <- vec2par(c(0,1),type='nor') # standard normal
TL00 <- theoTLmoms(para) # compute ordinary L-moments
TL30 <- theoTLmoms(para,leftrim=3,rightrim=0) # trim 3 smallest samples

# Lets look at the difference from simulation to theoretrical using
# L-kurtosis and asymmetrical trimming for generalized Lambda dist.
P <- vec2par(c(10000,10000,6,.4),type='gld')
Lkurt <- TLmoms(quagld(runif(100),P),rightrim=3,leftrim=0)$ratios[4]
theoLkurt <- theoTLmoms(P,rightrim=3,leftrim=0)$ratios[4]
Lkurt - theoLkurt # as the number for runif goes up, this
                  # difference goes to zero

# Example using the Generalized Pareto Distribution
# to verify computations from theoretical and sample stand point.
n      <- 100 # really a much larger sample should be used---for speed
P      <- vec2par(c(12,34,4),type='gpa')
theoTL <- theoTLmoms(P,rightrim=2,leftrim=4)
samTL  <- TLmoms(quagpa(runif(n),P),rightrim=2,leftrim=4)
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
}
\keyword{L-moment (theoretical)}
\keyword{L-moment (trimmed theoretical)}
