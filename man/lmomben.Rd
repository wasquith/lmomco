\encoding{utf8}
\name{lmomben}
\alias{lmomben}
\title{L-moments of the Benford Distribution}
\description{
\emph{Experimental}---This function returns previously numerical estimations of the L-moments of the Benford distribution (Benford's Law) given parameters defining the number of first M-significant digits and the numeric base.

For the first significant digits (\eqn{d \in 1, \cdots, 9}) (base 10) (designate as \eqn{m = 1}), the L-moments were estimated through very large sample-size simulation and sample L-moments computed ( \code{\link{lmoms}}), direct numerical integration (\code{\link{theoLmoms}}), and through numerical integration of the probability weighted moments and conversion to L-moments (\code{\link{pwm2lmom}}) as
\deqn{\lambda_1 = 3.43908699617500524\mbox{,}}
\deqn{\lambda_2 = 1.34518434179517077\mbox{,}}
\deqn{\tau_3 = 0.24794090889493661\mbox{, and}}
\deqn{\tau_4 = 0.01614509742647182\mbox{.}}

For the first two-significant digits (\eqn{d \in 10, \cdots, 99}) (base 10) (designate as \eqn{m = 2}), the L-moments were estimated through very large sample-size simulation, direct numerical integration (\code{\link{theoLmoms}}), and through numerical integration of the probability weighted moments and conversion to L-moments (\code{\link{pwm2lmom}}) as
\deqn{\lambda_1 = 38.59062918136093145\mbox{,}}
\deqn{\lambda_2 = 13.81767809210059283\mbox{,}}
\deqn{\tau_3 = 0.22237541787527126\mbox{, and}}
\deqn{\tau_4 = 0.03541037418894027\mbox{.}}

For the first three-significant digits (\eqn{d \in 100, \cdots, 999}) (base 10) (designate as \eqn{m = 3}), the L-moments were estimated through very large sample-size simulation, direct numerical integration (\code{\link{theoLmoms}}), and through numerical integration of the probability weighted moments and conversion to L-moments (\code{\link{pwm2lmom}}) as
\deqn{\lambda_1 = 390.36783537821605705\mbox{,}}
\deqn{\lambda_2 = 138.21917489739223583\mbox{,}}
\deqn{\tau_3 = 0.22192482374529940\mbox{, and}}
\deqn{\tau_4 = 0.03571514686148788\mbox{.}}

\bold{Source of the L-moments}---The script \code{inst/doc/benford/compLmomsBenford.R} in the \pkg{lmomco} package sources is the authoritative source of the computation of the L-moments shown. Three methods are used, and the arithmetic average of the three provides the L-moments: (1) Probability-weighted simulation of the probability mass function (PMF) is used in very large sample size and sample L-moments computed by \code{\link{lmoms}}, (2) direct numerical integration for the theoretical L-moments of the quantile function (\code{\link{quaben}}) of the distribution that itself is from the cumulative distribution function (\code{\link{cdfben}}) that itself is from the PMF (\code{\link{pmfben}}), and (3) direct numerical integration of the probability-weighted moments of the quantile function (\code{\link{quaben}}) and subsequent linear system of equations to compute the L-moments. Each of the aforementioned methods result in numerical differences say at about the fourth decimal. (No previous description of the L-moments of the Benford distribution appear extant in the literature in July 2024.)
}
\usage{
lmomben(para=list(para=c(1, 10)), ...)
}
\arguments{
    \item{para}{The number of first M-significant digits followed by the numerical base (only base10 supported) and the list structure mimics similar uses of the \pkg{lmomco} list structure. Default are the first significant digits and hence the digits 1 through 9.}
  \item{...}{Additional arguments to pass (not likely to be needed but changes in base handling might need this).}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomben}.}
}
\note{\bold{Hypothesis Testing}---Let the squared Euclidean distance of the L-moments (not the L-moment ratios) between the first four sample L-moments (\eqn{\hat{\lambda}_r}) and the theoretical versions (\eqn{\lambda}) provided by this function be defined as

\deqn{D^2 = (\hat{\lambda}_1 - \lambda_1)^2 + (\hat{\lambda}_2 - \lambda_2)^2 + (\hat{\lambda}_3 - \lambda_3)^2 + (\hat{\lambda}_4 - \lambda_4)^2\mbox{.}}

Let \eqn{\alpha \in 0.10, 0.05, 0.01, 0.005, 0.001} be upper tail probability levels (statistical significance thresholds). Let \eqn{m} denote the number of significant digits (\eqn{m \in 1, 2, 3}) in base 10 and \eqn{n} denote sample size. Let \eqn{\gamma = -\mathrm{log}(-\mathrm{log}(\alpha))} be a transformation (in the style of a Gumbel reduced variate) (\code{\link{prob2grv}}). Using extensive simulation for many sample sizes, the \eqn{\alpha} values, and computing \eqn{D^2(\alpha, m; n)}, it can be shown that the critical values for the \eqn{D^2} distances are

\deqn{D^2(\alpha) = \frac{1}{n}\,\mathrm{exp}\bigl[(-2.6607150 + 4.6154937m) - 1.217283\gamma\bigr]\mbox{,}}

wherein linear regression was used to estimate relation between each \eqn{D^2} and \eqn{n \ge 5} and the coefficients subsequently subjected to linear regression as functions of \eqn{\alpha}. The \bold{Examples} shows an implementation of the critical values.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfben}}, \code{\link{pmfben}}, \code{\link{quaben}} }
\examples{
lmomben(para=list(para=c(3, 10)))

\dontrun{
  # Code suitable for study of performance of Cho and Gaines D
  # against using the first for L-moments with controls for having
  # the Benford distribution as the true parent or alternative
  # distributions fit to the the L-moments of the Benford for the
  # first significant digit.
  # https://en.wikipedia.org/wiki/Benford\%27s_law#Statistical_tests
  ChoGainesD <- function(x) {
    n <- length(x)
    d <- sapply(1:9, function(d) (length(x[x == d])/n - log10(1+1/d))^2)
    return(sqrt(n * sum(d)))
  }
  CritChoGainesD <- function(alpha=c("0.1", "0.05", "0.01")) {
    alpha <- as.character(as.numeric( alpha ))
    alpha <-   as.numeric(match.arg(  alpha ))
    if(alpha == 0.10) return(1.212)
    if(alpha == 0.05) return(1.330)
    if(alpha == 0.01) return(1.569)
    return(NULL)
  }
  D2lmom <- function(x, theolmr=NULL) {
    lmr <- lmoms(x)
    sum((lmr$lambdas[1:4] - theolmr$lambdas[1:4])^2)
  }
  CritD2lmom <-
    function(m, n, alpha=c("0.1", "0.05", "0.01", "0.005", "0.001")) {
      alpha <- as.character(as.numeric( alpha ))
      alpha <-   as.numeric(match.arg(  alpha ))
      exp((-2.6607150 + 4.6154937*m) - 1.217283*(-log(-log(alpha))))/n
  }

  nsim <- 2E4; n <- 100; alpha <- 0.05
  is_Benford_parent <- FALSE

  CritCGD <- CritChoGainesD(  alpha=alpha )
  CritLMR <- CritD2lmom(1, n, alpha=alpha )
  bens <- 1:9; pmf <- log10(1 + 1/bens) # for the Benford being true
  benlmr <- lmomben(list(para=c(1, 10))); dtype <- "nor" # Normal (say)
  parent <- lmom2par(benlmr, type=dtype)

  DF <- NULL
  ix <- seq(1, n, by=2)
  for(i in 1:nsim) {
    if(is_Benford_parent) {
      x <- sample(bens, n, replace=TRUE, prob=pmf)
    } else {
      x <- rlmomco(n, parent) # simulate from the parent
      x <- unlist(strsplit(sprintf("\%2.0E", x), "E"))[ix]
      x <- as.integer(x) # complete extraction of the first digit
    }
    CGD    <- ChoGainesD(x)
    LMR    <- D2lmom(x, theolmr=benlmr)
    rejCGD <- ifelse(CGD > CritCGD, TRUE, FALSE)
    rejLMR <- ifelse(LMR > CritLMR, TRUE, FALSE)
    DF <- rbind(DF, data.frame(CGD=rejCGD, LMR=rejLMR))
  }
  print(summary(DF))
  if(is_Benford_parent) { # H0 is True
    CGDpct <- 100*(sum(as.numeric(DF$CGD)) / nsim - alpha) / alpha;
    LMRpct <- 100*(sum(as.numeric(DF$LMR)) / nsim - alpha) / alpha;
    message("The ChoGainesD rejection rate for alpha=", alpha,
            " is ", sum(as.numeric(DF$CGD)) / nsim,
            " (", round(CGDpct, digits=2), " percent difference).")
    message("The   D2lmom   rejection rate for alpha=", alpha,
            " is ", sum(as.numeric(DF$LMR)) / nsim,
            " (", round(LMRpct, digits=2), " percent difference).")
  } else { # H0 is False
    acceptH0_H0false_CDG <- sum(as.numeric(! DF$CGD)) / nsim
    acceptH0_H0false_LMR <- sum(as.numeric(! DF$LMR)) / nsim
    betaCDG <- round(1 - acceptH0_H0false_CDG, digits=2)
    betaLMR <- round(1 - acceptH0_H0false_LMR, digits=2)
    message("Power of ChoGainesD = ", betaCDG, ".")
    message("Power of   D2lmom   = ", betaLMR, ".")
  } #}
}
\keyword{distribution (discrete)}
\keyword{cumulative distribution function}
\keyword{Distribution: Benford}
\concept{Benford Law}
