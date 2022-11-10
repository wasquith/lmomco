\name{mps2par}
\alias{mps2par}
\title{Use Maximum Product of Spacings to Estimate the Parameters of a Distribution}
\description{
This function uses the method of maximum product of spacings (MPS; maximum spacing estimation or maximum product of spacings estimation) to estimate the parameters of a distribution. The method is based on maximization of the \emph{geometric mean} of probability spacings in the data where the spacings are defined as the differences between the values of the cumulative distribution function, \eqn{F(x)}, at sequential data indices.

MPS (Dey \emph{et al.}, 2016, pp. 13--14) is an optimization problem formed by maximizing the geometric mean of the spacing between consecutively ordered observations standardized to a U-statistic. Let \eqn{\Theta} represent a vector of parameters for a candidate fit of \eqn{F(x|\Theta)}, and let \eqn{U_i(\Theta) = F(X_{i:n}|\Theta)} be the nonexceedance probabilities of the observed values of the order statistics \eqn{x_{i:n}} for a sample of size \eqn{n}.  Define the differences
\deqn{D_i(\Theta) = U_i(\Theta) - U_{i-1}(\Theta)\mbox{\ for\ } i = 1, \ldots, n+1\mbox{,}}
with the additions to the vector \eqn{U} of \eqn{U_0(\Theta) = 0} and \eqn{U_{n+1}(\Theta) = 1}. The objective function is
\deqn{M_n(\Theta) = - \sum_{i=1}^{n+1} \log\, D_i(\Theta)\mbox{,}}
where the \eqn{\Theta} for a maximized \eqn{{-}M_n} represents the parameters fit by MPS. Some authors to keep with the idea of geometric mean include factor of \eqn{1/(n+1)} for the definition of \eqn{M_n}. Whereas other authors (Shao and Hahn, 1999, eq. 2.0), show
\deqn{S_n(\Theta) = (n+1)^{-1} \sum_{i=1}^{n+1} \log[(n+1)D_i(\Theta)]\mbox{.}}
So it seems that some care is needed when considering the implementation when the value of \dQuote{the summation of the logarithms} is to be directly interpreted. Wong and Li (2006) provide a salient review of MPS in regards to an investigation of maximum likelihood (MLE), MPS, and probability-weighted moments (\code{\link{pwm}}) for the GEV (\code{\link{quagev}}) and GPA (\code{\link{quagpa}}) distributions. Finally, Soukissian and Tsalis (2015) also study MPS, MLE, L-moments, and several other methods for GEV fitting.

If the initial parameters have a support inside the range of the data, infinity is returned immediately by the optimizer and further action stops and the parameters returned are \code{NULL}. For the implementation here, if \code{check.support} is true, and the initial parameter estimate (if not provided and acceptable by \code{para.int}) by default will be seeded through the method of L-moments (unbiased, \code{\link{lmoms}}), which should be close and convergence will be fairly fast if a solution is possible. If these parameters can not be used for spinup, the implementation will then attempt various probability-weighted moment by plotting position (\code{\link{pwm.pp}}) converted to L-moments (\code{\link{pwm2lmom}}) as part of an extended attempt to find a support of the starting distribution encompass the data. Finally, if that approach fails, a last ditch effort using starting parameters from maximum likelihood computed by a default call to \code{\link{mle2par}} is made. Sometimes data are pathological and user supervision is needed but not always successful---MPS can show failure for certain samples and(or) choice of distribution.

It is important to remark that the support of a fitted distribution is not checked within the loop for optimization once spun up. The reasons are twofold: (1) The speed hit by repeated calls to \code{\link{supdist}}, but in reality (2) PDFs in \pkg{lmomco} are supposed to report zero density for outside the support of a distribution (see NEWS) and for the \eqn{-\log(D_i(\Theta)\rightarrow 0) \rightarrow \infty} and hence infinity is returned for that state of the optimization loop and alternative solution will be tried.

As a note, if all \eqn{U} are equally spaced, then \eqn{|M(\Theta)| = I_o = (n+1)\log(n+1)}. This begins the concept towards goodness-of-fit. The \eqn{M_n(\Theta)} is a form of the Moran-Darling statistic for goodness-of-fit. The \eqn{M_n(\Theta)} is a Normal distribution with
\deqn{\mu_M \approx (n+1)[\log(n+1)+\gamma{}] - \frac{1}{2} - \frac{1}{12(n+1)}\mbox{,}}
\deqn{\sigma_M \approx (n+1)\biggl(\frac{\pi^2}{6\,{}} - 1\biggr)-\frac{1}{2} - \frac{1}{6(n+1)}\mbox{,}}
where \eqn{\gamma \approx 0.577221} (Euler--Mascheroni constant, \code{-digamma(1)}) or as the definite integral

\deqn{\gamma^\mathrm{Euler}_{\mathrm{Mascheroni}} = -\int_0^\infty \mathrm{exp}(-t) \log(t)\; \mathrm{d}{t}\mbox{,}}

An extension into small samples using the Chi-Square distribution is
\deqn{A = C_1 + C_2\times\chi^2_n\mbox{,}}
where
\deqn{C_1 = \mu_M - \sqrt{\frac{\sigma^2_M\,n}{2}}\mbox{\ and\ }C_2 = \sqrt{\frac{\sigma^2_M}{2n}}\mbox{,}}
and where \eqn{\chi^2_n} is the Chi-Square distribution with \eqn{n} degrees of freedom. A test statistic is
\deqn{T(\Theta) = \frac{M_n(\Theta) - C_1 + \frac{p}{2}}{C_2}\mbox{,}}
where the term \eqn{p/2} is a bias correction based on the number of fitted distribution parameters \eqn{p}. The null hypothesis that the fitted distribution is correct is to be rejected if \eqn{T(\Theta)} exceeds a critical value from the Chi-Square distribution. The MPS method has a relation to maximum likelihood (\code{\link{mle2par}}) and the two are asymptotically equivalent.

\bold{Important Remark Concerning Ties}---Ties in the data cause \emph{instant degeneration} with MPS and must be mitigated for and thus attention to this documentation and even the source code itself is required.
}
\usage{
mps2par(x, type, para.int=NULL, ties=c("bernstein", "rounding", "density"),
            delta=0, log10offset=3, get.untied=FALSE, check.support=TRUE,
            moran=TRUE, silent=TRUE, null.on.not.converge=TRUE,
            ptransf=  function(t) return(t),
            pretransf=function(t) return(t),
            mle2par=TRUE, ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{type}{Three character (minimum) distribution type (for example, \code{type="gev"}, see \code{\link{dist.list}}).}
  \item{para.int}{Initial parameters as a vector \eqn{\Theta} or as an \pkg{lmomco} parameter \dQuote{object} from say \code{\link{vec2par}}. If a vector is given, then internally \code{\link{vec2par}} is called with distribution equal to \code{type}.}
  \item{ties}{Ties cause degeneration in the computation of \eqn{M(\Theta)}:\cr
  Option \code{bernstein} triggers a smoothing of only the ties using the \code{\link{dat2bernqua}} function---Bernstein-type smoothing for ties is likely near harmless when ties are near the center of the distribution, but of course caution is advised if ties exist near the extremal values; the settings for \code{log10offset} and \code{delta} are ignored if \code{bernstein} is selected Also for a tie-run having an odd number of elements, the middle tied value is left as original data.\cr
  Option \code{rounding} triggers two types of adjustment: if \code{delta > 0} then a round-off error approach inspired by Cheng and Stephens (1989, eq. 4.1) is used (see \bold{Note}) and \code{log10offset} is ignored, but if \code{delta=0}, then \code{log10offset} is picked up as an order of magnitude offset (see \bold{Note}). Use of options \code{log10offset} and \code{delta} are likely to not keep a middle unmodified in an odd-length, tie-run in contrast to use of \code{bernstein}.\cr
  Option \code{density} triggers the substitution of the probability density \eqn{g(x_{i:n}|\Theta)} at the \eqn{i}th tie from the current fit of the distribution. \bold{Warning}---It appears that inference is lost almost immediately because the magnitude of \eqn{M_n} losses meaning because probability densities are not in the same scale as changes in probabilities exemplified by the \eqn{D_i}. This author has not yet found literature discussing this, but density substitution is a recognized strategy.}
  \item{delta}{The optional \eqn{\delta} value if \eqn{\delta > 0} and if \code{ties=}\code{"rounding"}.}
  \item{log10offset}{The optional base-10 logarithmic offset approach to roundoff errors if \code{delta=0} and if \code{ties=}\code{"rounding"}.}
  \item{get.untied}{A logical to populate a \code{ties} element in the returned \code{list} with the untied-pseudo data as it was made available to the optimizer and the number of iternations required to exhaust all ties. An emergency break it implemented if the number of iterations appears to be blowing up.}
  \item{check.support}{A logical to trigger a call to \code{\link{supdist}} to compute the support of the distribution at the initial parameters. As mentioned, MPS degenerates if \code{min(x)} \eqn{<} the lower support or if \code{max(x)} \eqn{>} the upper support. Regardless of the setting of \code{check.support} and \code{NULL} will be returned because this is what the optimizer will do anyway.}
  \item{moran}{A logical to trigger the goodness-of-fit test described previously.}
  \item{silent}{A logical to silence the \code{try()} function wrapping the \code{optim()} function and to provide a returned list of the optimization output.}
  \item{null.on.not.converge}{A logical to trigging simple return of \code{NULL} if the \code{optim()} function returns a nonzero convergence status.}
  \item{ptransf}{An optional parameter transformation function (see \bold{Examples}) that is useful to guide the optimization run. For example, suppose the first parameter of a three parameter distribution resides in the positive domain, then \cr \code{ptransf(t) = }
\code{function(t) c(log(t[1]), t[2], t[3])}.}
  \item{pretransf}{An optional parameter retransformation function (see \bold{Examples}) that is useful to guide the optimization run. For example, suppose the first parameter of a three parameter distribution resides in the positive domain, then \cr \code{pretransf(t) = } \code{function(t) c(exp(t[1]), t[2], t[3])}.}
  \item{mle2par}{A logical to turn off the potential last attempt at maximum likelihood estimates of a valid seed as part of \code{check.support=TRUE}.}
  \item{...}{Additional arguments for the \code{optim()} function and other uses.}
}
\value{
  An \R \code{list} is returned.  This list should contain at least the following items, but some distributions such as the \code{revgum} have extra.
  \item{type}{The type of distribution in three character (minimum) format.}
  \item{para}{The parameters of the distribution.}
  \item{source}{Attribute specifying source of the parameters.}
  \item{para.int}{The initial parameters. Warning to users, when inspecting returned values make sure that one is referencing the MPS parameters in \code{para} and not those shown in \code{para.int}!}
  \item{optim}{An optional \code{list} of returned content from the optimizer if not \code{silent}.}
  \item{ties}{An optional \code{list} of untied-pseudo data and number of iterations required to achieve no ties (usually unity!) if and only if there were ties in the original data, \code{get.untied} is true, and \code{ties != "density"}.}
  \item{MoranTest}{An optional \code{list} of returned values that will include both diagnostics and statistics. The diagnostics are the computed \eqn{\mu_M(n)}, \eqn{\sigma^2_M(n)}, \eqn{C_1}, \eqn{C_2}, and \eqn{n}. The statistics are the minimum value \eqn{I_o} theoretically attainable \eqn{|M_n(\Theta)|} for equally spaced differences, the minimized value \eqn{M_n(\Theta)}, the \eqn{T(\Theta)}, and the corresponding \code{p.value} from the upper tail of the \eqn{\chi^2_n} distribution.}
}
\note{
During optimization, the objective function requires evaluation at the initial parameters and must be finite. If \code{Inf} is returned on first call to the objective function, then a warning like this
\preformatted{
  optim() attempt is NULL
}
should be seen. The \code{silent} by default though will silence this error. Error trapping for the estimated support of the distribution from the initial parameter values is made by \code{check.support=TRUE} and verbose warnings given to help remind the user. Considerable attempt is made internally to circumvent the appearance of the above error.

More specifically, an MPS solution degenerates when the fitted distribution has a narrower support than the underlying data and artificially \dQuote{ties} show up within the objective function even if the original data lacked ties or were already mitigated for. The user's only real recourse is to try fitting another distribution either by starting parameters or even distribution type. Situations could arise for which carefully chosen starting parameters could permit the optimizer to keep its simplex within the viable domain. The MPS method is sensitive to tails of a distribution having asymoptic limits as \eqn{F \rightarrow 0^{+}} or \eqn{F \rightarrow 1^{-}}.

The Moran test can be quickly checked with highly skewed and somewhat problematic data by
\preformatted{
  gev <- vec2par(c(4,0.3,-0.2), type="gev"); nsim <- 5000
  G <- replicate(nsim, mps2par(rlmomco(100, gev), # extract the p-values
                               type="gev")$MoranTest$statistics[4])
  G <- unlist(G) # unlisting required if NULLs came back from mps2par()
  length(G[G <= 0.05])/length(G) # 0.0408 (!=0.05 but some fits not possible)
  V <- replicate(nsim, mps2par(rlmomco(100, gev),
                               type="nor")$MoranTest$statistics[4])
  V <- unlist(V) # A test run give 4,518 solutions
  length(V[V <= 0.05])/length(V) # 0.820 higher because not gev used
  W <- replicate(nsim, mps2par(rlmomco(100, gev),
                               type="glo")$MoranTest$statistics[4])
  W <- unlist(W)
  length(W[W <= 0.05])/length(W) # 0.0456 higher because not gev used but
  # very close because of the proximity of the glo to the gev for the given
  # L-skew of the parent: lmomgev(gev)$ratios[3] = 0.3051
}

Concerning round-off errors, the Cheng and Stephens (1989, eq. 4.1) approach is to assume that the round-off errors are \eqn{x \pm \delta}, compute the upper and lower probabilities \eqn{f} for \eqn{f_L \mapsto x - \delta} and \eqn{f_U \mapsto x + \delta}, and then prorate the \eqn{D_i} in even spacings of \eqn{1/(r-1)} where \eqn{r} is the number of tied values in a given tie-run. The approach for \code{mps2par} is similar but simplies the algorithm to evenly prorate the \eqn{x} values in a tie-run. In other words, the current implementation is to actually massage the data before passage into the optimizer.  If the \eqn{\delta = 0}, a base-10 logarithmic approach will be used in which, the order of magnitude of the value in a tie-run is computed and the \code{log10offset} subtracted to approximate the roundoff but recognize that for skewed data the roundoff might be scale dependent. The default treats a tie of three \eqn{x_i = 15{,}000} as \eqn{x_{i|r}=14{,}965.50; 15{,}000.00; 15{,}034.58}. In either approach, an iterative loop is present to continue looping until no further ties are found---this is made to protect against the potential for the algorithm to create new ties. A sorted vector of the final data for the optimize is available in the \code{ties} element of the returned list if and only if ties were originally present, \code{get.untied=TRUE}, and \code{ties != "density"}. Ties and compensation likely these prorations can only make \eqn{M(\Theta)} smaller, and hence the test becomes conservative.

A note of other MPS implementations in \R is needed. The \pkg{fBasics} and \pkg{gld} packages both provide for MPS estimation for the generalized lambda distribution. The salient source files and code chunks are shown. First, consider package \pkg{fBasics}:
\preformatted{
  fBasics --> dist-gldFit.R --> .gldFit.mps -->
            f = try(-typeFun(log(DH[DH > 0])), silent = TRUE)
}
where it is seen that \eqn{D_i = 0} are ignored! Such a practice does not appear efficacious during development and testing of the implementation in \pkg{lmomco}, parameter solutions very substantially different than reason can occur or even failure of convergence by the \pkg{fBasics} implementation. Further investigation is warranted. Second, consider package \pkg{gld}:
\preformatted{
  gld --> fit_fkml.R --> fit_fkml.c --> method.id == 2:
  # If F[i]-F[i-1] = 0, replace by f[i-1]
  #                      (ie the density at smaller observation)
}
which obviously make the density substitution for ties as well \code{ties=}\code{"density"} for the implementation here. Testing indicates that viable parameter solutions will result with direct insertion of the density in the case of ties. Interference, however, of the \eqn{M_n} is almost assuredly to be greatly weakened or destroyed depending on the shape of the probability density function or a large number of ties. The problem is that the sum of the \eqn{D_i} are no longered ensured to sum to unity. The literature appears silent on this particular aspect of MPS, and further investigation is warranted.

The \pkg{eva} package provides MPS for GEV and GPD. The approach there does not appear to replace changes of zero by density but to insert a \dQuote{smallness} in conjunction with other conditioning checking (only the \code{cond3} is shown below) and a curious penalty of \code{1e6}. The point is that different approaches have been made by others.
\preformatted{
  eva --> gevrFit --> method="mps"
  cdf[(is.nan(cdf) | is.infinite(cdf))] <- 0
  cdf <- c(0, cdf, 1); D <- diff(cdf); cond3 <- any(D < 0)
  ## Check if any differences are zero due to rounding and adjust
  D <- ifelse(D <= 0, .Machine$double.eps, D)
  if(cond1 | cond2 | cond3) { abs(sum(log(D))) + 1e6 } else { -sum(log(D)) }
}
Let us conclude with an example for the GEV between \pkg{eva} and \pkg{lmomco} and note sign difference in definition of the GEV shape but otherwise a general similarity in results:
\preformatted{
  X <- rlmomco(97, vec2par(c(100,12,-.5), type="gev"))
  pargev(lmoms(X))$para
       #                  xi                alpha                kappa
       #         100.4015424           12.6401335           -0.5926457
  eva::gevrFit(X, method="mps")$par.ests
       #Location (Intercept)    Scale (Intercept)    Shape (Intercept)
       #         100.5407709           13.5385491            0.6106928
}
}
\author{W.H. Asquith}
\references{
Cheng, R.C.H., Stephens, M.A., 1989, A goodness-of-fit test using Moran's statistic with estimated parameters: Biometrika, v. 76, no. 2, pp. 385--392.

Dey, D.K., Roy, Dooti, Yan, Jun, 2016, Univariate extreme value analysis, chapter 1, \emph{in} Dey, D.K., and Yan, Jun, eds., Extreme value modeling and risk analysis---Methods and applications: Boca Raton, FL, CRC Press, pp. 1--22.

Shao, Y., and Hahn, M.G., 1999, Strong consistency of the maximum product of spacings estimates with applications in nonparametrics and in estimation of unimodal densities: Annals of the Institute of Statistical Mathematics, v. 51, no. 1, pp. 31--49.

Soukissian, T.H., and Tsalis, C., 2015, The effect of the generalized extreme value distribution parameter estimation methods in extreme wind speed prediction: Natural Hazards, v. 78, pp. 1777--1809.

Wong, T.S.T., and Li, W.K., 2006, A note on the estimation of extreme value distributions using maximum product of spacings: IMS Lecture Notes, v. 52, pp. 272--283.
}
\seealso{\code{\link{lmom2par}}, \code{\link{mle2par}}, \code{\link{tlmr2par}}}
\examples{
\dontrun{
pe3 <- vec2par(c(4.2, 0.2, 0.6), type="pe3") # Simulated values should have at least
X <- rlmomco(202, pe3); Xr  <- round(sort(X), digits=3) # one tie-run after rounding,
mps2par(X,  type="pe3")$para      # and the user can observe the (minor in this case)
mps2par(Xr, type="pe3")$para      # effect on parameters.
# Another note on MPS is needed. It is not reflection symmetric.
mps2par( X, type="pe3")$para
mps2par(-X, type="pe3")$para}

\dontrun{
# Use 1,000 replications for sample size of 75 and estimate the bias and variance of
# the method of L-moments and maximum product spacing (MPS) for the 100-year event
# using the Pearson Type III distribution.
set.seed(1596)
nsim <- 1000; n <- 75; Tyear <- 100; type <- "pe3"
parent.lmr <- vec2lmom(c(5.5, 0.15, 0.03))   # L-moments of the "parent"
parent  <- lmom2par(parent.lmr, type="pe3")  # "the parent"
Q100tru <- qlmomco(T2prob(Tyear), parent)    # "true value"
Q100lmr <- Q100mps <- rep(NA, nsim)          # empty vectors
T3lmr <- T4lmr <- T3mps <- T4mps <- rep(NA, nsim)
for(i in 1:nsim) { # simulate from the parent, compute L-moments
   tmpX <- rlmomco(n, parent); lmrX <- lmoms(tmpX)
   if(! are.lmom.valid(lmrX)) { # quiet check on viability
     lmrX <- pwm2lmom(pwms.pp(tmpX)) # try a pwm by plotting positions instead
     if(! are.lmom.valid(lmrX)) next
   }
   lmrpar <- lmom2par(lmrX, type=type)                  # Method of L-moments
   mpspar <-  mps2par(tmpX, type=type, para.int=lmrpar) # Method of MPS
   if(! is.null(lmrpar)) {
      Q100lmr[i] <- qlmomco(T2prob(Tyear), lmrpar); lmrlmr <- par2lmom(lmrpar)
      T3lmr[i] <- lmrlmr$ratios[3]; T4lmr[i] <- lmrlmr$ratios[4]
   }
   if(! is.null(mpspar)) {
      Q100mps[i] <- qlmomco(T2prob(Tyear), mpspar); mpslmr <- par2lmom(mpspar)
      T3mps[i] <- mpslmr$ratios[3]; T4mps[i] <- mpslmr$ratios[4]
   }
}
print(summary(Q100tru - Q100lmr)) # Method of L-moment   (mean = -0.00176)
print(summary(Q100tru - Q100mps)) # Method of MPS        (mean = -0.02746)
print(var(Q100tru - Q100lmr, na.rm=TRUE)) # Method of L-moments (0.009053)
print(var(Q100tru - Q100mps, na.rm=TRUE)) # Method of MPS       (0.009880)
# CONCLUSION: MPS is very competitive to the mighty L-moments.

LMR <- data.frame(METHOD=rep("Method L-moments",        nsim), T3=T3lmr, T4=T4lmr)
MPS <- data.frame(METHOD=rep("Maximum Product Spacing", nsim), T3=T3mps, T4=T4mps)
ZZ <- merge(LMR, MPS, all=TRUE)
boxplot(ZZ$T3~ZZ$METHOD, data=ZZ); mtext("L-skew Distributions")
boxplot(ZZ$T4~ZZ$METHOD, data=ZZ); mtext("L-kurtosis Distributions") #}

\dontrun{
# Data shown in Cheng and Stephens (1989). They have typesetting error on their
# "sigma." Results mu=34.072 and sigma=sqrt(6.874)=2.6218
H590 <- c(27.55, 31.82, 33.74, 34.15, 35.32, 36.78,
          29.89, 32.23, 33.74, 34.44, 35.44, 37.07,
          30.07, 32.28, 33.86, 34.62, 35.61, 37.36,
          30.65, 32.69, 33.86, 34.74, 35.61, 37.36,
          31.23, 32.98, 33.86, 34.74, 35.73, 37.36,
          31.53, 33.28, 34.15, 35.03, 35.90, 40.28,
          31.53, 33.28, 34.15, 35.03, 36.20) # breaking stress MPAx1E6 of carbon block.
mps2par(H590, type="nor", ties="rounding", delta=0.005)$para
mps2par(H590, type="nor", ties="rounding" )$para
mps2par(H590, type="nor", ties="bernstein")$para
#        mu     sigma
# 34.071424  2.622484 # using a slight variant on their eq. 4.1.
# 34.071424  2.622614 # using log10offset=3
# 34.088769  2.690781 # using Bernstein smooth and unaffecting middle of odd tie runs
# The MoranTest show rejection of the Normal distribution at alpha=0.05, with the
# "rounding" and "delta=0.005"" and T=63.8 compared to their result of T=63.1,
# which to be considered that the strategy here is not precisely the same as theirs.}

\dontrun{
# Demonstration of parameter transformation and retransformation
set.seed(9209) # same seed used under mle2par() in parallel example
x <- rlmomco(500, vec2par(c(1,1,3), type="gam")) # 3-p Generalized Gamma
guess <- lmr2par(x, type="gam", p=3) # By providing a 3-p guess the 3-p
# Generalized Gamma will be triggered internally. There are problems passing
# "p" argument to optim() if that function is to pick up the ... argument.
mps2par(x, type="gam", para.int=guess, silent=FALSE,
           ptransf=  function(t) { c(log(t[1]), log(t[2]), t[3])},
           pretransf=function(t) { c(exp(t[1]), exp(t[2]), t[3])})$para
# Reports:       mu     sigma        nu   for some simulated data.
#         0.9997019 1.0135674 3.0259012 }
}
\keyword{distribution (parameters)}
\keyword{parameter (maximum product spacing to)}
\keyword{distribution (utility)}
\keyword{maximum product spacing}
\keyword{maximum product of spacings}
\keyword{parameter (estimation)}
\concept{maximum spacing estimation}
\concept{product of spacing estimation}
\keyword{goodness-of-fit}
\keyword{Moran-Darling statistic}
\keyword{Ensemble list of major parameter estimation methods}
