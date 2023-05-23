\name{mle2par}
\alias{mle2par}
\title{Use Maximum Likelihood to Estimate Parameters of a Distribution}
\description{
This function uses the method of maximum likelihood (MLE) to estimate the parameters of a distribution.

MLE is a straightforward optimization problem that is formed by maximizing the sum of the logarithms of probability densities. Let \eqn{\Theta} represent a vector of parameters for a candidate fit to the specified probability density function \eqn{g(x|\Theta)} and \eqn{x_i} represent the observed data for a sample of size \eqn{n}. The objective function is
\deqn{\mathcal{L}(\Theta) = -\sum_{i=1}^{n} \log\, g(x_i|\Theta)\mbox{,}}
where the \eqn{\Theta} for a maximized \eqn{{-}\mathcal{L}} (note the 2nd negation for the adjective \dQuote{maximized}, \code{optim()} defaults as a minimum optimizer) represents the parameters fit by MLE. The initial parameter estimate by default will be seeded by the method of L-moments.
}
\usage{
mle2par(x, type, para.int=NULL, silent=TRUE, null.on.not.converge=TRUE,
                 ptransf=  function(t) return(t),
                 pretransf=function(t) return(t), ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{type}{Three character (minimum) distribution type (for example, \code{type="gev"}), see \code{\link{dist.list}}.}
  \item{para.int}{Initial parameters as a vector \eqn{\Theta} or as an \pkg{lmomco} parameter \dQuote{object} from say \code{\link{vec2par}}. If a vector is given, then internally \code{\link{vec2par}} is called with distribution equal to \code{type}.}
  \item{silent}{A logical to silence the \code{try()} function wrapping the \code{optim()} function.}
  \item{null.on.not.converge}{A logical to trigging simple return of \code{NULL} if the \code{optim()} function returns a nonzero convergence status.}
  \item{ptransf}{An optional parameter transformation function (see \bold{Examples}) that is useful to guide the optimization run. For example, suppose the first parameter of a three parameter distribution resides in the positive domain, then \cr \code{ptransf(t) = }
\code{function(t) c(log(t[1]), t[2], t[3])}.}
  \item{pretransf}{An optional parameter retransformation function (see \bold{Examples}) that is useful to guide the optimization run. For example, suppose the first parameter of a three parameter distribution resides in the positive domain, then \cr \code{pretransf(t) = } \code{function(t) c(exp(t[1]), t[2], t[3])}.}
  \item{...}{Additional arguments for the \code{optim()} function and other uses.}
}
\value{
  An \R \code{list} is returned.  This list should contain at least the following items, but some distributions such as the \code{revgum} have extra.
  \item{type}{The type of distribution in three character (minimum) format.}
  \item{para}{The parameters of the distribution.}
  \item{source}{Attribute specifying source of the parameters.}
  \item{AIC}{The Akaike information criterion (AIC).}
  \item{optim}{The returned \code{list} of the \code{optim()} function.}
}
\note{
During the optimization process, the function requires evaluation at the initial parameters. The following error rarely will be seen:
\preformatted{
  Error in optim(para.int$para, afunc) :
    function cannot be evaluated at initial parameters
}
if \code{Inf} is returned on first call to the objective function. The \code{silent} by default though will silence this error.  Alternative starting parameters might help.  This function is not built around subordinate control functions to say keep the parameters within distribution-specific bounds.  However, in practice, the L-moment estimates should already be fairly close and the optimizer can take it from there.  More sophisticated MLE for many distributions is widely available in other \R packages. The \pkg{lmomco} package uses its own probability density functions.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmom2par}}, \code{\link{mps2par}}, \code{\link{tlmr2par}}}
\examples{
\dontrun{
# This example might fail on mle2par() or mps2par() depending on the values
# that stem from the simulation. Trapping for a NULL return is not made here.
father <- vec2par(c(37,25,114), type="st3"); FF <- nonexceeds(); qFF <- qnorm(FF)
X <- rlmomco(78, father) # rerun if MLE and MPS fail to get a solution
plot(qFF,  qlmomco(FF, father), type="l", xlim=c(-3,3),
     xlab="STANDARD NORMAL VARIATE", ylab="QUANTILE") # parent (black)
lines(qFF, qlmomco(FF, lmr2par(X, type="gev")), col="red"  ) # L-moments (red)
lines(qFF, qlmomco(FF, mps2par(X, type="gev")), col="green") #     MPS (green)
lines(qFF, qlmomco(FF, mle2par(X, type="gev")), col="blue" ) #     MLE  (blue)
points(qnorm(pp(X)), sort(X)) # the simulated data}

\dontrun{
# REFLECTION SYMMETRY
set.seed(451)
X <- rlmomco(78, vec2par(c(2.12, 0.5, 0.6), type="pe3"))
# MLE and MPS are almost reflection symmetric, but L-moments always are.
mle2par( X, type="pe3")$para #  2.1796827 0.4858027  0.7062808
mle2par(-X, type="pe3")$para # -2.1796656 0.4857890 -0.7063917
mps2par( X, type="pe3")$para #  2.1867551 0.5135882  0.6975195
mps2par(-X, type="pe3")$para # -2.1868252 0.5137325 -0.6978034
parpe3(lmoms( X))$para       #  2.1796630 0.4845216  0.7928016
parpe3(lmoms(-X))$para       # -2.1796630 0.4845216 -0.7928016}

\dontrun{
Ks <- seq(-1,+1,by=0.02); n <- 100; MLE <- MPS <- rep(NA, length(Ks))
for(i in 1:length(Ks)) {
  sdat   <- rlmomco(n, vec2par(c(1,0.2,Ks[i]), type="pe3"))
  mle    <- mle2par(sdat, type="pe3")$para[3]
  mps    <- mps2par(sdat, type="pe3")$para[3]
  MLE[i] <- ifelse(is.null(mle), NA, mle) # A couple of failures expected as NA's.
  MPS[i] <- ifelse(is.null(mps), NA, mps) # Some amount fewer failures than MLE.
}
plot( MLE, MPS, xlab="SKEWNESS BY MLE", ylab="SKEWNESS BY MPS")#}

\dontrun{
# Demonstration of parameter transformation and retransformation
set.seed(9209) # same seed used under mps2par() in parallel example
x <- rlmomco(500, vec2par(c(1,1,3), type="gam")) # 3-p Generalized Gamma
guess <- lmr2par(x, type="gam", p=3) # By providing a 3-p guess the 3-p
# Generalized Gamma will be triggered internally. There are problems passing
# "p" argument to optim() if that function is to pick up the ... argument.
mle2par(x, type="gam", para.int=guess, silent=FALSE,
           ptransf=  function(t) { c(log(t[1]), log(t[2]), t[3])},
           pretransf=function(t) { c(exp(t[1]), exp(t[2]), t[3])})$para
# Reports:       mu     sigma        nu   for some simulated data.
#         1.0341269 0.9731455 3.2727218 }

\dontrun{
# Demonstration of parameter estimation with tails of density zero, which
# are intercepted internally to maintain finiteness. We explore the height
# distribution for male cats of the cats dataset from the MASS package and
# fit the generalized lambda. The log-likelihood is shown by silent=FALSE
# to see that the algorithm converges slowly. and shows how we can control
# the relative tolerance of the optim() function as shown below and
# investigate the convergence by reviewing the five fits to the data.
FF <- nonexceeds(sig6=TRUE)
library(MASS); data(cats); x <- cats$Hwt[cats$Sex == "M"]
p2 <- mle2par(x, type="gld", silent=FALSE, control=list(reltol=1E-2))
p3 <- mle2par(x, type="gld", silent=FALSE, control=list(reltol=1E-3))
p4 <- mle2par(x, type="gld", silent=FALSE, control=list(reltol=1E-4))
p5 <- mle2par(x, type="gld", silent=FALSE, control=list(reltol=1E-5))
p6 <- mle2par(x, type="gld", silent=FALSE, control=list(reltol=1E-6))
plot( FF,  quagld(FF, p2), type="l", col="black")
points(pp(x), sort(x), lwd=0.6, col=grey(0.6))
lines(FF,  quagld(FF, p3), col="red"    )
lines(FF, par2qua(FF, p4), col="green"  )
lines(FF,  quagld(FF, p5), col="blue"   )
lines(FF, par2qua(FF, p6), col="magenta") #}
}
\keyword{distribution (parameters)}
\keyword{parameter (maximum likelihood to)}
\keyword{distribution (utility)}
\keyword{parameter (estimation)}
\keyword{Ensemble list of major parameter estimation methods}
