\name{pargam}
\alias{pargam}
\title{Estimate the Parameters of the Gamma Distribution}
\description{
This function estimates the parameters of the Gamma distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. Both the two-parameter Gamma and three-parameter Generalized Gamma distributions are supported based on the desired choice of the user, and numerical-hybrid methods are required. The \code{\link{pdfgam}} documentation provides further details.
}
\usage{
pargam(lmom, p=c("2", "3"), checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{A L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{p}{The number of parameters to estimate for the 2-p Gamma or 3-p Generalized Gamma.} 
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gam}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{pargam}.}
}
\note{
The two-parameter Gamma is supported by Hosking's code-based approximations to avoid direct numerical techniques. The three-parameter version is based on a dual approach to parameter optimization. The \eqn{\log(\sigma)} and \eqn{\sqrt{\log(\lambda_1/\lambda_2)}} conveniently has a relatively narrow range of variation. A polynomial approximation to provide a first estimate of \eqn{\sigma} (named \eqn{\sigma'}) is used through the \code{optim()} function to isolated the best estimates of \eqn{\mu'} and \eqn{\nu'} of the distribution holding \eqn{\sigma} constant at \eqn{\sigma = \sigma'}---a 2D approach is thus involved. Then, the initial parameter for a second three-dimensional optimization is made using the initial parameter estimates as the tuple \eqn{\mu', \sigma', \nu'}. This 2D approach seems more robust and effectively canvases more of the Generalized Gamma parameter domain, though a doubled-optimization is not quite as fast as a direct 3D optimization. The following code was used to derive the polynomial coefficients used for the first approximation of \eqn{sigma'}:
\preformatted{
  nsim <- 10000; mu <- sig <- nu <- l1 <- l2 <- t3 <- t4 <- rep(NA, nsim)
  for(i in 1:nsim) {
    m <- exp(runif(1, min=-4, max=4)); s <- exp(runif(1, min=-8, max=8))
    n <- runif(1, min=-14, max=14); mu[i] <- m; sig[i] <- s; nu[i] <- n
    para <- vec2par(c(m,s,n), type="gam"); lmr <- lmomgam(para)
    if(is.null(lmr)) next
    lam <- lmr$lambdas[1:2]; rat <- lmr$ratios[3:4]
    l1[i]<-lam[1]; l2[i]<-lam[2];t3[i]<-rat[1]; t4[i]<-rat[2]
  }
  ZZ <- data.frame(mu=mu, sig=sig, nu=nu, l1=l1, l2=l2, t3=t3, t4=t4)
  ZZ$ETA <- sqrt(log(ZZ$l1/ZZ$l2)); ZZ <- ZZ[complete.cases(ZZ), ]
  ix <- 1:length(ZZ$ETA);  ix <- ix[(ZZ$ETA < 0.025 & log(ZZ$sig) < 1)]
  ZZ <- ZZ[-ix,]
  with(ZZ, plot(ETA, log(sig), xlim=c(0,4), ylim=c(-8,8)))
  LM <- lm(log(sig)~
           I(1/ETA^1)+I(1/ETA^2)+I(1/ETA^3)+I(1/ETA^4)+I(1/ETA^5)+
               ETA   +I(  ETA^2)+I(  ETA^3)+I(  ETA^4)+I(  ETA^5), data=ZZ)
  ETA <- seq(0,4,by=0.002) # so the line of fit can be seen
  lines(ETA, predict(LM, newdata=list(ETA=ETA)), col=2)
  The.Coefficients.In.pargam.Function <- LM$coefficients
}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomgam}},
\code{\link{cdfgam}}, \code{\link{pdfgam}}, \code{\link{quagam}}
}
\examples{
pargam(lmoms(abs(rnorm(20, mean=10))))

\dontrun{
pargam(lmomgam(vec2par(c(0.3,0.4,+1.2), type="gam")), p=3)$para
pargam(lmomgam(vec2par(c(0.3,0.4,-1.2), type="gam")), p=3)$para
#        mu      sigma         nu 
# 0.2999994  0.3999990  1.1999696
# 0.2999994  0.4000020 -1.2000567}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Gamma}
\keyword{Distribution: Generalized Gamma}
