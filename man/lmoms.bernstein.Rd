\encoding{utf8}
\name{lmoms.bernstein}
\alias{lmoms.bernstein}
\title{Numerically Integrated L-moments of Smoothed Quantiles from Bernstein or Kantorovich Polynomials }
\description{
Compute the L-moment by numerical integration of the smoothed quantiles from Bernstein or Kantorovich polynomials (see \code{\link{dat2bernqua}}). Letting \eqn{\tilde{X}_n(F)} be the smoothed quantile function for nonexceedance probability \eqn{F} for a sample of size \eqn{n}, from Asquith (2011) the first five L-moments in terms of quantile function integration are
 \deqn{\lambda_1 = \int_0^1 \tilde{X}_n(F)\;\mathrm{d}F \mbox{,}}
 \deqn{\lambda_2 = \int_0^1 \tilde{X}_n(F)\times(2F - 1)\;\mathrm{d}F\mbox{,}}
 \deqn{\lambda_3 = \int_0^1 \tilde{X}_n(F)\times(6F^2 - 6F + 1)\;\mathrm{d}F\mbox{,}}
 \deqn{\lambda_4 = \int_0^1 \tilde{X}_n(F)\times(20F^3 - 30F^2 + 12F - 1)\;\mathrm{d}F\mbox{, and}}
 \deqn{\lambda_5 = \int_0^1 \tilde{X}_n(F)\times(70F^4 - 140F^3 + 90F^2 - 20F + 1)\;\mathrm{d}F\mbox{.}}
}
\usage{
lmoms.bernstein(x, bern.control=NULL,
                   poly.type=c("Bernstein", "Kantorovich", "Cheng"),
                   bound.type=c("none", "sd", "Carv", "either"),
                   fix.lower=NULL, fix.upper=NULL, p=0.05)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{bern.control}{A \code{list} that holds \code{poly.type}, \code{bound.type}, \code{fix.lower}, and \code{fix.upper}. And this list will supersede the respective
values provided as separate arguments.}
  \item{poly.type}{Same argument as for \code{\link{dat2bernqua}}.}
  \item{bound.type}{Same argument as for \code{\link{dat2bernqua}}.}
  \item{fix.lower}{Same argument as for \code{\link{dat2bernqua}}.}
  \item{fix.upper}{Same argument as for \code{\link{dat2bernqua}}.}
  \item{p}{The \dQuote{p-factor} is the same argument as for \code{\link{dat2bernqua}}.}
}
\value{
  An \R \code{vector} is returned.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{dat2bernqua}}, \code{\link{pfactor.bernstein}}, \code{\link{lmoms}} }
\examples{
\dontrun{
X <- exp(rnorm(100))
lmoms.bernstein(X)$ratios
lmoms.bernstein(X, fix.lower=0)$ratios
lmoms.bernstein(X, fix.lower=0, bound.type="sd")$ratios
lmoms.bernstein(X, fix.lower=0, bound.type="Carv")$ratios
lmoms(X)$ratios

lmoms.bernstein(X, poly.type="Kantorovich")$ratios
lmoms.bernstein(X, fix.lower=0, poly.type="Kantorovich")$ratios
lmoms.bernstein(X, fix.lower=0, bound.type="sd", poly.type="Kantorovich")$ratios
lmoms.bernstein(X, fix.lower=0, bound.type="Carv", poly.type="Kantorovich")$ratios
lmoms(X)$ratios
}

\dontrun{
lmr <- vec2lmom(c(1,.2,.3))
par <- lmom2par(lmr, type="gev")
lmr <- lmorph(par2lmom(par))
lmT <- c(lmr$lambdas[1:2], lmr$ratios[3:5])
ns  <- 200; nsim <- 1000; empty <- rep(NA, nsim)

sink("ChengLmomentTest.txt")
cat(c("N errmeanA  errlscaleA  errtau3A  errtau4A  errtau5A",
        "errmeanB  errlscaleB  errtau3B  errtau4B  errtau5B\n"))
for(n in 1:ns) {
   message(n);
   SIM <- data.frame(errmeanA=empty, errlscaleA=empty,   errtau3A=empty, errtau4A=empty,
                     errtau5A=empty,   errmeanB=empty, errlscaleB=empty, errtau3B=empty,
                     errtau4B=empty,   errtau5B=empty)
   for(i in 1:nsim) {
      X <- rlmomco(30, par)
      lmrA <- lmoms(X)
      lmA <- c(lmrA$lambdas[1:2], lmrA$ratios[3:5])
      lmrB <- lmoms.bernstein(X, poly.type="Cheng")
      lmB <- c(lmrB$lambdas[1:2], lmrB$ratios[3:5])
      EA <- lmA - lmT; EB <- lmB - lmT
      SIM[i,] <- c(EA,EB)
   }
   MeanErr <- sapply(1:length(SIM[1,]), function(x) { return(mean(SIM[,x])) })
   line <- paste(c(n, round(MeanErr, digits=6), "\n"), sep=" ")
   cat(line)
}
sink()
}
}
\keyword{L-moment (theoretical)}
\keyword{smoothing (Bernstein-Kantorovich)}

