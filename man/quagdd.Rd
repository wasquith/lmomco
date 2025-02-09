\name{quagdd}
\alias{quagdd}
\title{Quantile Function of the Gamma Difference Distribution}
\description{
This function computes the quantiles of the Gamma Difference distribution (Klar, 2015) given parameters (\eqn{\alpha_1 > 0}, \eqn{\beta_1 > 0}, \eqn{\alpha_2 > 0}, \eqn{\beta_2 > 0}) computed by \code{\link{pargdd}}. The quantile function requires numerical rooting of the cumulative distribution function \code{\link{cdfgdd}}.
}
\usage{
quagdd(f, para, paracheck=TRUE, silent=TRUE, ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargdd}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity.}
  \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}.}
  \item{...}{Additional arguments to pass.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\author{W.H. Asquith}
\references{
Klar, B., 2015, A note on gamma difference distributions: Journal of Statistical Computation and Simulation v. 85, no. 18, pp. 1--8, \doi{10.1080/00949655.2014.996566}.
}
\seealso{\code{\link{cdfgdd}}, \code{\link{pdfgdd}}, \code{\link{lmomgdd}}, \code{\link{pargdd}}}
\examples{
\dontrun{
  para <- list(para=c(3, 0.1, 0.1, 4), type="gdd")
  quagdd(0.5, para) # [1] 26.71568 }

\dontrun{
  p <- c(3, 1, 0.2, 2)
  NEP  <- seq(0.001, 0.999, by=0.001)
  para <- list(para=p, type="gdd")
  F1 <- runif(10000); F2 <- runif(10000)
  XX  <- sort(qgamma(F1, p[1], p[2]) - qgamma(F2, p[3], p[4])); FF  <- pp(XX)
  plot(NEP, quagdd(NEP, para), type="l", col=grey(0.8), lwd=6,
       xlab="Nonexceedance probability", ylab="Gamma difference quantile")
  lines(FF, XX, col="red") # }

\dontrun{
  para <- list(para=c(3, 2, 0.2, 2), type="gdd")
  nsam <- 50; nsim <- 10
  F1 <- runif(1000); F2 <- runif(1000)
  plot(c(1,4), c(0, 2), type="n", xlim=c(0.5, 4.5),
       xlab="Lmoment order", ylab="Lambda value")
  for(i in seq_len(nsim)) {
    X <- quagdd(runif(nsam), para)
    afunc <- function(par, lmr=NA) { p <- exp(par)
      tlmr <- pwm2lmom(pwm(qgamma(F1, p[1], p[2]) -
                           qgamma(F2, p[3], p[4])))
      sum((lmr$lambdas[1:4] - tlmr$lambdas[1:4])^2)
    }
    slmr <- lmoms(X, nmom=4); init.para <- c(0, 0, 0, 0)
    sara <- NULL
    try( sara <- optim( init.para, afunc, lmr=slmr ) )
    if(is.null(sara)) next
    sara$para <- exp(sara$par); sara$type <- "gdd"
    mara <- pargdd(slmr)
    points(1:4+0.1, lmomgdd(mara)$lambdas[1:4], col="red" )
    points(1:4-0.1, lmomgdd(sara)$lambdas[1:4], col="blue")
    print(lmomgdd(mara)$lambdas[1:4])
    print(lmomgdd(sara)$lambdas[1:4])
  }
  lines( 1:4, lmomgdd(para)$lambdas[1:4], col="darkgreen")
  points(1:4, lmomgdd(para)$lambdas[1:4], pch=16, col="darkgreen") # }
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Gamma Difference}
