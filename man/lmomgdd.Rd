\name{lmomgdd}
\alias{lmomgdd}
\title{L-moments of the Gamma Difference Distribution}
\description{
This function estimates the L-moments of the Gamma Difference distribution (Klar, 2015) given the parameters (\eqn{\alpha_1 > 0}, \eqn{\beta_1 > 0}, \eqn{\alpha_2 > 0}, \eqn{\beta_2 > 0}) from \code{\link{pargdd}}. The L-moments in terms of the parameters higher than the mean are complex and numerical methods are required. The mean is

\deqn{\lambda_1 = \frac{\alpha_1}{\beta_1} - \frac{\alpha_2}{\beta_2} \mbox{.}}

The product moments, however, have simple expressions, the variance and skewness, respectively are

\deqn{\sigma^2 = \frac{\alpha_1}{\beta_2^2} + \frac{\alpha_2}{\beta_2^2}\mbox{,}}

and

\deqn{\gamma = \frac{2\bigl(\alpha_1{\beta_2^3} + \alpha_2{\beta_2^2}\bigr)}
                     {\bigl(\alpha_2{\beta_1^2} + \alpha_2{\beta_1^2}\bigr)^{3/2}}\mbox{.}}
}
\usage{
lmomgdd(para, nmom=6, paracheck=TRUE, silent=TRUE, ...)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{nmom}{The number of L-moment to numerically compute for the distribution.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity.}
  \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}.}
  \item{...}{Additional argument to pass.}
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
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomgdd}.}
}
\note{
\bold{Experimental Summer 2024}---For a symmetrical version of the distribution, the relation between \eqn{\tau_4} and \eqn{\tau_6} and other coupling to \eqn{\lambda_2} can be computed using the following recipe:
\preformatted{
  LMR <- NULL
  plotlmrdia46(lmrdia46(), autolegend=TRUE, xleg="topleft")
  for(i in 1:4000) {
    if(length(grep("00$", as.character(i)))) message(i)
    para <- 10^runif(2, min=-3, max=3)
    para <- list(para=c(para[1], para[2], para[1], para[2]), type="gdd")
    lmr  <- lmomgdd(para, nmom=6, subdivisions=200)
    points(lmr$ratios[4], lmr$ratios[6], pch=1, cex=0.8, lwd=0.8)
    LMR <- rbind(LMR, data.frame(A12=para$para[1], B12=para$para[2],
                  L1=lmr$lambdas[1], L2=lmr$lambdas[2], T3=lmr$ratios[3],
                  T4=lmr$ratios[ 4], T5=lmr$ratios[ 3], T6=lmr$ratios[6]))
  }
  LMR <- LMR[complete.cases(LMR),  ]
  LMR <- LMR[abs(  LMR$T3) < 0.01, ]
  LMR <- LMR[order(LMR$T4),        ]
}

We have swept through, hopefully, a sufficiently large span of viable parameter values under a constraint of symmetry. The following recipe continues in post-processing with the goal of producing a polynomial approximation between \eqn{\tau_4} and \eqn{\tau_6} for \code{\link{lmrdia46}}.

\preformatted{
  plotlmrdia46(lmrdia46(), autolegend=TRUE, xleg="topleft")
  points(LMR$T4, LMR$T6, pch=1, cex=0.8, lwd=0.8)
  LM <- lm(T6~I(T4  ) + I(T4^2) + I(T4^3) + I(T4^4) +
              I(T4^5) + I(T4^6) + I(T4^7) + I(T4^8), data=LMR)
  lines(LMR$T4, fitted.values(LM), col="blue", lwd=3)

  res <- residuals(LM)
  plot(fitted.values(LM), res, ylim=c(-0.02, 0.02))
  abline(h=c(-0.002, 0.002), col="red")
  LMRthin <- LMR[abs(res) < 0.002, ]

  LM <- lm(T6~I(T4  ) + I(T4^2) + I(T4^3) + I(T4^4)+
              I(T4^5) + I(T4^6) + I(T4^7) + I(T4^8), data=LMRthin)

  plot(  LMRthin$T4, fitted.values(LM), col="blue", type="l", lwd=3  )
  points(LMRthin$T4, LMRthin$T6,        col="red",  cex=0.4,  lwd=0.5)

  tau4 <- c(lmrdia46()$nor$tau4, 0.1227, 0.123, 0.125, seq(0.13, 1, by=0.01))
  tau6 <- predict(LM, newdata=data.frame(T4=tau4))
  names(tau6) <- NULL
  gddsymt46   <- data.frame(tau4=tau4, tau6=tau6)

  gddsymt46f <- function(t4) { # print(coefficients(LM))
    coe <- c( -0.0969112,    2.1743687, -12.8878580,  47.8931168, -108.0871549,
             156.9200440, -139.5599813,  69.3492358, -14.7052424)
    ix <- seq_len(length(coe))-1
    sapply(t4, function(t) sum(coe[ix+1]*t^ix))
  } # This function is inserted into the lmrdia46() for deployment as symgdd.

  plotlmrdia46(lmrdia46(), autolegend=TRUE, xleg="topleft")
  lines(          tau4, gddsymt46f(tau4), lwd=3, col="deepskyblue3")
  lines(gddsymt46$tau4, gddsymt46$tau6,   lwd=3, col="deepskyblue3")
  legend("bottomright", "Symmetrical Gamma Difference distribution",
         bty="n", cex=0.9, lwd=3, col="deepskyblue3")
}

This is the first known derivation of the relation between these two L-moment ratios for the symmetrical version of this distribution. The quantities recorded in the \code{LMR} data frame in the recipe can be useful for additional study of the quality of numerical implementation of the distribution by \pkg{lmomco}. Next, for purposes of helping parameter estimation for \eqn{\alpha_1 = \alpha_2} and \eqn{\beta_1 = \beta_2} and \eqn{\tau_4}, let us build a polynomial for \eqn{\alpha} estimation from \eqn{\tau_4}:

\preformatted{
  tlogit <- function(x) log(x/(1-x))
  ilogit <- function(x) 1/(1+exp(-x))
  A12l <-    log(LMRthin$A12)
  T4l  <- tlogit(LMRthin$T4 )
  A12p <- exp( approx(T4l, y=A12l, xout=tlogit(tau4))$y )

  plot(  tlogit(tau4), A12p, log="y",     col="blue", type="l", lwd=3  )
  points(LMRthin$T4, LMRthin$A12, col="red",  cex=0.4,  lwd=0.5)
}
}
\author{W.H. Asquith}
\seealso{\code{\link{pargdd}}, \code{\link{cdfgdd}}, \code{\link{pdfgdd}}, \code{\link{quagdd}} }
\examples{
#
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Gamma Difference}
