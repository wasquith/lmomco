\name{lmomtri}
\alias{lmomtri}
\title{L-moments of the Asymmetric Triangular Distribution}
\description{
This function estimates the L-moments of the Asymmetric Triangular distribution given the parameters (\eqn{\nu}, \eqn{\omega}, and \eqn{\psi}) from \code{\link{partri}}. The first three L-moments in terms of the parameters are
\deqn{\lambda_1 = \frac{(\nu+\omega+\psi)}{3}\mbox{,}}
\deqn{\lambda_2 = \frac{1}{15}\biggl[\frac{(\nu-\omega)^2}{(\psi-\nu)^{\phantom{1}}} - (\nu+\omega) + 2\psi\biggr] \mbox{, and}}
\deqn{\lambda_3 = G + H_1 + H_2 + J \mbox{,}}
where \eqn{G} is dependent on the integral definining the L-moments in terms of the quantile function (Asquith, 2011, p. 92) with limits of integration of \eqn{[0,P]}, \eqn{H_1} and \eqn{H_2} are dependent on the integral defining the L-moment in terms of the quantile function with limits of integration of \eqn{[P,1]}, and \eqn{J} is dependent on the \eqn{\lambda_2} and \eqn{\lambda_1}. Finally, the variables \eqn{G}, \eqn{H_1}, \eqn{H_2}, and \eqn{J} are
\deqn{G = \frac{2}{7}\frac{(\nu+6\omega)(\omega-\nu)^3}{(\psi-\nu)^3}\mbox{,}}
\deqn{H_1 = \frac{12}{7}\frac{(\omega-\psi)^4}{(\nu-\psi)^3} - 2\psi\frac{(\nu-\omega)^3}{(\nu-\psi)^3} + 2\psi\mbox{,}}
\deqn{H_2 = \frac{4}{5}\frac{(5\nu-6\omega+\psi)(\omega-\psi)^2}{(\nu-\psi)^2}\mbox{, and}}
\deqn{J = -\frac{1}{15}\biggl[\frac{3(\nu-\omega)^2}{(\psi-\nu)} + 7(\nu+\omega) + 16\psi\biggl]\mbox{.}}
The higher L-moments are even more ponderous and simple expressions for the L-moment ratios appear elusive.


The L-moments of the Symmetrical Triangular distribution for \eqn{\tau_3 = 0} are considered by Nagaraja (2013) and therein for such a distribution having \eqn{\lambda_1 = 0.5} then \eqn{\lambda_4 = 0.0105} and \eqn{\tau_4 = 0.09}. These L-kurtosis values agree with results of this function that are based on the \code{\link{theoLmoms.max.ostat}} function. Continuing, for the symmetrical case of \eqn{\omega = (\psi + \nu)/2}, then
\deqn{\lambda_1 = \frac{(\nu+\psi)}{2}\mbox{\ and\ } \lambda_2 = \frac{7}{60}\biggl[\psi - \nu\biggr]\mbox{,}}
which might be useful for initial parameter estimation through
\deqn{\psi = \lambda_1 + \frac{30}{7}\lambda_2 \mbox{\ and\ } \nu  = \lambda_1 - \frac{30}{7}\lambda_2 \mbox{.}}
}
\usage{
lmomtri(para, paracheck=TRUE)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity. Overriding of this check might help in numerical optimization of parameters for modes near either the minimum or maximum. The argument here makes code base within \code{\link{partri}} a little shorter.}
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
  \item{tau3poly}{A cross check on the polynomial approximation scheme, and \code{tau3poly} is \eqn{\tau_3} from the polynomial and not the direct solution show above, and \code{tau3poly} is not used to populate the \code{ratios[3]}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomtri}.}
}
\note{
Study of some simulation output of \code{inst/newcoes/makelmomtri.R}, bounds for \eqn{\tau_3} and \eqn{\tau_4} are \eqn{|\tau_3| \le 0.14285717} or \eqn{\tau_3 = 1/7} (numerically at tolerance 1E-7) and \eqn{1/21 < \tau_4 < 233 / 2585} (numerically at tolerances 1E-7). An approximation for \eqn{\tau_4} is

\deqn{\tau_4 = 0.090128216 - 1.7893336\tau_3^2 - 14.486666\tau_3^4 + 619.08286\tau_3^6 - 29658.04718201972\tau_3^8 \mbox{,}}

where the residual standard error is \eqn{{<}1.750\times 10^{-5}} and the absolute value of the maximum residual is \eqn{<9.338\times 10^{-5}}. Lastly, an expression for \eqn{\tau_4} in terms of \eqn{\tau_3} by \code{\link{theoLmoms}} numerical integration is
\preformatted{
  "tau4tri" <- function(t3) {
     t3[t3 < -0.1428571650715543 | t3 >  0.1428571650715543] <- NA
     b <- 9.012821605998006e-02
     a <- c(0, -1.789333565397665e+00, 0,
               -1.448666575025029e+01, 0,
               66.190828635274557e+02, 0, -2.965804718201972e+04)
     t4 <- b + a[2]*t3^2 + a[4]*t3^4 + a[6]*t3^6 + a[8]*t3^8
     return(t4)
  }
}

\bold{Deprecated May 2026 in favor of built-in polynomial approximations}---The 4th and 5th L-moments \eqn{\lambda_4} and \eqn{\lambda_5}, respectively, are computed using expectations of order statistic maxima\cr (\code{\link{expect.max.ostat}}) and are defined (Asquith, 2011, p. 95) as
\deqn{\lambda_4 = 5\mathrm{E}[X_{4:4}] - 10\mathrm{E}[X_{3:3}] +  6\mathrm{E}[X_{2:2}] - \mathrm{E}[X_{1:1}]}
and
\deqn{\lambda_5 = 14\mathrm{E}[X_{5:5}] - 35\mathrm{E}[X_{4:4}] + 30\mathrm{E}[X_{3:3}] - 10\mathrm{E}[X_{2:2}] + \mathrm{E}[X_{1:1}]\mbox{.}}
These expressions are solved using the \code{\link{expect.max.ostat}} function to compute the \eqn{\mathrm{E}[X_{r:r}]}. The historical logic is retained in \code{lmomtri} but is below the \code{return()}.
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Nagaraja, H.N., 2013, Moments of order statistics and L-moments for the symmetric triangular distribution: Statistics and Probability Letters, v. 83, no. 10, pp. 2357--2363.
}
\author{W.H. Asquith}
\seealso{\code{\link{partri}}, \code{\link{cdftri}}, \code{\link{pdftri}}, \code{\link{quatri}} }
\examples{
lmr <- lmoms(c(46, 70, 59, 36, 48, 46, 63, 35, 52))
print(                lmr)
print( lmomtri(partri(lmr)) )

\dontrun{
  par <- vec2par(c(-405, -390, -102), type="tri")
  lmomtri(par)$lambdas
  # -299.0000000   39.4495050    5.5670228    1.9317517    0.8007754    0.5057529
  theoLmoms.max.ostat(para=par, qua=quatri, nmom=6)$lambdas
  # -299.0000126   39.4494885    5.5670486    1.9318732    0.8002989    0.5070410
  # The -299 is the correct by exact solution as are 39.4495050 and 5.5670228, the
  # 4-6th L-moments slightly diverge from theoLmoms.max.ostat() and its numerical
  # integration relative the theoLmoms() used to form basis on which the built-in
  # polynomial expressions for the 3-6 Lmoment ratios were made. }

\dontrun{
   # At extreme limit of Tau3 for the triangular distribution, L-moment ratio diagram
   # shows convergence to the trajectory of the Generalized Pareto distribution.
  FF   <- seq(0, 1, by=0.001)
  lmr  <- vec2lmom( c(10, 9, 1/7, 233 / 2585) )
  parA <- partri(lmr); parB <- pargpa(lmr)
  xA   <- qlmomco(FF,  parA); xB <- qlmomco(FF, parB); x <- sort(unique(c(xA, xB)))
  plot(x,  pdftri(x, parA), type="l", col=8, lwd=4) # Compare Asym. Tri. to
  lines(x, pdfgpa(x, parB),           col=2)        # Gen. Pareto }

\dontrun{
  # Alternative view of the higher L-moment ratios for the triangular distribution.
  df <- NULL
  for(i in seq_len(2000)) {
    l <- -runif(1, min=0, max=1000); r <- runif(1, min=0, max=1000)
    lmr <- theoLmoms(vec2par(c(l,0, r), type="tri"), nmom=6, subdivisions=200)
    t3 <- lmr$ratios[3]; t4 <- lmr$ratios[4]; t5 <- lmr$ratios[5]; t6 <- lmr$ratios[6]
    df <- rbind(df, data.frame(l=l, r=r, rho=abs(l/r), t3=t3, t4=t4, t5=t5, t6=t6))
  }
  df <- df[order(df$t3),]

  rho <- 10^seq(-3, 3, by=0.01)
  t3x <- approx(log10(df$rho), df$t3, xout=log10(rho), rule=2)$y
  plot(rho, t3x, log="x", type="l")
  t4x <- approx(log10(df$rho), df$t4, xout=log10(rho), rule=2)$y
  plot(rho, t4x, log="x", type="l")
  t5x <- approx(log10(df$rho), df$t5, xout=log10(rho), rule=2)$y
  plot(rho, t5x, log="x", type="l")
  t6x <- approx(log10(df$rho), df$t6, xout=log10(rho), rule=2)$y
  plot(rho, t6x, log="x", type="l")

  dfo <- data.frame(rho, t3=t3x, t4=t4x, t5=t5x, t6=t6x)
  t3s <- unique(c(seq(min(dfo$t3), max(dfo$t3), by=0.001), max(dfo$t3)))
  t4s <- approx(dfo$t3, dfo$t4, xout=t3s, rule=2)$y
  plotlmrdia(lmrdia(), xlim=c(-0.2, .2), ylim=c(0,0.2))
  lines(t3s, t4s) # }
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}
