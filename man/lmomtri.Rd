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
The higher L-moments are even more ponderous and simpler expressions for the L-moment ratios appear elusive. Bounds for \eqn{\tau_3} and \eqn{\tau_4} are \eqn{|\tau_3| \le 0.14285710} and \eqn{0.04757138 < \tau_4 < 0.09013605}. An approximation for \eqn{\tau_4} is
\deqn{\tau_4 = 0.09012180 - 1.777361\tau_3^2 - 17.89864\tau_3^4 + 920.4924\tau_3^6 - 37793.50\tau_3^8 \mbox{,}}
where the residual standard error is \eqn{{<}1.750\times 10^{-5}} and the absolute value of the maximum residual is \eqn{<9.338\times 10^{-5}}. The L-moments of the Symmetrical Triangular distribution for \eqn{\tau_3 = 0} are considered by Nagaraja (2013) and therein for a symmetric triangular distribution having \eqn{\lambda_1 = 0.5} then \eqn{\lambda_4 = 0.0105} and \eqn{\tau_4 = 0.09}. These L-kurtosis values agree with results of this function that are based on the \code{\link{theoLmoms.max.ostat}} function. The 4th and 5th L-moments \eqn{\lambda_4} and \eqn{\lambda_5}, respectively, are computed using expectations of order statistic maxima (\code{\link{expect.max.ostat}}) and are defined (Asquith, 2011, p. 95) as
\deqn{\lambda_4 = 5\mathrm{E}[X_{4:4}] - 10\mathrm{E}[X_{3:3}] +  6\mathrm{E}[X_{2:2}] - \mathrm{E}[X_{1:1}]}
and
\deqn{\lambda_5 = 14\mathrm{E}[X_{5:5}] - 35\mathrm{E}[X_{4:4}] + 30\mathrm{E}[X_{3:3}] - 10\mathrm{E}[X_{2:2}] + \mathrm{E}[X_{1:1}]\mbox{.}}
These expressions are solved using the \code{\link{expect.max.ostat}} function to compute the \eqn{\mathrm{E}[X_{r:r}]}.

For the symmetrical case of \eqn{\omega = (\psi + \nu)/2}, then 
\deqn{\lambda_1 = \frac{(\nu+\psi)}{2}\mbox{\ and}}
\deqn{\lambda_2 = \frac{7}{60}\biggl[\psi - \nu\biggr]\mbox{,}}
which might be useful for initial parameter estimation through
\deqn{\psi = \lambda_1 + \frac{30}{7}\lambda_2 \mbox{\ and}}
\deqn{\nu  = \lambda_1 - \frac{30}{7}\lambda_2 \mbox{.}}
}
\usage{
lmomtri(para, paracheck=TRUE, nmom=c("3", "5"))
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity. Overriding of this check might help in numerical optimization of parameters for modes near either the minimum or maximum. The argument here makes code base within \code{\link{partri}} a little shorter.}
  \item{nmom}{The L-moments greater the \eqn{r > 3} require numerical integration using the expectations of the maxima order statistics of the fitted distribution. If this argument is set to \code{"3"} then executation of \code{lmomtri} is stopped at \eqn{r = 3} and the first three L-moments returned, otherwise the 4th and 5th L-moments are computed.}
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
  \item{E33err}{A percent error between the expectation of the \eqn{X_{3:3}} order statistic by analytical expression versus a theoretical by numerical integration using the \cr \code{\link{expect.max.ostat}} function. This will be \code{NA} if \code{nmom == "3"}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomtri}.}
}
\note{
The expression for \eqn{\tau_4} in terms of \eqn{\tau_3} is
\preformatted{
  "tau4tri" <- function(t3) {
     t3[t3 < -0.14285710 | t3 >  0.14285710] <- NA
     b <- 0.09012180
     a <- c(0, -1.777361, 0, -17.89864, 0,  920.4924, 0, -37793.50)
     t4 <- b + a[2]*t3^2 + a[4]*t3^4 + a[6]*t3^6 + a[8]*t3^8
     return(t4)
  }
}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Nagaraja, H.N., 2013, Moments of order statistics and L-moments for the symmetric triangular distribution: Statistics and Probability Letters, v. 83, no. 10, pp. 2357--2363.
}
\author{ W.H. Asquith}
\seealso{\code{\link{partri}}, \code{\link{cdftri}}, \code{\link{pdftri}}, \code{\link{quatri}} }
\examples{
lmr <- lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52))
lmr
lmomtri(partri(lmr), nmom="5")

par <- vec2par(c(-405, -390, -102), type="tri")
lmomtri(par, nmom="5")$lambdas
# -299           39.4495050    5.5670228    1.9317914    0.8007511
theoLmoms.max.ostat(para=par, qua=quatri, nmom=5)$lambdas
# -299.0000126   39.4494885    5.5670486    1.9318732    0.8002989
# The -299 is the correct by exact solution as are 39.4495050 and 5.5670228, the 4th and
# 5th L-moments diverge from theoLmoms.max.ostat() because the exact solutions and not
# numerical integration of the quantile function was used for E11, E22, and E33.
# So although E44 and E55 come from expect.max.ostat() within both lmomtri() and
# theoLmoms.max.ostat(), the Lambda4 and Lambda5 are not the same because the E11, E22,
# and E33 values are different.

\dontrun{
# At extreme limit of Tau3 for the triangular distribution, L-moment ratio diagram
# shows convergence to the trajectory of the Generalized Pareto distribution.
"tau4tri" <- function(t3) { t3[t3 < -0.14285710 | t3 >  0.14285710] <- NA
   b <- 0.09012180; a <- c(0, -1.777361, 0, -17.89864, 0,  920.4924, 0, -37793.50)
   t4 <- b + a[2]*t3^2 + a[4]*t3^4 + a[6]*t3^6 + a[8]*t3^8; return(t4)
}
F <- seq(0,1, by=0.001)
lmr  <- vec2lmom(c(10,9,0.142857, tau4tri(0.142857)))
parA <- partri(lmr); parB <- pargpa(lmr)
xA <- qlmomco(F,  parA); xB <- qlmomco(F, parB); x <- sort(unique(c(xA,xB)))
plot(x,  pdftri(x,parA), type="l", col=8, lwd=4) # Compare Asym. Tri. to 
lines(x, pdfgpa(x,parB),           col=2)        # Gen. Pareto
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}
