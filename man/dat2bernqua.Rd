\encoding{utf8}
\name{dat2bernqua}
\alias{dat2bernqua}
\title{Observed Data to Empirical Quantiles through Bernstein or Kantorovich Polynomials }
\description{
The empirical quantile function can be \dQuote{smoothed} (\enc{Hernández-Maldonado}{Hernandez-Maldonado} and others, 2012, p. 114) through the Kantorovich polynomial (\enc{Muñoz-Pérez}{Munoz-Perez} and \enc{Fernández-Palacín}{Fernandez-Palacin}, 1987) for the sample order statistics \eqn{x_{k:n}} for a sample of size \eqn{n} by
\deqn{\tilde{X}_n(F) = \frac{1}{2}\sum_{k=0}^n (x_{k:n} + x_{(k+1):n}) {n \choose k} F^k (1-F)^{n-k}\mbox{,}}
where \eqn{F} is nonexceedance probability, and \eqn{(n\:k)} are the binomial coefficients from the \R function \code{choose()}, and the special situations for \eqn{k=0} and \eqn{k=n} are described within the Note section. The form for the Bernstein polynomial is
\deqn{\tilde{X}_n(F) = \sum_{k=0}^{n+1} (x_{k:n}) {n+1 \choose k} F^k (1-F)^{n+1-k}\mbox{.}} There are subtle differences between the two and \code{\link{dat2bernqua}} function supports each. Readers are also directed to the \emph{Special Attention} section.

Turnbull and Ghosh (2014) consider through the direction of a referee and recommendation of \eqn{p=0.05} by that referee (and credit to ideas by de Carvalho [2012]) that the support of the probability density function for the Turnbull and Ghosh (2014) study of Bernstein polynomials can be computed letting \eqn{\alpha = (1 - p)^{-2} - 1} by
\deqn{ \biggl(x_{1:n} - (x_{2:n} - x_{1:n})/\alpha,\: x_{n:n} + (x_{n:n} - x_{n-1:n})/\alpha\biggr)\mbox{,}}
for the minimum and maximum, respectively. Evidently, the original support considered by Turnbull and Ghosh (2014) was
\deqn{ \biggl(x_{1:n} - \lambda_2\sqrt{\pi/n},\: x_{n:n} +  \lambda_2\sqrt{\pi/n}\biggr)\mbox{,}}
for the minimum and maximum, respectively and where the standard deviation is estimated in the function using the 2nd L-moment as \eqn{s = \lambda\sqrt{\pi}}.

The \eqn{p} is referred to by this author as the \dQuote{p-factor} this value has great influence in the estimated support of the distribution and therefore distal-tail estimation or performance is sensitive to the value for \eqn{p}. General exploratory analysis suggests that the \eqn{p} can be optimized based on information external or internal to the data for shape restrained smoothing. For example, an analyst might have external information as to the expected L-skew of the phenomenon being studied or could use the sample L-skew of the data (internal information) for shape restraint (see \code{\link{pfactor.bernstein}}).

An alternative formula for smoothing is by Cheng (1995) and is
\deqn{\tilde{X}^{\mathrm{Cheng}}_n(F) = \sum_{k=1}^n x_{k:n}\:{n - 1 \choose k-1}\: F^{k-1}(1-F)^{n-k}\mbox{.}}
}
\usage{
dat2bernqua(f, x, bern.control=NULL,
                  poly.type=c("Bernstein", "Kantorovich", "Cheng", "Parzen",
                              "bernstein", "kantorovich", "cheng", "parzen"),
                  bound.type=c("none", "sd", "Carv", "either", "carv"),
                  fix.lower=NULL, fix.upper=NULL, p=0.05, listem=FALSE)
}
\arguments{
  \item{f}{A vector of nonexceedance probabilities \eqn{F}.}
  \item{x}{A vector of data values.}
  \item{bern.control}{A \code{list} that holds \code{poly.type}, \code{bound.type}, \code{fix.lower}, and \code{fix.upper}. And this list will supersede the respective
values provided as separate arguments.}
  \item{poly.type}{The Bernstein or Kantorovich polynomial will be used. The two are quite closely related. Or the formula by Cheng (1995) will be used and \code{bound.type}, \code{fix.lower}, \code{fix.upper}, and \code{p} are not applicable. Or the formula credited by Nair et al. (2013, p. 17) to Parzen (1979) will be used.}
  \item{bound.type}{Triggers to the not involve alternative supports (\code{"none"}) then the minimum and maximum are used unless already provided by the \code{fix.lower} or \code{fix.upper}, the support based \code{"sd"} on the standard deviation, the support \code{"Carv"} based on the arguments of de Carvalho (2012), or \code{"either"} method.}
  \item{fix.lower}{For \eqn{k = 0}, either the known lower bounds is used if provided as non \code{NULL} or the observed minimum of the data. If the minimum of the data is less than the \code{fix.lower}, a warning is triggered and \code{fix.lower} is set to the minimum. Following Turnbull and Ghosh (2014) to avoid bounds that are extremely lower than the data, it will use the estimated lower bounds by the method \code{"sd"}, \code{"Carv"}, or \code{"either"} if these bounds are larger than the provided \code{fix.lower}.}
  \item{fix.upper}{For \eqn{k = n}, either the known upper bounds is used if provided as non \code{NULL} or the observed maximum of the data; If the maximum of the data is less than the \code{fix.upper}, a warning is triggered and \code{fix.upper} is set to the maximum.}
  \item{p}{A small probability value to serve as the \eqn{p} in the \code{"Carv"} support computation. The default is recommended as mentioned above. The program will return \code{NA} if \eqn{10^{-6} < p \ge (1-10^{-6})} is not met. The value \code{p} is the \dQuote{p-factor} \eqn{p}.}
  \item{listem}{A logical controlling whether (1) a vector of \eqn{\tilde{X}_n(F)} is returned or (2) a list containing \eqn{\tilde{X}_n(F)}, the \code{f}, original sample size \eqn{n} of the data, the de Carvalho probability \code{p} (whether actually used internally or not), and both \code{fix.lower} and \code{fix.upper} as computed within the function or provided (less likely) by the function arguments.}
}
\details{
Yet another alternative formula for smoothing if by Parzen (1979) and known as the \dQuote{Parzen weighting method} is
\deqn{\tilde{X}^{\mathrm{Parzen}}_n(F) = n\left(\frac{r}{n} - F\right)x_{r-1:n} + n\left(F - \frac{r-1}{n}\right)x_{r:n}\mbox{,}}
where \eqn{(r-1)/n \le F \le (r/n)} for \eqn{r = 1, 2, \ldots, n} and \eqn{x_{0:n}} is taken as either the minimum of the data (\eqn{\mathrm{min}(x)}) or the lower bounds \code{fix.lower} as externally set by the user.  For protection, the minimum of \eqn{(\mathrm{min}(x),} \code{fix.lower}\eqn{)} is formally used. If the Parzen method is used, the only arguments considered are \code{poly.type} and \code{fix.lower}; all others are ignored including the \code{f} (see Value section). The user does not actually have to provide \code{f} in the arguments but a place holder such as \code{f=NULL} is required; internally the Parzen method takes over full control. The Parzen method in general is not smooth and not recommended like the others that rely on a polynomial basis function. Further the Parzen method has implicit asymmetry in the estimated \eqn{F}. The method has \eqn{F=0} and \eqn{F < 1} on output, but if the data are reversed, then the method has \eqn{F > 0} and \eqn{F=1}. Data reversal is made in \code{-X} as this example illustrates:
\preformatted{
X <- sort(rexp(30))
P <- dat2bernqua(f=NULL,  X, poly.type="Parzen")
R <- dat2bernqua(f=NULL, -X, poly.type="Parzen")
plot(pp(X, a=0.5), X, xlim=c(0, 1)) # Hazen plotting position to
lines(  P$f,  P$x, col="red" )      # basically split the horizontal
lines(1-R$f, -R$x, col="blue")      # differences between blue and red.
}
}
\value{
An \R \code{vector} is returned unless the Parzen weighting method is used and in that case an \R \code{list} is returned with elements \code{f} and \code{x}, which respectively are the \eqn{F} values as shown in the formula and the \eqn{\tilde{X}^{\mathrm{Parzen}}_n(F)}.
}
\note{
\enc{Muñoz-Pérez}{Munoz-Perez} and \enc{Fernández-Palacín}{Fernandez-Palacin} (1987, p. 391) describe what to do with the condition of \eqn{k = 0} but seemingly do not comment on the condition of \eqn{k = n}. There is no 0th-order statistic nor is there a \eqn{k > n} order statistic. \enc{Muñoz-Pérez}{Munoz-Perez} and \enc{Fernández-Palacín}{Fernandez-Palacin} (1987) bring up the notion of a natural minimum for the data (for example, data that must be positive, \code{fix.lower = 0} could be set). Logic dictates that a similar argument must be made for the maximum to keep a critical error from occurring if one tries to access the not plausible \code{x[n+1]}-order statistic. Lastly, the argument names \code{bound.type}, \code{fix.lower}, and \code{fix.upper} mimic, as revisions were made to this function in December 2013, the nomenclature of software for probability density function smoothing by Turnbull and Ghosh (2014). The \code{dat2bernqua} function was originally added to \pkg{lmomco} in May 2013 prior to the author learning about Turnbull and Ghosh (2014).

Lastly, there can be many practical situations in which transformation is desired. Because of the logic structure related to how \code{fix.lower} and \code{fix.upper} are determined or provided by the user, it is highly recommended that this function not internally handle transformation and detransformation. See the second example for use of logarithms.
}
\section{Special Attention}{
The limiting properties of the Bernstein and Kantorovich polynomials differ. The Kantorovich polynomial uses the average of the largest (smallest) value and the respective outer order statistics (\eqn{x_{n+1:n}} or \eqn{x_{0:n}}) unlike the Bernstein polynomial whose \eqn{F = 0} or \eqn{F = 1} are purely a function of the outer order statistics. Thus, the Bernstein polynomial can attain the \code{fix.lower} and(or) \code{fix.upper} whereas the Kantorovich fundamentally can not.  For a final comment, the function \code{\link{dat2bernquaf}} is an inverse of \code{dat2bernqua}.
}
\references{
Cheng, C., 1995, The Bernstein polynomial estimator of a smooth quantile function: Statistics and Probability Letters, v. 24, pp. 321--330.

de Carvalho, M., 2012, A generalization of the Solis-Wets method: Journal of Statistical Planning and Inference, v. 142, no. 3, pp. 633--644.

\enc{Hernández-Maldonado}{Hernandez-Maldonado}, V., \enc{Díaz-Viera}{Diaz-Viera}, M., and Erdely, A., 2012, A joint stochastic simulation method using the Bernstein copula as a flexible tool for modeling nonlinear dependence structures between petrophysical properties: Journal of Petroleum Science and Engineering, v. 90--91, pp. 112--123.

\enc{Muñoz-Pérez}{Munoz-Perez}, J., and \enc{Fernández-Palacín}{Fernandez-Palacin}, A., 1987, Estimating the quantile function by Bernstein polynomials: Computational Statistics and Data Analysis, v. 5, pp. 391--397.

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.

Turnbull, B.C., and Ghosh, S.K., 2014, Unimodal density estimation using Bernstein polynomials: Computational Statistics and Data Analysis, v. 72, pp. 13--29.

Parzen, E., 1979, Nonparametric statistical data modeling: Journal American Statistical Association, v. 75, pp. 105--122.
}
\section{Implentation Note}{
The function makes use of \R functions \code{lchoose} and \code{exp} and logarithmic expressions, such as \eqn{(1-F)^{(n-k)} \rightarrow (n-k)\log(1-F)}, for numerical stability for large sample sizes.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmoms.bernstein}}, \code{\link{pfactor.bernstein}}, \code{\link{dat2bernquaf}} }
\examples{
# Compute smoothed extremes, quartiles, and median
# The smoothing seems to extend to F=0 and F=1.
set.seed(1); X <- exp(rnorm(20)); F <- c(0, .25, .50, .75, 1)
dat2bernqua(F, X, bound.type="none",   listem=TRUE)$x
dat2bernqua(F, X, bound.type="Carv",   listem=TRUE)$x
dat2bernqua(F, X, bound.type="sd",     listem=TRUE)$x
dat2bernqua(F, X, bound.type="either", listem=TRUE)$x
dat2bernqua(F, X, bound.type="sd",     listem=TRUE, fix.lower=0)$x

\dontrun{
X <- sort(10^rnorm(20)); F <- nonexceeds(f01=TRUE)
plot(qnorm(pp(X)), X, xaxt="n", xlab="", ylab="QUANTILE", log="y")
add.lmomco.axis(las=2, tcl=0.5, side.type="NPP", twoside=TRUE)
lines(qnorm(F),     dat2bernqua(F,    X,  bound.type="sd"), col="red", lwd=2)
lines(qnorm(F), exp(dat2bernqua(F,log(X), bound.type="sd"))) # }

\dontrun{
X <- exp(rnorm(20)); F <- seq(0.001, 0.999, by=.001)
dat2bernqua(0.9, X, poly.type="Bernstein",   listem=TRUE)$x
dat2bernqua(0.9, X, poly.type="Kantorovich", listem=TRUE)$x
dat2bernqua(0.9, X, poly.type="Cheng",       listem=TRUE)$x
plot(pp(X), sort(X), log="y", xlim=range(F))
lines(F, dat2bernqua(F, X, poly.type="Bernstein"  ), col="red"  )
lines(F, dat2bernqua(F, X, poly.type="Kantorovich"), col="green")
lines(F, dat2bernqua(F, X, poly.type="Cheng"      ), col="blue" ) #}

\dontrun{
X <- exp(rnorm(20)); F <- nonexceeds()
plot(pp(X), sort(X))
lines(F, dat2bernqua(F,X, bound.type="sd", poly.type="Bernstein"))
lines(F, dat2bernqua(F,X, bound.type="sd", poly.type="Kantorovich"), col=2) #}

\dontrun{
X <- rnorm(25); F <- nonexceeds()
Q <- dat2bernqua(F, X) # the Bernstein estimates
plot( F, dat2bernqua(F, X, bound.type="Carv"), type="l"   )
lines(F, dat2bernqua(F, X, bound.type="sd"),   col="red"  )
lines(F, dat2bernqua(F, X, bound.type="none"), col="green")
points(pp(X),      sort(X), pch=16, cex=.75,   col="blue" ) #}

\dontrun{
set.seed(13)
par <- parkap(vec2lmom(c(1, .5, .4, .2)))
F <- seq(0.001, 0.999, by=0.001)
X <- sort(rlmomco(100, par))
pp <- pp(X)
pdf("lmomco_example_dat2bernqua.pdf")
plot(qnorm(pp(X)), dat2bernqua(pp, X), col="blue", pch=1,
     ylim=c(0,qlmomco(0.9999, par)))
lines(qnorm(F), dat2bernqua(F, sort(X)), col="blue")
lines(qnorm(F),     qlmomco(F,     par), col="red" )
sampar  <- parkap(lmoms(X))
sampar2 <- parkap(lmoms(dat2bernqua(pp, X)))
lines( qnorm(pp(F)), qlmomco(F, sampar ), col="black")
lines( qnorm(pp(F)), qlmomco(F, sampar2), col="blue", lty=2)
points(qnorm(pp(X)), X, col="black", pch=16)
dev.off() #}
}
\keyword{smoothing (Bernstein-Kantorovich)}
\keyword{smoothing (Bernstein polynomial)}
\keyword{smoothing (Kantorovich polynomial)}
\keyword{smoothing (Parzen weighting)}
\keyword{quantile function}

