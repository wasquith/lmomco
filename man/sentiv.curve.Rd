\name{sentiv.curve}
\alias{sentiv.curve}
\title{Compute the Sensitivity Curve for a Single Quantile}
\description{
The \emph{sensitivity curve} (\eqn{SC}) is a means to assess how sensitive a particular statistic \eqn{T_{n+1}} for a sample of size \eqn{n} is to an additional sample \eqn{x} to be included. For the implementation by this function, the statistic \eqn{T} is a specific quantile \eqn{x(F)} of interest set by a nonexceedance probability \eqn{F}.  The \eqn{SC} is
\deqn{SC_{n+1}(x,\,| F) = (n+1)(T_{n+1} - T_n)\mbox{,}}
where \eqn{T_n} represent the statistic for the sample of size \eqn{n}. The notation here follows that of Hampel (1974, p. 384) concerning \eqn{n} and \eqn{n+1}.
}
\usage{
sentiv.curve(f, x, method=c("bootstrap", "polynomial", "none"),
                   data=NULL, para=NULL, ...)
}
\arguments{
  \item{f}{The nonexceedance probability \eqn{F} of the quantile for which the sensitivity of its estimation is needed. Only the first value if a vector is given is used and a warning issued.}
  \item{x}{The \eqn{x} values representing the potential \emph{one more value} to be added to the original data.}
  \item{data}{A vector of mandatory sample data values. These will either be converted to (1) order statistic expectations exact analytical expressions or simulation (backup plan), (2) Bernstein (or similar) polynomials, or (3) the provided values treated as if they are the order statistic expectations.}
  \item{method}{A character variable determining how the statistics \eqn{T} are computed (see Details).}
  \item{para}{A distribution parameter list from a function such as \code{\link{vec2par}} or \code{\link{lmom2par}}.}
  \item{...}{Additional arguments to pass either to the \code{\link{lmoms.bootbarvar}} or to the \cr \code{\link{dat2bernqua}} function.}
}
\details{
The main features of this function involve how the statistics are computed and are controlled by the \code{method} argument. Three different approaches are provided.

\bold{Bootstrap:} Arguments \code{data} and \code{para} are \emph{mandatory}. If \code{boostrap} is requested, then the distribution type set by the \code{type} attribute in \code{para} is used along with the method of L-moments for \eqn{T(F)} estimation. The \eqn{T_n(F)} is directly computed from the distribution in \code{para}. And for each \code{x}, the \eqn{T_{n+1}(F)} is computed by \code{\link{lmoms}}, \code{\link{lmom2par}}, and the distribution type. The sample so fed to \code{\link{lmoms}} is denoted as \code{c(EX, x)}.

\bold{Polynomial:} Argument \code{data} is \emph{mandatory} and \code{para} is \emph{not} used. If \code{polynomial} is requested, then the Bernstein polynomial (likely) from the \code{\link{dat2bernqua}} is used. The \eqn{T_n(F)} is computed by the \code{data} sample. And for each \code{x}, the \eqn{T_{n+1}(F)} also is computed by \code{\link{dat2bernqua}}, but the sample so fed to \code{\link{dat2bernqua}} is denoted as \code{c(EX, x)}.

\bold{None:} Arguments \code{data} and \code{para} are \emph{mandatory}. If \code{none} is requested, then the distribution type set by the \code{type} attribute in \code{para} is used along with the method of L-moments. The \eqn{T_n(F)} is directly computed from the distribution in \code{para}. And for each \code{x}, the \eqn{T_{n+1}(F)} is computed by \code{\link{lmoms}}, \code{\link{lmom2par}}, and the distribution type. The sample so fed to \code{\link{lmoms}} is denoted as \code{c(EX, x)}.

The internal variable \code{EX} now requires discussion. If \code{method=none}, then the \code{data} are sorted and set into the internal variable \code{EX}. Conversely, if \code{method=}\code{bootstrap} or \code{method=}\code{polynomial}, then \code{EX} will contain the expectations of the order statistics from \code{\link{lmoms.bootbarvar}}.

Lastly, the Weibull plotting positions are used for the probability values for the data as provided by the \code{\link{pp}} function. Evidently, if \code{method} is either \code{parent} or \code{polynomial} then a \dQuote{stylized sensitivity curve} would created (David, 1981, p. 165) because the expectations of the sample order statistics and not the sample order statistics (the sorted sample) are used.
}
\value{
An \R \code{list} is returned.

\item{curve}{The value for \eqn{SC(x) = (n+1)(T_{n+1} - T_n)}.}
\item{curve.perchg}{The percent change sensitivity curve by \eqn{SC^{(\%)}(x) = 100\times (T_{n+1} - T_n)/T_n}.}
\item{Tnp1}{The values for \eqn{T_{n+1} = T_n + SC(x)/(n+1)}.}
\item{Tn}{The value (singular) for \eqn{T_n} which was estimated according to \code{method}.}
\item{color}{The curve potentially passes through a zero depending on the values for \eqn{x}. The \code{color} is set to distinquish between negatives and positives so that the user could use the absolute value of \code{curve} on logarithmic scales and use the color to distinquish the original negatives.}
\item{EX}{The values for the internal variable \code{EX}.}
\item{source}{An attribute identifying the computational source of the sensitivity curve: \dQuote{sentiv.curve}.}
}
\references{
David, H.A., 1981, Order statistics: John Wiley, New York.

Hampel, F.R., 1974, The influence curve and its role in robust estimation: Journal of the American Statistical Association, v. 69, no. 346, pp. 383--393.
}
\author{W.H. Asquith}
\seealso{\code{\link{expect.max.ostat}}
}
\examples{
\dontrun{
set.seed(50)
mean <- 12530; lscale <- 5033; lskew <- 0.4
n <- 46; type <- "gev"; lmr <- vec2lmom(c(mean,lscale,lskew))
F <- 0.90 # going to explore sensitivity on the 90th percentile
par.p <- lmom2par(lmr, type=type) # Parent distribution
TRUE.Q <- par2qua(F, par.p)
X <- sort(rlmomco(n, par.p)) # Simulate a small sample
par.s <- lmom2par(lmoms(X), type=type) # Now fit the distribution
SIM.Q <- par2qua(F, par.s); SIM.BAR <- par2lmom(par.s)$lambdas[1]
D <- log10(mean) - log10(lscale)
R <- as.integer(log10(mean)) + c(-D, D) # need some x-values to explore
Xs <- 10^(seq(R[1], R[2], by=.01)) # x-values to explore
# Sample estimate are the "parent" only to mimic a more real-world setting.
# where one "knows" the form of the parent but perhaps not the parameters.
SC1 <- sentiv.curve(F, Xs, data=X, para=par.s, method="bootstrap")
SC2 <- sentiv.curve(F, Xs, data=X, para=par.s, method="polynomial",
                              bound.type="Carv")
SC3 <- sentiv.curve(F, Xs, data=X, para=par.s, method="none")
xlim <- range(c(Xs,SC1$Tnp1,SC2$Tnp1,SC3$Tnp1))
ylim <- range(c(SC1$curve.perchg, SC2$curve.perchg, SC3$curve.perchg))
plot(xlim, c(0,0), type="l", lty=2, ylim=ylim, xaxs="i", yaxs="i",
     xlab=paste("Magnitude of next value added to sample of size",n),
     ylab=paste("Percent change fitted",F,"probability quantile"))
mtext(paste("Distribution",par.s$type,"with parameters",
      paste(round(par.s$para, digits=3), collapse=", ")))
lines(rep(TRUE.Q,  2), c(-10,10), lty=4, lwd=3)
lines(rep(SIM.BAR, 2), c(-10,10), lty=3, lwd=2)
lines(rep(SIM.Q,   2), c(-10,10), lty=2)
lines(Xs, SC1$curve.perchg, lwd=3, col=1)
lines(Xs, SC2$curve.perchg, lwd=2, col=2)
lines(Xs, SC3$curve.perchg, lwd=1, col=4)
rug(SC1$Tnp1, col=rgb(0,0,0,0.3))
rug(SC2$Tnp1, col=rgb(1,0,0,0.3))
rug(SC3$Tnp1, col=rgb(0,0,1,0.3), tcl=-.75) #}
}
\keyword{order statistics (misc.)}
\keyword{Sensitivity curve} 
\keyword{Graphics}
\concept{Sensitivity function}
\concept{sensitivity function}


