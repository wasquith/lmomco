\name{pargep}
\alias{pargep}
\title{Estimate the Parameters of the Generalized Exponential Poisson Distribution}
\description{
This function estimates the parameters of the Generalized Exponential Poisson distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomgep}}. However, the expectations of order statistic extrema are computed through numerical integration of the quantile function and the fundamental definition of L-moments (\code{\link{theoLmoms.max.ostat}}). The mean must be \eqn{\lambda_1 > 0}. The implementation here fits the first three L-moments. A distribution having two scale parameters produces more than one solution.  The higher L-moments are not consulted as yet in an effort to further enhance functionality. This function has deterministic starting points but on subsequent iterations the starting points do change. If a solution is not forthcoming, try running the whole function again.
}
\usage{
pargep(lmom, checklmom=TRUE, checkdomain=TRUE, maxit=10, verbose=FALSE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{checkdomain}{A logical controlling whether the empirically derived (approximated) boundaries of the GEP in the \eqn{\tau_2} and \eqn{\tau_3} domain are used for early exiting if the \code{lmom} do not appear compatible with the distribution.}
  \item{maxit}{The maximum number of iterations. The default should be about twice as big as necessary.}
  \item{verbose}{A logical controlling intermediate results, which is useful given the experimental nature of GEP parameter estimation and if the user is evaluating results at each iteration. The verbosity is subject to change.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.
  
  \item{type}{The type of distribution: \code{gep}.}
  \item{para}{The parameters of the distribution.}
  \item{convergence}{A numeric code on covergence, a value of 0 means solution looks ok.}
  \item{error}{Sum of relative error: \eqn{\epsilon = |(\lambda'_2 - \hat\lambda'_2)/\hat\lambda'_2|} \eqn{+} \eqn{|(\lambda_3 - \hat\lambda_3)/\hat\lambda_3|} for the fitted (prime) and sample (hat, given in \code{lmom}) 2nd and 3rd L-moments. A value of 10 means that the \eqn{\tau_2} and \eqn{\tau_3} values are outside the domain of the distribution as determined by brute force computations and custom polynomial fits.}
  \item{its}{Iteration count.}
  \item{source}{The source of the parameters: \dQuote{pargep}.}
}
\note{
There are various inequalities and polynomials demarcating the \eqn{\tau_2} and \eqn{\tau_3} of the distribution. These were developed during a protracted period of investigation into the numerical limits of the distribution with a specific implementation in \pkg{lmomco}. Some of these bounds may or may not be optimal as empirically-arrived estimates of theoretical bounds. The polynomials where carefully assembled however. The straight inequalities are a bit more ad hoc following supervision of domain exploration.  More research is needed but the domain constraint provided should generally produce parameter solutions.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomgep}}, \code{\link{cdfgep}}, \code{\link{pdfgep}}, \code{\link{quagep}}
}
\examples{
\dontrun{
# Two examples well inside the domain but known to produce difficulty in
# the optimization process; pargep() engineered with flexibility to usually
# hit the proper solutions.
mygepA <- pargep(vec2lmom(c(1,0.305,0.270), lscale=FALSE))
mygepB <- pargep(vec2lmom(c(1,0.280,0.320), lscale=FALSE))
}
\dontrun{
gep1 <- vec2par(c(2708, 3, 52), type="gep")
 lmr <- lmomgep(gep1);  print(lmr$lambdas)
gep2 <- pargep(lmr);    print(lmomgep(gep2)$lambdas)
# Note that we are close on matching the L-moments but we do
# not recover the parameters given because to shape parameters.
gep3 <- pargep(lmr, nk=1, nh=2);
x <- quagep(nonexceeds(), gep1)
x <- sort(c(x, quagep(nonexceeds(), gep2)))
plot(x, pdfgep(x, gep1), type="l", lwd=2)
lines(x, pdfgep(x, gep2), lwd=3, col=2)
lines(x, pdfgep(x, gep3), lwd=2, col=3)
}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Generalized Exponential Poisson}
