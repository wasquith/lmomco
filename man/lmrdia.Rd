\name{lmrdia}
\alias{lmrdia}
\title{L-moment Ratio Diagram Components}
\description{
This function returns a list of the L-skew and L-kurtosis (\eqn{\tau_3} and \eqn{\tau_4}, respectively) ordinates for construction of L-moment Ratio (L-moment diagrams) that are useful in selecting a distribution to model the data.
}
\usage{
lmrdia()
}
\value{
  An \R \code{list} is returned.

  \item{limits}{The theoretical limits of \eqn{\tau_3} and \eqn{\tau_4}; below \eqn{\tau_4} of the theoretical limits are theoretically
  not possible.}
  \item{aep4}{\eqn{\tau_3} and \eqn{\tau_4} lower limits of the Asymmetric Exponential Power distribution.}
  \item{cau}{\eqn{\tau^{(1)}_3 = 0} and \eqn{\tau^{(1)}_4 = 0.34280842} of the Cauchy distribution (TL-moment [trim=1])  (see \bold{Examples} \code{\link{lmomcau}} for source).}
  \item{exp}{\eqn{\tau_3} and \eqn{\tau_4} of the Exponential distribution.}
  \item{gev}{\eqn{\tau_3} and \eqn{\tau_4} of the Generalized Extreme Value distribution.}
  \item{glo}{\eqn{\tau_3} and \eqn{\tau_4} of the Generalized Logistic distribution.}
  \item{gpa}{\eqn{\tau_3} and \eqn{\tau_4} of the Generalized Pareto distribution.}
  \item{gum}{\eqn{\tau_3} and \eqn{\tau_4} of the Gumbel distribution.}
  \item{gno}{\eqn{\tau_3} and \eqn{\tau_4} of the Generalized Normal distribution.}
  \item{gov}{\eqn{\tau_3} and \eqn{\tau_4} of the Govindarajulu distribution.}
  \item{ray}{\eqn{\tau_3} and \eqn{\tau_4} of the Rayleigh distribution.}
  \item{lognormal}{\eqn{\tau_3} and \eqn{\tau_4} of the Generalized Normal (3-parameter Log-Normal) distribution.}
  \item{nor}{\eqn{\tau_3} and \eqn{\tau_4} of the Normal distribution.}
  \item{pe3}{\eqn{\tau_3} and \eqn{\tau_4} of the Pearson Type III distribution.}
  \item{pdq3}{\eqn{\tau_3} and \eqn{\tau_4} of the Polynomial Density-Quantile3 distribution.}
  \item{rgov}{\eqn{\tau_3} and \eqn{\tau_4} of the reversed Govindarajulu.}
  \item{rgpa}{\eqn{\tau_3} and \eqn{\tau_4} of the reversed Generalized Pareto.}
  \item{sla}{\eqn{\tau^{(1)}_3 = 0} and \eqn{\tau^{(1)}_4 = 0.30420472} of the Slash distribution (TL-moment [trim=1]) (see \bold{Examples} \code{\link{lmomsla}} for source).}
  \item{uniform}{\eqn{\tau_3} and \eqn{\tau_4} of the uniform distribution.}
  \item{wei}{\eqn{\tau_3} and \eqn{\tau_4} of the Weibull distribution (reversed Generalized Extreme Value).}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.

Hosking, J.R.M., 1986, The theory of probability weighted moments: Research Report RC12210, IBM Research Division, Yorkton Heights, N.Y.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2,870--2,891, \doi{10.1016/j.jspi.2006.10.010}.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{plotlmrdia}}}
\examples{
lratios <- lmrdia()
}
\keyword{L-moment ratio diagram}
\keyword{L-moment ratio diagram (Tau3 and Tau4)}
