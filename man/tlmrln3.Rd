\name{tlmrln3}
\alias{tlmrln3}
\title{ Compute Select TL-moment ratios of the 3-Parameter Log-Normal Distribution }
\description{
This function computes select TL-moment ratios of the Log-Normal3 distribution for defaults of \eqn{\zeta = 0} and \eqn{\mu_\mathrm{log} = 0}. This function can be useful for plotting the trajectory of the distribution on TL-moment ratio diagrams of \eqn{\tau^{(t_1,t_2)}_2}, \eqn{\tau^{(t_1,t_2)}_3}, \eqn{\tau^{(t_1,t_2)}_4}, \eqn{\tau^{(t_1,t_2)}_5}, and \eqn{\tau^{(t_1,t_2)}_6}. In reality, \eqn{\tau^{(t_1,t_2)}_2} is dependent on the values for \eqn{\zeta} and \eqn{\mu_\mathrm{log}}. If the message \preformatted{Error in integrate(XofF, 0, 1) : the integral is probably divergent} occurs then careful adjustment of the shape parameter \eqn{\sigma_\mathrm{log}} parameter range is very likely required. Remember that TL-moments with nonzero trimming permit computation of TL-moments into parameter ranges beyond those recognized for the usual (untrimmed) L-moments.
}
\usage{
tlmrln3(trim=NULL, leftrim=NULL, rightrim=NULL,
        zeta=0, mulog=0, sbeg=0.01, send=3.5, by=.1)
}
\arguments{
  \item{trim}{Level of symmetrical trimming to use in the computations.
Although \code{NULL} in the argument list, the default is 0---the usual L-moment ratios are returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample.}
  \item{rightrim}{Level of trimming of the right-tail of the sample.}
  \item{zeta}{Location parameter of the distribution.}
  \item{mulog}{Mean of the logarithms of the distribution.}
  \item{sbeg}{The beginning \eqn{\sigma_\mathrm{log}} value of the distribution.}
  \item{send}{The ending \eqn{\sigma_\mathrm{log}} value of the distribution.}
  \item{by}{The increment for the \code{seq()} between \code{sbeg} and \code{send}.}
}
\value{
  An \R \code{list} is returned.

  \item{tau2}{A vector of the \eqn{\tau^{(t_1,t_2)}_2} values.}
  \item{tau3}{A vector of the \eqn{\tau^{(t_1,t_2)}_3} values.}
  \item{tau4}{A vector of the \eqn{\tau^{(t_1,t_2)}_4} values.}
  \item{tau5}{A vector of the \eqn{\tau^{(t_1,t_2)}_5} values.}
  \item{tau6}{A vector of the \eqn{\tau^{(t_1,t_2)}_6} values.}
}
\note{
The function uses numerical integration of the quantile function of the distribution through the \code{\link{theoTLmoms}} function.
}
\author{ W.H. Asquith }
\seealso{\code{\link{qualn3}}, \code{\link{theoTLmoms}}, \code{\link{tlmrgno}}
}
\examples{
\dontrun{
  # Recalling that generalized Normal and log-Normal3 are
  # the same with the GNO being the more general.

  # Plot and L-moment ratio diagram of Tau3 and Tau4
  # with exclusive focus on the GNO distribution.
  plotlmrdia(lmrdia(), autolegend=TRUE, xleg=-.1, yleg=.6,
             xlim=c(-.8, .7), ylim=c(-.1, .8),
             nolimits=TRUE, noglo=TRUE, nogpa=TRUE, nope3=TRUE,
             nogev=TRUE, nocau=TRUE, noexp=TRUE, nonor=TRUE,
             nogum=TRUE, noray=TRUE, nouni=TRUE)

  LN3 <- tlmrln3(sbeg=.001, mulog=-1)
  lines(LN3$tau3, LN3$tau4) # See how it overplots the GNO
  # for right skewness. So only part of the GNO is covered.

  # Compute the TL-moment ratios for trimming of one
  # value on the left and four on the right.
  J <- tlmrgno(kbeg=-3.5, kend=3.9, leftrim=1, rightrim=4)
  lines(J$tau3, J$tau4, lwd=2, col=2) # RED CURVE

  LN3 <- tlmrln3(, leftrim=1, rightrim=4, sbeg=.001)
  lines(LN3$tau3, LN3$tau4) # See how it again over plots
  # only part of the GNO
}
}
\keyword{L-moment (trimmed distribution)}
\keyword{Distribution: Log-Normal (3-parameter)}
