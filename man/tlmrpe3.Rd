\name{tlmrpe3}
\alias{tlmrpe3}
\title{ Compute Select TL-moment ratios of the Pearson Type III }
\description{
This function computes select TL-moment ratios of the Pearson Type III distribution for defaults of \eqn{\xi = 0} and \eqn{\beta = 1}. This function can be useful for plotting the trajectory of the distribution on TL-moment ratio diagrams of \eqn{\tau^{(t_1,t_2)}_2}, \eqn{\tau^{(t_1,t_2)}_3}, \eqn{\tau^{(t_1,t_2)}_4}, \eqn{\tau^{(t_1,t_2)}_5}, and \eqn{\tau^{(t_1,t_2)}_6}. In reality, \eqn{\tau^{(t_1,t_2)}_2} is  dependent on the values for \eqn{\xi} and \eqn{\alpha}. If the message \preformatted{Error in integrate(XofF, 0, 1) : the integral is probably divergent} occurs then careful adjustment of the shape parameter \eqn{\beta} parameter range is very likely required. Remember that TL-moments with nonzero trimming permit computation of TL-moments into parameter ranges beyond those recognized for the usual (untrimmed) L-moments. The function uses numerical integration of the quantile function of the distribution through the \code{\link{theoTLmoms}} function.
}
\usage{
tlmrpe3(trim=NULL, leftrim=NULL, rightrim=NULL,
        xi=0, beta=1, abeg=-.99, aend=0.99, by=.1)
}
\arguments{
  \item{trim}{Level of symmetrical trimming to use in the computations.
Although \code{NULL} in the argument list, the default is 0---the usual L-moment ratios are returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample.}
  \item{rightrim}{Level of trimming of the right-tail of the sample.}
  \item{xi}{Location parameter of the distribution.}
  \item{beta}{Scale parameter of the distribution.}
  \item{abeg}{The beginning \eqn{\alpha} value of the distribution.}
  \item{aend}{The ending \eqn{\alpha} value of the distribution.}
  \item{by}{The increment for the \code{seq()} between \code{abeg} and \code{aend}.}
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
\seealso{\code{\link{quape3}}, \code{\link{theoTLmoms}}}
\examples{
\dontrun{
tlmrpe3(leftrim=2, rightrim=4, xi=0, beta=2)
tlmrpe3(leftrim=2, rightrim=4, xi=100, beta=20) # another slow example
  # Plot and L-moment ratio diagram of Tau3 and Tau4
  # with exclusive focus on the PE3 distribution.
  plotlmrdia(lmrdia(), autolegend=TRUE, xleg=-.1, yleg=.6,
             xlim=c(-.8, .7), ylim=c(-.1, .8),
             nolimits=TRUE, nogev=TRUE, nogpa=TRUE, noglo=TRUE,
             nogno=TRUE, nocau=TRUE, noexp=TRUE, nonor=TRUE,
             nogum=TRUE, noray=TRUE, nouni=TRUE)

  # Compute the TL-moment ratios for trimming of one
  # value on the left and four on the right. Notice the
  # expansion of the alpha parameter space from
  # -1 < a < -1 to something larger based on manual
  # adjustments until blue curve encompassed the plot.
  J <- tlmrpe3(abeg=-15, aend=6, leftrim=1, rightrim=4)
  lines(J$tau3, J$tau4, lwd=2, col=2) # RED CURVE

  # Compute the TL-moment ratios for trimming of four
  # values on the left and one on the right.
  J <- tlmrpe3(abeg=-6, aend=10, leftrim=4, rightrim=1)
  lines(J$tau3, J$tau4, lwd=2, col=4) # BLUE CURVE

  # The abeg and aend can be manually changed to see how
  # the resultant curve expands or contracts on the
  # extent of the L-moment ratio diagram.
}
\dontrun{
  # Following up, let us plot the two quantile functions
  LM  <- vec2par(c(0,1,0.99), type='pe3', paracheck=FALSE)
  TLM <- vec2par(c(0,1,3.00), type='pe3', paracheck=FALSE)
  F <- nonexceeds()
  plot(qnorm(F),  quape3(F, LM), type="l")
  lines(qnorm(F), quape3(F, TLM, paracheck=FALSE), col=2)
  # Notice how the TLM parameterization runs off towards
  # infinity much much earlier than the conventional
  # near limits of the PE3.
}
}
\keyword{L-moment (trimmed, distribution)}
\keyword{Distribution: Pearson Type III}
