\name{lmomgld}
\alias{lmomgld}
\title{L-moments of the Generalized Lambda Distribution}
\description{
This function estimates the L-moments of the Generalized Lambda distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h})  from \code{\link{vec2par}}. The L-moments in terms of the parameters are complicated; however, there are analytical solutions. There are no simple expressions of the parameters in terms of the L-moments. The first L-moment or the mean is
\deqn{\lambda_1 = \xi + \alpha
                              \left(\frac{1}{\kappa+1} -
			            \frac{1}{h+1} \right) \mbox{.}}

The second L-moment or L-scale in terms of the parameters and the mean is
\deqn{\lambda_2 = \xi + \frac{2\alpha}{(\kappa+2)} -
                      2\alpha
                      \left( \frac{1}{h+1} -
		             \frac{1}{h+2} \right) - \xi \mbox{.}}

The third L-moment in terms of the parameters, the mean, and L-scale is
\deqn{Y = 2\xi + \frac{6\alpha}{(\kappa+3)} -
                                         3(\alpha+\xi) + \xi \mbox{, and}}
\deqn{\lambda_3 = Y + 6\alpha
                      \left(\frac{2}{h+2} -
		            \frac{1}{h+3} -
			    \frac{1}{h+1}\right) \mbox{.}}

The fourth L-moment in termes of the parameters and the first three L-moments is
\deqn{Y = \frac{-3}{h+4}\left(\frac{2}{h+2} -
                                \frac{1}{h+3} -
				\frac{1}{h+1}\right) \mbox{,}}
\deqn{Z = \frac{20\xi}{4} + \frac{20\alpha}{(\kappa+4)} -
                                              20 Y\alpha \mbox{, and}}
\deqn{\lambda_4 = Z -
                         5(\kappa + 3(\alpha+\xi) - \xi) +
                                6(\alpha + \xi) - \xi \mbox{.}}

It is conventional to express L-moments in terms of only the parameters and not the other L-moments. Lengthy algebra and further manipulation yields such a system of equations. The L-moments are
\deqn{\lambda_1 = \xi + \alpha
                              \left(\frac{1}{\kappa+1} -
			            \frac{1}{h+1} \right) \mbox{,}}
\deqn{\lambda_2 = \alpha \left(\frac{\kappa}{(\kappa+2)(\kappa+1)} +
                                           \frac{h}{(h+2)(h+1)}\right) \mbox{,}}
\deqn{\lambda_3 =  \alpha \left(\frac{\kappa (\kappa - 1)}
                                                {(\kappa+3)(\kappa+2)(\kappa+1)} -
					   \frac{h (h - 1)}
					        {(h+3)(h+2)(h+1)} \right) \mbox{, and}}
\deqn{\lambda_4 = \alpha \left(\frac{\kappa (\kappa - 2)(\kappa - 1)}
                                                {(\kappa+4)(\kappa+3)(\kappa+2)(\kappa+1)} +
					   \frac{h (h - 2)(h - 1)}
					        {(h+4)(h+3)(h+2)(h+1)} \right) \mbox{.}}

The L-moment ratios are
\deqn{\tau_3 = \frac{\kappa(\kappa-1)(h+3)(h+2)(h+1) -
                     h(h-1)(\kappa+3)(\kappa+2)(\kappa+1)}
		    {(\kappa+3)(h+3) \times [\kappa(h+2)(h+1) +
		                                h(\kappa+2)(\kappa+1)]
		    }
		    \mbox{, and}}
\deqn{\tau_4 = \frac{\kappa(\kappa-2)(\kappa-1)(h+4)(h+3)(h+2)(h+1) +
                     h(h-2)(h-1)(\kappa+4)(\kappa+3)(\kappa+2)(\kappa+1)}
		    {(\kappa+4)(h+4)(\kappa+3)(h+3) \times [\kappa(h+2)(h+1) +
		                                h(\kappa+2)(\kappa+1)]
		    }
		    \mbox{.}}

The pattern being established through symmetry, even higher L-moment ratios are readily obtained. Note the alternating substraction and addition of the two terms in the numerator of the L-moment ratios (\eqn{\tau_r}). For odd \eqn{r \ge 3} substraction is seen and for even \eqn{r \ge 3} addition is seen. For example, the fifth L-moment ratio is
\deqn{N1 = \kappa(\kappa-3)(\kappa-2)(\kappa-1)(h+5)(h+4)(h+3)(h+2)(h+1) \mbox{,}}
\deqn{N2 = h(h-3)(h-2)(h-1)(\kappa+5)(\kappa+4)(\kappa+3)(\kappa+2)(\kappa+1) \mbox{,}}
\deqn{D1 = (\kappa+5)(h+5)(\kappa+4)(h+4)(\kappa+3)(h+3) \mbox{,}}
\deqn{D2 = [\kappa(h+2)(h+1) + h(\kappa+2)(\kappa+1)] \mbox{, and}}
\deqn{\tau_5 = \frac{N1 - N2}{D1 \times D2} \mbox{.}}

By inspection the \eqn{\tau_r} equations are not applicable for negative integer values \eqn{k=\{-1, -2, -3, -4, \dots \}} and \eqn{h=\{-1, -2, -3, -4, \dots \}} as division by zero will result. There are additional, but difficult to formulate, restrictions on the parameters both to define a valid Generalized Lambda distribution as well as valid L-moments. Verification of the parameters is conducted through \code{\link{are.pargld.valid}}, and verification of the L-moment validity is conducted through \code{\link{are.lmom.valid}}.
}
\usage{
lmomgld(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R list is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational source
of the L-moments: \dQuote{lmomgld}.}
}
\source{
Derivations conducted by W.H. Asquith on February 11 and 12, 2006.
}
\references{
Asquith, W.H., 2007, L-moments and TL-moments of the generalized lambda distribution: Computational Statistics and Data Analysis, v. 51, no. 9, pp. 4484--4496.

Karvanen, J., Eriksson, J., and Koivunen, V., 2002, Adaptive score functions for maximum likelihood ICA: Journal of VLSI Signal Processing, v. 32, pp. 82--92.

Karian, Z.A., and Dudewicz, E.J., 2000, Fitting statistical distibutions---The generalized lambda distribution and generalized bootstrap methods: CRC Press, Boca Raton, FL, 438 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pargld}}, \code{\link{cdfgld}}, \code{\link{pdfgld}},  \code{\link{quagld}}}
\examples{
\dontrun{
lmomgld(vec2par(c(10,10,0.4,1.3),type='gld'))
}

\dontrun{
PARgld <- vec2par(c(0,1,1,.5), type="gld")
theoTLmoms(PARgld, nmom=6)
lmomgld(PARgld)
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Generalized Lambda}
