\name{lmomTLgld}
\alias{lmomTLgld}
\title{Trimmed L-moments of the Generalized Lambda Distribution}
\description{
This function estimates the symmetrical trimmed L-moments (TL-moments) for \eqn{t=1} of the Generalized Lambda distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h}) from \code{\link{parTLgld}}. The TL-moments in terms of the parameters are complicated; however, there are analytical solutions. There are no simple expressions of the parameters in terms of the L-moments. The first four TL-moments (trim = 1) of the distribution are
\deqn{\lambda^{(1)}_1 = \xi + 6\alpha
                              \left(\frac{1}{(\kappa+3)(\kappa+2)} -
               \frac{1}{(h+3)(h+2)} \right) \mbox{,}}
\deqn{\lambda^{(1)}_2 = 6\alpha \left(\frac{\kappa}{(\kappa+4)(\kappa+3)(\kappa+2)} +
                                           \frac{h}{(h+4)(h+3)(h+2)}\right) \mbox{,}}
\deqn{\lambda^{(1)}_3 =  \frac{20\alpha}{3} \left(\frac{\kappa (\kappa - 1)}
                                                {(\kappa+5)(\kappa+4)(\kappa+3)(\kappa+2)} -
        \frac{h (h - 1)}
             {(h+5)(h+4)(h+3)(h+2)} \right) \mbox{,}}
\deqn{\lambda^{(1)}_4 = \frac{15\alpha}{2} \left(\frac{\kappa (\kappa - 2)(\kappa - 1)}
                                                {(\kappa+6)(\kappa+5)(\kappa+4)(\kappa+3)(\kappa+2)} +
        \frac{h (h - 2)(h - 1)}
             {(h+6)(h+5)(h+4)(h+3)(h+2)} \right) \mbox{,}}
\deqn{\lambda^{(1)}_5 = \frac{42\alpha}{5} \left(N1 - N2 \right) \mbox{,}}
where
\deqn{N1 = \frac{\kappa (\kappa - 3)(\kappa - 2)(\kappa - 1) }
                {(\kappa+7)(\kappa+6)(\kappa+5)(\kappa+4)(\kappa+3)(\kappa+2)} \mbox{ and}}
\deqn{N2 = \frac{h (h - 3)(h - 2)(h - 1)}{(h+7)(h+6)(h+5)(h+4)(h+3)(h+2)} \mbox{.}}

The TL-moment (\eqn{t=1}) for \eqn{\tau^{(1)}_3} is
\deqn{\tau^{(1)}_3 = \frac{10}{9} \left( \frac{\kappa(\kappa-1)(h+5)(h+4)(h+3)(h+2) -
                     h(h-1)(\kappa+5)(\kappa+4)(\kappa+3)(\kappa+2)}
      {(\kappa+5)(h+5) \times [\kappa(h+4)(h+3)(h+2) +
                                  h(\kappa+4)(\kappa+3)(\kappa+2)]
      } \right)
      \mbox{.}}

The TL-moment (\eqn{t=1}) for \eqn{\tau^{(1)}_4} is
\deqn{N1 = \kappa(\kappa-2)(\kappa-1)(h+6)(h+5)(h+4)(h+3)(h+2) \mbox{,}}
\deqn{N2 = h(h-2)(h-1)(\kappa+6)(\kappa+5)(\kappa+4)(\kappa+3)(\kappa+2) \mbox{,}}
\deqn{D1 = (\kappa+6)(h+6)(\kappa+5)(h+5) \mbox{,}}
\deqn{D2 = [\kappa(h+4)(h+3)(h+2) + h(\kappa+4)(\kappa+3)(\kappa+2)] \mbox{, and}}
\deqn{\tau^{(1)}_4 = \frac{5}{4} \left( \frac{N1 +  N2}{D1 \times D2} \right) \mbox{.}}

Finally the TL-moment (\eqn{t=1}) for \eqn{\tau^{(1)}_5} is
\deqn{N1 = \kappa(\kappa-3)(\kappa-2)(\kappa-1)(h+7)(h+6)(h+5)(h+4)(h+3)(h+2) \mbox{,}}
\deqn{N2 = h(h-3)(h-2)(h-1)(\kappa+7)(\kappa+6)(\kappa+5)(\kappa+4)(\kappa+3)(\kappa+2) \mbox{,}}
\deqn{D1 = (\kappa+7)(h+7)(\kappa+6)(h+6)(\kappa+5)(h+5) \mbox{,}}
\deqn{D2 = [\kappa(h+4)(h+3)(h+2) + h(\kappa+4)(\kappa+3)(\kappa+2)] \mbox{, and}}
\deqn{\tau^{(1)}_5 = \frac{7}{5} \left( \frac{N1 -  N2}{D1 \times D2} \right)\mbox{.}}

By inspection the \eqn{\tau_r} equations are not applicable for negative integer values \eqn{k=\{-2, -3, -4, \dots \}} and \eqn{h=\{-2, -3, -4, \dots \}} as division by zero will result. There are additional, but difficult to formulate, restrictions on the parameters both to define a valid Generalized Lambda distribution as well as valid L-moments. Verification of the parameters is conducted through \code{\link{are.pargld.valid}}, and verification of the L-moment validity is conducted through \code{\link{are.lmom.valid}}.
}
\usage{
lmomTLgld(para, nmom=6, trim=1, leftrim=NULL, rightrim=NULL, tau34=FALSE)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{nmom}{Number of L-moments to compute.}
  \item{trim}{Symmetrical trimming level set to unity as the default.}
  \item{leftrim}{Left trimming level, \eqn{t_1}.}
  \item{rightrim}{Right trimming level, \eqn{t_2}.}
  \item{tau34}{A logical controlling the level of L-moments returned by the function. If true, then this function returns only \eqn{\tau_3} and \eqn{\tau_4}; this feature might be useful in certain research applications of the Generalized Lambda distribution associated with the multiple solutions possible for the distribution.}
}
\details{
The opening comments in the description pertain to single and symmetrical endpoint trimming, which has been extensively considered by Asquith (2007).  Deriviations backed by numerical proofing of variable arrangement in March 2011 led the the inclusion of the following generalization of the L-moments and TL-moments of the Generalized Lambda shown in Asquith (2011) that was squeezed in late ahead of the deadlines for that monograph.
\deqn{
\lambda^{(t_1,t_2)}_{r} = \alpha (r^{-1}) (r+t_1+t_2) \sum_{j=0}^{r-1} (-1)^{r}{r-1 \choose j}{r+t_1+t_2-1 \choose r+t_1-j-1} \times A\mbox{,}}
where \eqn{A} is
\deqn{A = \biggl(\frac{\Gamma(\kappa+r+t_1-j)\Gamma(t_2+j+1)}{\Gamma(\kappa+r+t_1+t_2+1)} - \frac{\Gamma(r+t_1-j)\Gamma(h+t_2+j+1)}{\Gamma(h+r+t_1+t_2+1)}\biggr)\mbox{,}}
where for the special condition of \eqn{r = 1}, the real mean is
\deqn{
\mathrm{mean} = \xi + \lambda^{(t_1,t_2)}_{1}\mbox{,}}
but for \eqn{r \ge 2} the \eqn{\lambda^{(t_1,t_2)}} provides correct values.  So care is needed algorithmically also when \eqn{\tau^{(t1, t2)}_2} is computed. Inspection of the \eqn{\Gamma(\cdot)} arguments, which must be \eqn{ > 0}, shows that
\deqn{
\kappa > -(1+t_1)}
and
\deqn{h > -(1+t_2)
\mbox{.}}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the TL-moments. First element is
\eqn{\lambda^{(t_1,t_2)}_1}, second element is \eqn{\lambda^{(t_1,t_2)}_2}, and so on.}
  \item{ratios}{Vector of the TL-moment ratios. Second element is
\eqn{\tau^{(1)}}, third element is \eqn{\tau^{(1)}_3} and so on. }
  \item{trim}{Trim level = left or right values if they are equal. The default for this function is \code{trim = 1} because the \code{\link{lmomgld}} provides for \code{trim = 0}. }
  \item{leftrim}{Left trimming level}
  \item{rightrim}{Right trimming level}
  \item{source}{An attribute identifying the computational source of the TL-moments: \dQuote{lmomTLgld}.}
}
\source{
Derivations conducted by W.H. Asquith on February 18 and 19, 2006 and others in early March 2011.
}
\references{
Asquith, W.H., 2007, L-moments and TL-moments of the generalized lambda distribution: Computational Statistics and Data Analysis, v. 51, no. 9, pp. 4484--4496.

Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational statistics and data analysis, v. 43, pp. 299--314.

Karian, Z.A., and Dudewicz, E.J., 2000, Fitting statistical distributions---The generalized lambda distribution and generalized bootstrap methods: CRC Press, Boca Raton, FL, 438 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomgld}}, \code{\link{parTLgld}}, \code{\link{pargld}}, \code{\link{cdfgld}}, \code{\link{quagld}}}
\examples{

\dontrun{
lmomgld(vec2par(c(10,10,0.4,1.3), type='gld'))

PARgld <- vec2par(c(15,12,1,.5), type="gld")
theoTLmoms(PARgld, leftrim=0, rightrim=0, nmom=6)
lmomTLgld(PARgld, leftrim=0, rightrim=0)

theoTLmoms(PARgld, trim=2, nmom=6)
lmomTLgld(PARgld, trim=2)

theoTLmoms(PARgld, trim=3, nmom=6)
lmomTLgld(PARgld, leftrim=3, rightrim=3)

theoTLmoms(PARgld, leftrim=10, rightrim=2, nmom=6)
lmomTLgld(PARgld, leftrim=10, rightrim=2)
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Generalized Lambda}
