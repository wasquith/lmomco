\name{lmrdia46}
\alias{lmrdia46}
\title{L-moment Ratio Diagram Components of Tau4 and Tau6}
\description{
This function returns a list of the L-kurtosis (\eqn{\tau_4} and sixth L-moment ratio \eqn{\tau_6}, respectively) ordinates for construction of L-moment Ratio (L-moment diagrams) that are useful in selecting a distribution to model the data.
}
\usage{
lmrdia46()
}
\value{
  An \R \code{list} is returned.

  \item{aep4}{\eqn{\tau_4} and \eqn{\tau_6} of the 4-parameter Asymmetric Exponential Power (AEP4) distribution given L-skew set as \eqn{\tau_3 = 0}. This becomes then the (Symmetrical) Exponential Power. The complementary entry \code{pwrexp} are the effectively the same curve for the power exponential distribution based on lookup table archived in the \pkg{lmomco} package. The table stems from \code{inst/doc/SysDataBuilder02.R}. The \code{aep4} and not \code{pwrexp} is used in the line drawing by \code{\link{plotlmrdia46}}.}
  \item{gld}{\eqn{\tau_4} and \eqn{\tau_6} of the Generalized Lambda distribution given L-skew set as \eqn{\tau_3 = 0}. The table stems from \code{inst/doc/SysDataBuilder02.R}. The \code{gld} is used in the line drawing by \code{\link{plotlmrdia46}}.}
  \item{nor}{\eqn{\tau_4} and \eqn{\tau_6} of the Normal distribution. The table stems from \code{inst/doc/}\cr\code{SysDataBuilder02.R}. The \code{nor} is used in the point drawing by \cr\code{\link{plotlmrdia46}}.}
  \item{pdq4}{\eqn{\tau_4} and \eqn{\tau_6} of the Polynomial Density-Quantile4 distribution, which implicitly is symmetrical, and therefore L-skew set as \eqn{\tau_3 = 0}. The table stems from \code{inst/doc/SysDataBuilder02.R}. The \code{pdq4} is used in the line drawing by\cr\code{\link{plotlmrdia46}}.}
  \item{pwrexp}{\eqn{\tau_4} and \eqn{\tau_6} of the Power Exponential distribution of which the Asymmetric Exponential Power distribution (see also \code{\link{lmomaep4}}). The lookup table archive in the \pkg{lmomco} package for the Power Exponential (\code{PowerExponential.txt}) is confirmed to match the computation in \code{aep4} based on the AEP4 instead. The table stems from \code{inst/doc/}\cr\code{SysDataBuilder02.R}.}
  \item{st2}{\eqn{\tau_4} and \eqn{\tau_6} of the well-known Student t distribution. The lookup table archive in the \pkg{lmomco} package for the Student t (\code{StudentT.txt}) is confirmed to match the computation in \code{st3} based on the ST3 instead. The table stems from \cr\code{inst/doc/SysDataBuilder02.R}. The \code{st3} and not \code{st2} is used in the line drawing by \code{\link{plotlmrdia46}}.}
  \item{st3}{\eqn{\tau_4} and \eqn{\tau_6} of the Student 3t distribution (\code{\link{lmomst3}}). The table stems from \cr\code{inst/doc/SysDataBuilder02.R}. The \code{st3} and not \code{st2} is used in the line drawing by \code{\link{plotlmrdia46}}.}
  \item{symstable}{\eqn{\tau_4} and \eqn{\tau_6} of the Stable distribution, which is not otherwise supported in \pkg{lmomco}. The lookup table archive in the \pkg{lmomco} package for the Symmetrical Stable distribution is \code{StableDistribution.txt}. The table stems from \cr\code{inst/doc/SysDataBuilder02.R}. The \code{symstable} is used in the line drawing by \code{\link{plotlmrdia46}}.}
  \item{tukeylam}{\eqn{\tau_4} and \eqn{\tau_6} of the Tukey Lambda distribution (\url{https://en.wikipedia.org/wiki/Tukey_lambda_distribution}) that is not supported in \pkg{lmomco} but the Generalized Lambda distribution is instead. The \code{SymTukeyLambda.txt} is the lookup table archive in the \pkg{lmomco} package for the Tukey Lambda distribution confirmed to match the mathematics shown herein. The measure \eqn{L-scale} or the second L-moment is not constant for the Symmetric Tukey Lambda as formulated. So, the trajectory of this distribution is not for a constant L-scale, which is unlike that for the Generalized Lambda. The table stems from \code{inst/doc/SysDataBuilder02.R}. The \code{tukeylam} is used in the line drawing by \code{\link{plotlmrdia46}}.}
}
\details{
The \code{\link{lmrdia46}} returns a list of the tables for drawing the trajectories of the distributions by its access of \code{.lmomcohash$t46list} created by the \code{inst/doc/SysDataBuilder02.R} script for \code{sysdata.rda} construction used by the \pkg{lmomco} package itself. A lookup table for the Exponential Power distribution is provided as \code{inst/doc/t4t6/PowerExponential.txt} (\code{.lmomcohash$}\code{tau46list$}\code{pwrexp}), and this distribution is a special case of the Asymmetric Exponential Power4 (\code{\link{lmomaep4}})\cr(\code{.lmomcohash$}\code{tau46list$}\code{aep4}). A lookup table for the Symmetric Stable distribution is provided as \code{inst/doc/t4t6/StableDistribution.txt} (\code{.lmomcohash$}\code{tau46list$}\code{symstable}). A lookup table for the Student t distribution is provided as \code{inst/doc/t4t6/StudentT.txt}\cr(\code{.lmomcohash$}\code{tau46list$}\code{st2}), and this distribution is the same as the Student 3t (\code{\link{lmomst3}}) (\code{.lmomcohash$}\code{tau46list$}\code{st3}). A lookup table for the Tukey Lamda  distribution is provided as \code{inst/doc/t4t6/SymTukeyLambda.txt} (\code{.lmomcohash$}\code{tau46list$}\code{tukeylam}), and this distribution is not quite the same as the Generalized Lambda distribution (\code{\link{lmomgld}}) (\code{.lmomcohash$}\cr\code{tau46list$}\code{gld}). The normal distribution plots as a point in a Tau4-Tau6 L-moment ratio diagram as \code{.lmomcohash$}\code{tau46list$}\code{nor} for which \eqn{\tau_4^\mathrm{nor} = 30/\pi\times \mathrm{atan}(\sqrt{2}) - 9} \eqn{= 0.1226017} and \eqn{\tau_6^\mathrm{nor} = 0.04365901} (numerical integration).

The Tukey Lambda distribution is a simpler formulation than the Generalized Lambda.

\deqn{Q(F) = \frac{1}{\lambda} \biggl[F^\lambda - (1-F)^\lambda \biggr]}

for nonexceedance probability \eqn{F} and \eqn{\lambda \ne 0} and

\deqn{Q(F) = \mathrm{log}\biggl(\frac{F}{1-F}\biggr)}

for \eqn{\lambda = 0} using the natural logarithm.

Inspection of the distribution formulae inform us that the variation in the distribution, the scaling factor \eqn{1/\lambda} to far left in the first definition, for instance, implies that the L-scale (\eqn{\lambda_2}) is not constant and varies with \eqn{\lambda}. The second L-moment of the Tukey Lambda (all odd order L-moments are zero) is

\deqn{\lambda_2 = \frac{2}{\lambda}\biggl[ -\frac{1}{1+\lambda} + \frac{2}{2+\lambda}\biggr]}

and the fourth and sixth L-moments are

\deqn{\lambda_4 = \frac{2}{\lambda}\biggl[ -\frac{1}{1+\lambda} + \frac{12}{2+\lambda} - \frac{30}{3+\lambda} + \frac{20}{4+\lambda}\biggr]}

\deqn{\lambda_6 = \frac{2}{\lambda}\biggl[ -\frac{1}{1+\lambda} + \frac{30}{2+\lambda} - \frac{210}{3+\lambda} + \frac{560}{4+\lambda} - \frac{630}{5+\lambda} + \frac{252}{6+\lambda}\biggr]}

and \eqn{\tau_4 = \lambda_4 / \lambda_2} and \eqn{\tau_6 = \lambda_6 / \lambda_2}. The Tukey Lambda is not otherwise present or implemented in the \pkg{lmomco} package. It is provided here for theoretical completeness.
}
\author{W.H. Asquith}
\seealso{\code{\link{plotlmrdia46}}, \code{\link{lmrdia}}}
\examples{
lratios <- lmrdia46()
}
\keyword{L-moment ratio diagram}
\keyword{L-moment ratio diagram (Tau4 and Tau6)}
