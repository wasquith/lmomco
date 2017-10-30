\name{are.pargld.valid}
\alias{are.pargld.valid}
\title{Are the Distribution Parameters Consistent with the Generalized Lambda Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfgld}}, \code{\link{pdfgld}}, \code{quagld}, and \code{\link{lmomgld}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively.  These functions internally use the \code{\link{are.pargld.valid}} function.
}
\usage{
are.pargld.valid(para, verbose=FALSE, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{pargld}} or \code{\link{vec2par}}.}
  \item{verbose}{A logical switch on additional output to the user---default is \code{FALSE}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\details{
Karian and Dudewicz (2000) outline valid parameter space of the Generalized Lambda distribution. First, according to Theorem 1.3.3 the distribution is valid if and only if
\deqn{\alpha(\kappa F^{\kappa - 1} + h(1-F)^{h -1 }) \ge 0 \mbox{.}}
for all \eqn{F \in [0,1]}. The \code{are.pargld.valid} function tests against this condition by incrementing through \eqn{[0,1]} by \eqn{dF = 0.0001}. This is a brute force method of course. Further, Karian and Dudewicz (2002) provide a diagrammatic representation of regions in \eqn{\kappa} and \eqn{h} space for suitable \eqn{\alpha} in which the distribution is valid. The \code{\link{are.pargld.valid}} function subsequently checks against the 6 valid regions as a secondary check on Theorem 1.3.3. The regions of the distribution are defined for suitably choosen \eqn{\alpha} by

\deqn{\mbox{Region 1:  } \kappa \le -1 \mbox{ and } h \ge 1 \mbox{,}}
\deqn{\mbox{Region 2:  } \kappa \ge 1 \mbox{ and } h \le -1 \mbox{,}}
\deqn{\mbox{Region 3:  } \kappa \ge 0 \mbox{ and } h \ge 0 \mbox{,}}
\deqn{\mbox{Region 4:  } \kappa \le 0 \mbox{ and } h \le 0 \mbox{,}}
\deqn{\mbox{Region 5:  } h \ge (-1/\kappa) \mbox{ and } -1 \ge \kappa \le 0 \mbox{, and}}
\deqn{\mbox{Region 6:  } h \le (-1/\kappa) \mbox{ and } h \ge -1 \mbox{ and } \kappa \ge 1 \mbox{.}}
}
\value{
  \item{TRUE}{If the parameters are \code{gld} consistent.}
  \item{FALSE}{If the parameters are not \code{gld} consistent.}
}
\references{
Asquith, W.H., 2007, L-moments and TL-moments of the generalized lambda distribution: Computational Statistics and Data Analysis, v. 51, no. 9, pp. 4484--4496.

Karian, Z.A., and Dudewicz, E.J., 2000, Fitting statistical distributions---The generalized lambda distribution and generalized bootstrap methods: CRC Press, Boca Raton, FL, 438 p.
}
\note{
This function calls \code{\link{is.gld}} to verify consistency between the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.gld}}, \code{\link{pargld}} }
\examples{
\dontrun{
para <- vec2par(c(123,34,4,3),type='gld')
if(are.pargld.valid(para)) Q <- quagld(0.5,para)

# The following is an example of inconsistent L-moments for fitting but
# prior to lmomco version 2.1.2 and untrapped error was occurring.
lmr <- lmoms(c(33, 37, 41, 54, 78, 91, 100, 120, 124))
para <- pargld(lmr); are.pargld.valid(para)}
}
\keyword{utility (distribution)}
\keyword{Distribution: Generalized Lambda}
\keyword{utility (distribution, parameter validation)}

