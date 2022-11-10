\encoding{UTF-8}
\name{tlmr2par}
\alias{tlmr2par}
\title{Sample Trimmed L-moments to Fitted Distribution}
\description{
Parameter estimation of a distribution given initial estimate of the parameters of the distribution to the sample trimmed L-moment (TL-moment) using numerical optimization. Thought the TL-moments can be used with substantial depth into either tail and need not be symmetrically trimmed, the TL-moments do not appear as useful when substantial tail trimming is needed, say for mix population mitigation. Then censored or truncation methods might be preferred. The \code{\link{x2xlo}} family of operations can be used for conditional left-tail truncation, which is not uncommon in frequency analyses of rail-tail interest water resources phenomena.
}
\usage{
tlmr2par(x, type, para.int=NULL, trim=NULL, leftrim=NULL, rightrim=NULL, ...)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{type}{Three character (minimum) distribution type (for example, \code{type="gev"}, see \code{\link{dist.list}}.}
  \item{para.int}{Initial parameters as a vector \eqn{\Theta} or as an \pkg{lmomco} parameter \dQuote{object} from say \code{\link{vec2par}}. If a vector is given, then internally \code{\link{vec2par}} is called with distribution equal to \code{type}.}
  \item{trim}{Level of symmetrical trimming to use in the computations. Although \code{NULL} is in the argument list, the default is 0---the usual L-moment is returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
  \item{rightrim}{Level of trimming of the right-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
  \item{...}{Other arguments to pass to the \code{optim()} function.}
}
\value{
  An \R \code{list} is returned.  This list should contain at least the following items, but some distributions such as the \code{revgum} have extra.

  \item{type}{The type of distribution in three character (minimum) format.}
  \item{para}{The parameters of the distribution.}
  \item{text}{Optional material. If the solution fails but the optimization appears to converge, then this element is inserted into the list and the \code{para} will be all \code{NA}.}
  \item{source}{Attribute specifying source of the parameters.}
  \item{rt}{The list from the \code{optim()} function.}
  \item{para.int}{A copy of the initial parameters given.}
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{W.H. Asquith}
\seealso{\code{\link{theoTLmoms}}, \code{\link{TLmoms}}, \code{\link{lmr2par}}}
\examples{
\donttest{
# (1) An example to check that trim(0,0) should recover whole sample
the.data <- rlmomco(140, vec2par(c(3, 0.4, -0.1), type="pe3"))
wild.guess <- vec2par(c(mean(the.data), 1, 0),    type="pe3")
pe3whole <- lmom2par(lmoms(the.data),             type="pe3")
pe3trimA  <- tlmr2par(the.data, "pe3", para.int=wild.guess, leftrim=0,  rightrim=0)
pe3trimB  <- tlmr2par(the.data, "pe3", para.int=wild.guess, leftrim=10, rightrim=3)
message("PE3 parent       = ", paste0(pe3whole$para, sep=" "))
message("PE3 whole sample = ", paste0(pe3whole$para, sep=" "))
message("PE3 trim( 0, 0)  = ", paste0(pe3trimA$para, sep=" "))
message("PE3 trim(10, 3)  = ", paste0(pe3trimB$para, sep=" ")) #}

\donttest{
# (2) An example with "real" outliers
FF <- lmomco::nonexceeds(); qFF <- qnorm(FF); type <- "gev"
the.data <- c(3.064458, 3.139879, 3.167317, 3.225309, 3.324282, 3.330414,
             3.3304140, 3.340444, 3.357935, 3.376577, 3.378398, 3.392697,
             3.4149730, 3.421604, 3.424882, 3.434569, 3.448706, 3.451786,
             3.4517860, 3.462398, 3.465383, 3.469822, 3.491362, 3.501059,
             3.5224440, 3.523746, 3.527630, 3.527630, 3.531479, 3.546543,
             3.5932860, 3.597695, 3.600973, 3.614897, 3.620136, 3.660865,
             3.6848450, 3.820858, 4.708421)
the.data <- sort(the.data) # though already sorted, backup for plotting needs

# visually, looks like 4 outliers to the left and one outlier to the right
# perhaps the practical situation is that we do not wan the left tail to
# mess up the right when fitting a distribution because maybe the practical
# aspects are the that right tail is of engineering interest, but then we
# have some idea that the one very large event is of questionable suitability
t1 <- 4; t2 <- 1 # see left and right trimming and then estimation parameters
whole.para <- lmom2par(lmoms(the.data), type=type)
trim.para  <- tlmr2par(the.data, type, para.int=whole.para, leftrim=t1, rightrim=t2)

n <- length(the.data)
cols <- rep(grey(0.5), n)
pchs <- rep(1, n)
if(t1 != 0) {
  cols[      1 :t1] <- "red"
  cols[(n-t2+1):n ] <- "purple"
}
if(t2 != 0) {
  pchs[      1 :t1] <- 16
  pchs[(n-t2+1):n ] <- 16
}
plot( qFF, qlmomco(FF, whole.para), type="l", lwd=2, ylim=c(3.1,4.8),
           xlab="Standard normal variate",
           ylab="Some phenomena, log10(cfs)")
lines(qFF, qlmomco(FF, trim.para), col=4, lwd=3)
points(qnorm(pp(the.data)), sort(the.data), pch=pchs, col=cols)
legend("topleft", c("L-moments",
                   paste0("TL-moments(", t1, ",", t2,")")), bty="n",
                  lty=c(1,1), lwd=c(2,3), col=c(1,4))
# see the massive change from the whole sample to the trim(t1,t2) curve}
}
\keyword{L-moment (sample)}
\keyword{L-moment (trimmed sample)}
\keyword{Ensemble list of parameter estimation methods}
\keyword{Ensemble list of major parameter estimation methods}
