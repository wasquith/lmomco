\name{quape3}
\alias{quape3}
\title{Quantile Function of the Pearson Type III Distribution}
\description{
This function computes the quantiles of the Pearson Type III distribution given parameters (\eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma})  computed by \code{\link{parpe3}}. The quantile function  has no explicit form (see \code{\link{cdfpe3}}).

For the implementation in the \pkg{lmomco} package, the three parameters are \eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma} for the mean, standard deviation, and skew, respectively. Therefore, the Pearson Type III distribution is of considerable theoretical interest to this package because the parameters, which are estimated via the L-moments, are in fact the product moments, although, the values fitted by the method of L-moments will not be numerically equal to the sample product moments. Further details are provided in the Examples section under \code{\link{pmoms}}.
}
\usage{
quape3(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parpe3}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfpe3}}, \code{\link{pdfpe3}}, \code{\link{lmompe3}}, \code{\link{parpe3}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  quape3(0.5,parpe3(lmr))

\dontrun{
  # Let us run an experiment on the reflection symmetric PE3.
  # Pick some parameters suitable for hydrologic applications in log.
  para_neg <- vec2par(c(3,.3,-1), type="pe3") # Notice only the
  para_pos <- vec2par(c(3,.3,+1), type="pe3") # sign change of skew.

  nsim <- 1000 # Number of simulations
  nsam <- 70   # Reasonable sample size in hydrology
  neg <- pos <- rep(NA, nsim)
  for(i in 1:nsim) {
    ff <- runif(nsam) # Ensure that each qlmomco()-->quape3() has same probs.
    neg[i] <- lmoms.cov(qlmomco(ff, para_neg), nmom=3, se="lmrse")[3]
    pos[i] <- lmoms.cov(qlmomco(ff, para_pos), nmom=3, se="lmrse")[3]
    # We have extracted the sample standard error of L-skew from the sample
    # This is not the same as the standard error of so computed PE3 
    # parameters, but for the illustration here, it does not matter much.
  }
  zz <- data.frame(setau3=c(neg,pos), # preserve to make grouping boxplot
                   sign=c(rep("negskew", nsim), rep("posskew", nsim)))
  boxplot(zz$setau3~zz$sign, xlab="Sign of a '1' PE3 skew",
                             ylab="Standard error of L-skew")
  mtext("Standard Errors of 1,000 PE3 Parents (3,0.3,+/-1) (n=70)")
  # Notice that the distribution of the standard errors of L-skew are 
  # basically the same whether or no the sign of the skew is reversed.
  # Finally, we make a scatter plot as a check that for any given sample
  # derived from same probabilities that the standard errors are indeed,
  # that is, remain sample specific.
  plot(neg, pos, xlab="Standard error of -1 skew simulation",
                 ylab="Standard error of +1 skew simulation")
  mtext("Standard Errors of 1,000 PE3 Parents (3,0.3,+/-1) (n=70)") # }
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Pearson Type III}
