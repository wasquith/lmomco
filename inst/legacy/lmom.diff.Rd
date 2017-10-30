\name{lmom.diff}
\alias{lmom.diff}
\title{ Difference Between L-moments of the Distribution and the L-moments of the Data }
\description{
This function computes the difference between the L-moments derived from a
parameterized distribution and the L-moments as computed from the data.  This
function is useful to characterize the bias that develops between the theoretical
L-moments of a distribution and the L-moments of the data.  This function also
is an important test on the algorithms that fit distributions to the L-moments.
The difference is computed as the L-moment from the distribution minus the
L-moment of the data. The \code{\link{lmorph}} function is used internally to get the L-moment objects into the appropriate format.
}
\usage{
lmom.diff(lmomparm, lmomdata, verbose=TRUE, digits=4)
}
\arguments{
  \item{lmomparm}{L-moments of a distribution such as from \code{par2lmom}}
  \item{lmomdata}{L-moments of the data such as from \code{lmom.ub}}
  \item{verbose}{Logical switch on verbosity of output.  Default is \code{TRUE}.}
  \item{digits}{Number of digits to pass internally to the \code{signif} function for rounding of results.}
}
\value{

\preformatted{"THE FIVE DIFFERENCES BETWEEN L-MOMENTS OF DISTRIBUTION AND DATA"}

\preformatted{  L1diff    L2diff  T3diff  T4diff  T5diff}

\preformatted{1      0 -1.11e-16 -0.7345 -0.4005 -0.3942}

where the five values are the differences between the theoretical L-moments of the fitted distribution and
the sample L-moments of the data (theoretical minus sample) in the titled column. These values are also returned as part of an R \code{list}. If a \code{list} element has \code{numeric(0)} then likely one of the L-moments is \code{NULL} or \code{NA} for the distribution such as the results reported for the example involving the Generalized Logistic distribution (\code{glo}).

A dataframe of the difference is the returned value.

}
\author{ W.H. Asquith}
\seealso{ \code{\link{par2lmom}}, \code{\link{lmom2par}}, \code{\link{lmorph}}}
\examples{
  # The first three moment differences are zero because the GLO is only 
  # fit to these and not the higher moments.
  lmr <- lmom.ub(rnorm(40))
  para <- lmom2par(lmr, type = 'glo')
  lmom.diff(par2lmom(para),lmr)

  # The first two moment differences are zero because the Normal is only 
  # fit to these and not the higher moments.
  lmr <- lmom.ub(rnorm(40))
  lmr <- lmorph(lmr)
  para <- parnor(lmr)
  lmom.diff(lmomnor(para),lmr)

}
\keyword{distribution}
