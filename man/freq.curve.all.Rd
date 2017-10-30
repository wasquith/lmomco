\name{freq.curve.all}
\alias{freq.curve.all}
\title{ Compute Frequency Curve for Almost All Distributions}
\description{
This function is dispatcher on top of a select suite of \code{quaCCC} functions that compute frequency curves for the L-moments. The term \dQuote{frequency curves} is common in hydrology and is a renaming of the more widenly known by statisticians term the \dQuote{quantile function.} The notation \code{CCC} represents the character notation for the distribution: \code{exp}, \code{gam}, \code{gev}, \code{gld}, \code{glo}, \code{gno}, \code{gpa}, \code{gum}, \code{kap}, \code{nor}, \code{pe3}, \code{wak}, and \code{wei}. The nonexceedance probabilities to construct the curves are derived from \code{\link{nonexceeds}}.
}
\usage{
freq.curve.all(lmom, aslog10=FALSE, asprob=TRUE,
                     no2para=FALSE, no3para=FALSE,
                     no4para=FALSE, no5para=FALSE,
                     step=FALSE, show=FALSE,
                     xmin=NULL, xmax=NULL, xlim=NULL,
                     ymin=NULL, ymax=NULL, ylim=NULL,
                     aep4=FALSE, exp=TRUE, gam=TRUE, gev=TRUE, gld=FALSE,
                     glo=TRUE, gno=TRUE, gpa=TRUE, gum=TRUE, kap=TRUE,
                     nor=TRUE, pe3=TRUE, wak=TRUE, wei=TRUE,...)
}
\arguments{
  \item{lmom}{A L-moment object from \code{\link{lmoms}}, \code{\link{lmom.ub}}, or \code{\link{vec2lmom}}.}
  \item{aslog10}{Compute \code{log10} of quantiles---note that \preformatted{NaNs produced in: log(x, base)} will be produced for less than zero values.}
  \item{asprob}{The \R \code{qnorm} function is used to convert nonexceedance probabilities, which are produced by \code{\link{nonexceeds}}, to standard normal variates. The Normal distribution will plot as straight line when this argument is \code{TRUE}}
  \item{no2para}{If \code{TRUE}, do not run the 2-parameter distributions: \code{exp}, \code{gam}, \code{gum}, and \code{nor}.}
  \item{no3para}{If \code{TRUE}, do not run the 3-parameter distributions: \code{gev}, \code{glo}, \code{gno}, \code{gpa}, \code{pe3}, and \code{wei}.}
  \item{no4para}{If \code{TRUE}, do not run the 4-parameter distributions: \code{kap}, \code{gld}, \code{aep4}.}
  \item{no5para}{If \code{TRUE}, do not run the 5-parameter distributions: \code{wak}.}
  \item{step}{Shows incremental processing of each distribution.}
  \item{show}{Plots all the frequency curves in a simple (crowded) \code{plot}.}
  \item{xmin}{Minimum x-axis value to use instead of the automatic value determined from the nonexceedance probabilities. This argument is only used is \code{show=TRUE}.}
  \item{xmax}{Maximum x-axis value to use instead of the automatic value determined from the nonexceedance probabilities. This argument is only used is \code{show=TRUE}.}
  \item{xlim}{Both limits of the x-axis. This argument is only used is \code{show=TRUE}.}
  \item{ymin}{Minimum y-axis value to use instead of the automatic value determined from the nonexceedance probabilities. This argument is only used is \code{show=TRUE}.}
  \item{ymax}{Maximum y-axis value to use instead of the automatic value determined from the nonexceedance probabilities. This argument is only used is \code{show=TRUE}.}
  \item{ylim}{Both limits of the y-axis. This argument is only used is \code{show=TRUE}.}
  \item{aep4}{A logical switch on computation of corresponding distribution---default is \code{FALSE}.}
  \item{exp}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{gam}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{gev}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{gld}{A logical switch on computation of corresponding distribution---default is \code{FALSE}.}
  \item{glo}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{gno}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{gpa}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{gum}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{kap}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{nor}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{pe3}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{wak}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{wei}{A logical switch on computation of corresponding distribution---default is \code{TRUE}.}
  \item{...}{Additional parameters are passed to the parameter estimation routines such as \code{parexp}.}
}
\value{
  An extensive \R \code{data.frame} of frequency curves. The nonexceedance probability values, which are provided by \code{\link{nonexceeds}}, are the first item in the \code{data.frame} under the heading of \code{\link{nonexceeds}}. If a particular distribution could not be fit to the L-moments of the data; this particular function returns zeros.
}
\note{
	The distributions selected for this function represent a substantial fraction of, but not all, distributions supported by \pkg{lmomco}. The \code{all} and \dQuote{all} in the function name and the title of this documentation is a little misleading. The selection process was made near the beginning of \pkg{lmomco} availability and distributions available in the earliest versions. Further the selected distributions are frequently encountered in hydrology and because these are also those considered in length by Hosking and Wallis (1997) and the \pkg{lmom} package.
}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{
\code{\link{quaaep4}},
\code{\link{quaexp}},
\code{\link{quagam}},
\code{\link{quagev}},
\code{\link{quagld}},
\code{\link{quaglo}},
\code{\link{quagno}},
\code{\link{quagpa}},
\code{\link{quagum}},
\code{\link{quakap}},
\code{\link{quanor}},
\code{\link{quape3}},
\code{\link{quawak}}, and
\code{\link{quawei}}.
}
\examples{
L <- vec2lmom(c(35612,23593,0.48,0.21,0.11))
Qtable1 <- freq.curve.all(L, step=TRUE, no2para=TRUE, no4para=TRUE)
\dontrun{
Qtable2 <- freq.curve.all(L, gld=TRUE, show=TRUE)}
}
\keyword{misc}
