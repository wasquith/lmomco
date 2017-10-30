\name{check.pdf}
\alias{check.pdf}
\title{ Check and Potentially Graph Probability Density Functions }
\description{
This convenience function checks that a given probability density function (\code{pdf}) from \pkg{lmomco} appears to numerically be valid. By definition a \code{pdf} function must integrate to unity. This function permits some flexibility in the limits of integration and provides a high-level interface from graphical display of the \code{pdf}.
}
\usage{
check.pdf(pdf, para, lowerF=0.001, upperF=0.999,
               eps=0.02, verbose=FALSE, plot=FALSE, plotlowerF=0.001,
               plotupperF=0.999, ...)
}
\arguments{
  \item{pdf}{A probability density function from \pkg{lmomco}.}
  \item{lowerF}{The lower bounds of nonexceedance probability for the numerical integration.}
  \item{upperF}{The upper bounds of nonexceedance probability for the numerical integration.}
  \item{para}{The parameters of the distribution.}
  \item{eps}{An error term expressing allowable error (deviation) of the numerical integration from unity. (If that is the objective of the call to the \code{\link{check.pdf}} function.)}
  \item{verbose}{Is verbose output desired?}
  \item{plot}{Should a plot (polygon) of the \code{pdf} integration be produce?}
  \item{plotlowerF}{Alternative lower limit for the generation of the curve depicting the \code{pdf} function.}
  \item{plotupperF}{Alternative upper limit for the generation of the curve depicting the \code{pdf} function.}
  \item{...}{Additional arguments that are passed onto the \R function \code{integration} function.}
}
\value{
An \R \code{list} structure is returned
  \item{isunity}{Given the \code{eps} is \code{F} close enough.}
  \item{F}{The numerical integration of \code{pdf} from \code{lowerF} to \code{upperF}.}
}
\author{ W.H. Asquith}
\examples{
lmrg <- vec2lmom(c( 100, 40, 0.1)) # Arbitrary L-moments
lmrw <- vec2lmom(c(-100, 40,-0.1)) # Reversed Arbitrary L-moments
gev  <- pargev(lmrg) # parameters of Generalized Extreme Value distribution
wei  <- parwei(lmrw) # parameters of Weibull distribution

# The Weibull is a reversed GEV and plots in the following examples show this.
# Two examples that should integrate to "unity" given default parameters.
layout(matrix(c(1,2), 2, 2, byrow = TRUE), respect = TRUE)
check.pdf(pdfgev,gev,plot=TRUE)
check.pdf(pdfwei,wei,plot=TRUE)
}
\keyword{utility (distribution)}
\keyword{Graphics}
