\name{plotradarlmr}
\alias{plotradarlmr}
\title{Plot L-moment Radar Plot (Chart) Graphic}
\description{
Plot a L-moment radar plots (charts).  This graphic is somewhat experimental and of unknown application benefit as no known precedent seems available. L-moment ratio diagrams (\code{\link{plotlmrdia}}) are incredibly useful but have generally been restricted to the 2-D domain. The graphic supported here attempts to provide a visualization of \eqn{\tau_r} for an arbitrary \eqn{(r-2) > 3} number of axes in the form of a radar plot. The angle of the axes is uninformative but the order of the axes is for \eqn{\tau_r} for \eqn{r = 3, 4, \cdots}. The radar plot is essentially a line graph but mapped to a circular space at the expense of more \emph{ink} being used. The radar plot is primarily intended to be a mechansim in \pkg{lmomco} for which similarity between other radar plots or presence of outlier combinations of \eqn{\tau_r} can be judged when seen amongst various samples.
}
\usage{
plotradarlmr(lmom, num.axis=4, plot=TRUE, points=FALSE, poly=TRUE, tag=NA,
             title="L-moment Ratio Radar Plot", make.zero.axis=FALSE,
             minrat=NULL, maxrat=NULL, theomins=TRUE, rot=0,
             labadj=1.2, lengthadj=1.75, offsetadj=0.25, scaleadj=2.2,
     axis.control  = list(col=1, lty=2, lwd=0.5, axis.cex=0.75, lab.cex=0.95),
     point.control = list(col=8, lwd=0.5, pch=16),
     poly.control  = list(col=rgb(0,0,0,.1), border=1, lty=1, lwd=1), ...)
}
\arguments{
  \item{lmom}{L-moment object such as from \code{\link{lmoms}}.}
  \item{num.axis}{The number of axes. Some error trapping in axis count relative to the length of the \eqn{\tau_r} in \code{lmom} is made.}
  \item{plot}{A logical controlling whether \R function \code{plot} will be called.}
  \item{points}{A logical controlling whether the points of defined by the \eqn{\tau_r} in \code{lmom}.}
  \item{poly}{A logical controlling whether the polygon of defined by the \eqn{\tau_r} in \code{lmom}.}
  \item{tag}{A text tag plotted at the center of the plot. An \code{NA} will result in nothing being plotted.}
  \item{title}{The title of the plot. An \code{NA} will result in nothing being plotted.}
  \item{make.zero.axis}{A logical controlling whether polygon will be \dQuote{faked in} like as if \eqn{\tau_r} having all zeros are provided. This feature is to act as a mechanism to overlay only the zero axis such as might be needed when a lot of other material has been already been drawn on the plot.}
  \item{minrat}{A vector of the minimum values for the \eqn{\tau_r} axes in case the user desired to have some zoomability. The default is all \eqn{-1} values, and a scalar for \code{minrat} will be repeated for the \code{num.axis}.}
  \item{maxrat}{A vector of the maximum values for the \eqn{\tau_r} axes in case the user desired to have some zoomability. The default is all \eqn{+1} values, and a scalar for \code{maxrat} will be repeated for the \code{num.axis}.}
  \item{theomins}{The are some basic and fundamental lower limits other than -1 that if used provide for a better relative scaling of the axes on the plot. If \code{TRUE}, then some select overwritting of potential user-provided \code{minrat} is provided.}
  \item{rot}{The basic rotational offset for the angle of the first (\eqn{\tau_3}) axis.}
  \item{labadj}{An adjustment multiplier to help positions of the axis titles.}
  \item{lengthadj}{An adjustment multiplier characterize axis length.}
  \item{offsetadj}{An adjustment to help set the empty space in the middle of the plot for the \code{tag}.}
  \item{scaleadj}{An adjustment multiplier to help set the parent domain of the underlying (but hidden) x-y plot called by the \R function \code{plot}.}
  \item{axis.control}{A specially built and not error trapped \R \code{list} to hold the control elements of the axes.}
  \item{point.control}{A specially built and not error trapped \R \code{list} to hold the control elements for plotting of the points if \code{points=TRUE}.}
  \item{poly.control}{A specially built and not error trapped \R \code{list} to hold the control elements for plotting of the polygon if \code{poly=TRUE}.}
  \item{...}{Additional arguments passed on to the \R function \code{text} function for the \code{title} and \code{tag}. This argument is largely not intended for general use, unlike most idioms of \code{...} in \R, but is provided at the release of this function to help developers and avoid future backwards compatibility problems.}
}
\note{
This function has many implicit flexible features. The example below attempts to be reasonably comprehensive. Note that in the example that it is required to continue \dQuote{knowing} what \code{minrat} and \code{maxrat} where used with \code{plot=TRUE}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{plotlmrdia}}}
\examples{
\dontrun{
plotradarlmr(NULL, minrat=-0.6, maxrat=0.6, tag="2 GEVs") # create the plot base
gev  <- vec2par(c(1230,123,-.24), type="gev") # set first parent distribution
poly <- list(col=NA, border=rgb(0,0,1,.1))    # set up polygon handling (blue)
for(i in 1:100) { # perform 100 simulations of the GEV with a sample of size 36
   plotradarlmr(lmoms(rlmomco(36,gev), nmom=6), plot=FALSE,
                poly.control=poly, minrat=-0.6, maxrat=0.6)
}
poly <- list(col=NA, border=4, lwd=3) # set up parent polygon
plotradarlmr(theoLmoms(gev, nmom=6), plot=FALSE,
             poly.control=poly, minrat=-0.6, maxrat=0.6) # draw the parent
 gev <- vec2par(c(450,1323,.5), type="gev") # set second parent distribution
poly <- list(col=NA, border=rgb(0,1,0,.1))  # set up polygon handling (green)
for(i in 1:100) { # perform 100 simulations of the GEV with a sample of size 36
   plotradarlmr(lmoms(rlmomco(36,gev), nmom=6), plot=FALSE,
                poly.control=poly, minrat=-0.6, maxrat=0.6) # draw the parent
}
poly <- list(col=NA, border=3, lwd=3) # set up parent polygon
plotradarlmr(theoLmoms(gev, nmom=6), plot=FALSE,
             poly.control=poly, minrat=-0.6, maxrat=0.6)
poly <- list(col=NA, border=6, lty=1, lwd=2) # make the zeros purple to standout.
plotradarlmr(NULL, make.zero.axis=TRUE, plot=FALSE,
             poly.control=poly, minrat=-0.6, maxrat=0.6)
}
}
\keyword{L-moment radar plot}
\keyword{Graphics}
\concept{radar chart}
\concept{radar graphic}
