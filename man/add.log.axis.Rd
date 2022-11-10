\name{add.log.axis}
\alias{add.log.axis}
\title{Add a Polished Logarthimic Axis to a Plot}
\description{
This function provides special support for adding superior looking base-10 logarithmic axes relative to \bold{R}'s defaults, which are an embarassment. The \bold{Examples} section shows an overly elaborate version made by repeated calls to this function with a drawback that each call redraws the line of the axis so deletion in editing software might be required. This function is indexed under the \dQuote{lmomco functions} because of its relation to \code{\link{add.lmomco.axis}} and is not named \code{add.lmomcolog.axis} because such a name is too cumbersome.
}
\usage{
add.log.axis(make.labs=FALSE, logs=c(2, 3, 4, 5, 6, 7, 8, 9), side=1,
             two.sided=FALSE, label=NULL, x=NULL, col.ticks=1, ...)
}
\arguments{
  \item{make.labs}{A logical controlling whether the axis is labled according to the values in \code{logs}.}
  \item{logs}{A numeric vector of log-cycles for which ticking and (or) labeling is made. These are normalized to the first log-cycle, so a value of \eqn{3} would spawn values such as \eqn{\cdots, 0.03, 0.3, 3, 30, \cdots} through a range exceeding the axis limits. The default anticipates that a second call to the function will be used to make longer ticks at the even log-cycles; hence, the value 1 is not in the default vector. The \bold{Examples} section provides a thorough demonstration.}
  \item{side}{An integer specifying which side of the plot the axis is to be drawn on, and argument corresponds the axis side argument of the \code{axis} function. The axis is placed as follows: 1=below, 2=left, 3=above, and 4=right.}
  \item{two.sided}{A logical controlling whether the side oppose of \code{side} also is to be drawn.}
  \item{label}{The label (title) of the axis, which is placed by a call to function \code{mtext}, and thus either the \code{xlab} or \code{ylab} arguments for \code{plot} should be set to the empty string \code{""}.}
  \item{x}{This is an optional data vector (untransformed!), which will compute nice axis limits and return them. These limits will align with (snap to) the integers within a log10-cycle.}
  \item{col.ticks}{This is the same argument as the \code{axis} function.}
  \item{...}{Additional arguments to pass to \code{axis}.}
}
\value{
   No value is returned, except if argument \code{x} is given, for which nice axis limits are returned. By overall design, this function is used for its side effects.
}
\author{W.H. Asquith}
\seealso{\code{\link{add.lmomco.axis}} }
\examples{
\dontrun{
par(mgp=c(3,0.5,0)) # going to tick to the inside, change some parameters
X <- 10^sort(rnorm(65)); pp <- pp(X) # generate synthetic data
ylim <- add.log.axis(x=X) # snap to some nice integers within the cycle
plot(qnorm(pp), X, xaxt="n", yaxt="n", xlab="", ylab="", log="y",
     xlim=c(-2,3), ylim=ylim, pch=6, yaxs="i", col=4)
add.lmomco.axis(las=2, tcl=0.5, side.type="RI", otherside.type="NPP")
# Logarithmic axis: the base ticks to show logarithms
add.log.axis(side=2,      tcl=0.8*abs(par()$tcl), two.sided=TRUE)
#                   the long even-cycle tick, set to inside and outside
add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=2, two.sided=TRUE)
#                   now a micro tick at the 1.5 logs but only on the right
add.log.axis(logs=c(1.5), tcl=+0.5*abs(par()$tcl), side=4)
#                   and only label the micro tick at 1.5 on the right
add.log.axis(logs=c(1.5), side=4, make.labs=TRUE, las=3) # but weird rotate
#                   add the bulk tick labeling and axis label.
add.log.axis(logs=c(1, 2, 3, 4, 6), side=2, make.labs=TRUE, las=1, label="QUANTILE")
par(mgp=c(3,1,0)) # restore defaults}
}
\keyword{utility (logarithmic)}
\keyword{axis (utility)}
\keyword{Graphics}
\concept{logarithmic axis}
\concept{log axis}
\keyword{The lmomco functions}
