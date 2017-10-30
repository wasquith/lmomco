\name{add.lmomco.axis}
\alias{add.lmomco.axis}
\title{Add an lmomco Axis to a Plot}
\description{
This function provides special support for adding probability-like axes to an existing plot. The function supports a recurrence interval (RI) axis, normal probability axis (NPP), and standard normal variate (SNV) axis. The function is built around the interface model that standard normal transformation of the values for the respective axis controlled by this function are being plotted; this means that \code{qnorm()} should be wrapped on the values of nonexceedance probability. This is an ease oversight to make (see Examples section below and note use of \code{qnorm(pp)}).

The function provides a convenient interface for labeling and titling two axes, so adjustments to default margins might be desired. The pertinent control is achieved using the \code{par()} function, which might be of the form \code{par(mgp=c(3,0.5,0), mar=c(5,4,4,3))} say for plotting the \pkg{lmomco} axis both on the left and right (see \code{\link{z.par2cdf}} for an example).
}
\usage{
add.lmomco.axis(side=1, twoside=FALSE, side.type=c("NPP", "RI", "SNV"),
                otherside.type=c("NA", "RI", "SNV", "NPP"),
                alt.lab=NA, alt.other.lab=NA, npp.as.aep=FALSE,
                NPP.control=NULL, RI.control=NULL, SNV.control=NULL, ...)
}
\arguments{
  \item{side}{The side of the plot (1=bottom, 2=left, 3=top, 4=right).}
  \item{twoside}{A logical triggering whether the tick marks are echoed on the opposite side. This value is forced to \code{FALSE} if \code{otherside.type} is not \code{"NA"}.}
  \item{side.type}{The axis type for the primary \code{side}.}
  \item{otherside.type}{The optional axis type for the opposite side. The default is a literal not applicable.}
    \item{alt.lab}{A short-cut to change the axis label without having to specify a \code{*.control} argument and its \code{label} attribute. The label attribute of \code{alt.lab} is not \code{NA} is used instead of the defaults. This argument overrides behavior of the \code{otherside.type} labeling so use of \code{alt.lab} only makes sense if \code{otherside.type} is left as \code{NA}.}
   \item{alt.other.lab}{Similar to \code{alt.lab} but can house an alternative label (see \bold{Examples}.}
  \item{npp.as.aep}{Convert nonexceedance probability to exceedance probability, which is a que for \code{alt.other.lab} and nonexceedance probabilities are changed by \eqn{1-F}, but the real coordinates for plotting remain in the nonexceedance probability context.}
  \item{NPP.control}{An optional \R \code{list} used to influence the NPP axis.}
  \item{RI.control}{An optional \R \code{list} used to influence the RI axis.}
  \item{SNV.control}{An optional \R \code{list} used to influence the SNV axis.}
   \item{...}{Additional arguments that are passed to the \R function \code{Axis}.}
}
\value{
   No value is returned. This function is used for its side effects.
}
\note{
The \code{NPP.control}, \code{RI.control}, and \code{SNV.control} are \R \code{list} structures that can be populated (and perhaps someday extended) to feed various settings into the respective axis types.  In brief:

The \code{NPP.control} provides
\tabular{ll}{
  \code{label} \tab The title for the NPP axis---be careful with value of \code{as.exceed}. \cr
  \code{probs} \tab A vector of nonexceedance probabilities \eqn{F}. \cr
  \code{probs.lab} \tab A vector of nonexceedance probabilities \eqn{F} to label. \cr
  \code{digits} \tab The digits for the \R function \code{format} to enhance appearance. \cr
  \code{line}  \tab The line for the \R function \code{mtext} to place \code{label}. \cr
  \code{as.exceed} \tab A logical triggering \eqn{S = 1 - F}.
}

The \code{RI.control} provides
\tabular{ll}{
  \code{label} \tab The title for the RI axis. \cr
  \code{Tyear} \tab A vector of \eqn{T}-year recurrence intervals. \cr
  \code{line}  \tab The line for the \R function \code{mtext} to place \code{label}.
}

The \code{SNV.control} provides
\tabular{ll}{
  \code{label} \tab The title for the SNV axis. \cr
  \code{begin} \tab The beginning \dQuote{number of standard deviations}. \cr
  \code{end}   \tab The ending \dQuote{number of standard deviations}. \cr
  \code{by}    \tab The step between \code{begin} and \code{end}. \cr
  \code{line}  \tab The line for the \R function \code{mtext} to place \code{label}.
}

The user is responsible for appropriate construction of the \code{control} lists. Very little error trapping is made to keep the code base tight. The defaults when the function definition are likely good for many types of applications. Lastly, the manipulation of the \code{mgp} parameter in the example is to show how to handle the offset between the numbers and the ticks when the ticks are moved to pointing inward, which is opposite of the default in \R.
}
\author{ W.H. Asquith}
\seealso{\code{\link{prob2T}}, \code{\link{T2prob}}, \code{\link{add.log.axis}} }
\examples{
par(mgp=c(3,0.5,0)) # going to tick to the inside, change some parameters
X <- sort(rnorm(65)); pp <- pp(X) # generate synthetic data
plot(qnorm(pp), X, xaxt="n", xlab="", ylab="QUANTILE", xlim=c(-2,3))
add.lmomco.axis(las=2, tcl=0.5, side.type="RI", otherside.type="NPP")
par(mgp=c(3,1,0)) # restore defaults

\dontrun{
opts <- options(scipen=6); par(mgp=c(3,0.5,0))
X <- sort(rexp(65, rate=.0001))*100; pp <- pp(X) # generate synthetic data
plot(qnorm(pp), X, yaxt="n", xaxt="n", xlab="", ylab="", log="y")
add.log.axis(side=2,    tcl=+0.8*abs(par()$tcl),         two.sided=TRUE)
add.log.axis(logs=c(1), tcl=-0.5*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1), tcl=+1.3*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=1:8, side=2, make.labs=TRUE, las=1, label="QUANTILE")
add.lmomco.axis(las=2, tcl=0.5, side.type="NPP", npp.as.aep=TRUE)
options(opts)
par(mgp=c(3,1,0)) # restore defaults}
}
\keyword{utility (probability)}
\keyword{axis (utility)}
\keyword{Graphics}
\concept{probability axis}
\concept{normal probability paper}
\concept{probability paper}
\concept{recurrence interval}
\concept{return period}
\concept{standard normal axis}
\concept{standard normal paper}
\concept{standard normal variate axis}
\concept{standard normal variate paper}
\keyword{The lmomco functions}
