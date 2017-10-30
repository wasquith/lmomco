\name{supdist}
\alias{supdist}
\title{ The Support of a Distribution based on the Parameters }
\description{
This function takes a parameter object, such as that returned by \code{\link{lmom2par}}, and computes the support (the lower and upper bounds, \eqn{\{L, U\}}) of the distribution given by the parameters. The computation is based on two calls to \code{\link{par2qua}} for the parameters in argument \code{para} (\eqn{\Theta}) and nonexceedance probabilities \eqn{F \in \{0, 1\}}:
\preformatted{lower <- par2qua(0, para)
upper <- par2qua(1, para)}

The quality of \eqn{\{L, U\}} is dependent of the handling of \eqn{F \in \{0,1\}} internal to each quantile function. Across the suite of distributions supported by \pkg{lmomco}, potential applications, and parameter combinations, it difficult to ensure numerical results for the respective \eqn{\{L, U\}} are either very small, are large, or are (or should be) infinite. The distinction is sometimes difficult depending how fast the tail(s) of a distribution is (are) either approaching a limit as \eqn{F} respectively approaches \eqn{0^{+}} or \eqn{1^{-}}.

The intent of this function is to provide a unified portal for \eqn{\{L, U\}} estimation. Most of the time \R (and \pkg{lmomco}) do the right thing anyway and the further overhead within the parameter estimation suite of functions in \pkg{lmomco} is not implemented.

The support returned by this function might be useful in extended application development involving probability density functions \code{pdfCCC} (\eqn{f(x,\Theta)}, see \code{\link{dlmomco}}) and cumulative distribution functions \code{cdfCCC} (\eqn{F(x,\Theta)}, see \code{\link{plmomco}}) functions---both of these functions use as their primary argument a value \eqn{x} that exists along the real number line.
}
\usage{
supdist(para, trapNaN=FALSE, delexp=0.5, paracheck=TRUE, ...)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{trapNaN}{A logical influencing how \code{NaN} are handled (see Note).}
  \item{delexp}{The magnitude of the decrementing of the exponent to search down and up from.  A very long-tailed but highly peaked distribution might require this to be smaller than default.}
    \item{paracheck}{A logical controlling whether the parameters are checked for validity.}
  \item{...}{Additional arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{Three character (minimum) distribution type (for example, \code{type="gev"});}
  \item{support}{The support (or range) of the fitted distribution;}
  \item{nonexceeds}{The nonexceedance probabilities at the computed support.}
  \item{fexpons}{A vector indicating how the respective lower and upper boundaries were arrived at (see Note); and}
  \item{finite}{A logical on each entry of the \code{support} with a preemptive call by the \code{is.finite} function in \R.}
  \item{source}{An attribute identifying the computational source
of the distribution support: \dQuote{supdist}.}
}
\note{
Concerning \code{fexpons}, for the returned vectors of length 2, index 1 is for \eqn{\{L\}} and index 2 is for \eqn{\{U\}}. If an entry in \code{fexpons} is \code{NA}, then \eqn{F = 0} or \eqn{F = 1} for the respective bound was possible. And even if \code{trapNaN} is \code{TRUE}, no further refinement on the bounds was attempted.

On the otherhand, if \code{trapNaN} is \code{TRUE} and if the bounds \eqn{\{L\}} and (or) \eqn{\{U\}} is not \code{NA}, then an attempt was made to move away from \eqn{F \in \{0,1\}} in incremental integer exponent from \eqn{0^{+}} or \eqn{1^{-}} until a \code{NaN} was not encountered. The integer exponents are \eqn{i \in [-(\phi), -(\phi - 1), \ldots, -4]}, where \eqn{\phi} \code{= .Machine$sizeof.longdouble} and \eqn{-4} is a hardwired limit (1 part in 10,000). In the last example in the Examples section, the \eqn{\{U\}} for \eqn{F=1} quantile is \code{NaN} but \eqn{1 - 10^i} for which \eqn{i = -16}, which also is the \code{.Machine$sizeof.longdouble} on the author's development platform.

 At first release, it seems there was justification in triggering this to \code{TRUE} if a quantile function returns a \code{NA} when asked for \eqn{F = 0} or \eqn{F = 1}---some quantile functions partially trapped \code{NaN}s themselves. So even if \code{trapNaN == FALSE}, it is triggered to \code{TRUE} if a \code{NA} is discovered as described. \emph{Users are encouraged to discuss adaptions or changes to the implementation of \code{supdist} with the author.}

Thus it should be considered a feature of \code{supdist} that should a quantile function already trap errors at either \eqn{F = 0} or \eqn{F = 1} and return \code{NA}, then \code{trapNaN} is internally set to \code{TRUE} regardless of being originally \code{FALSE} and the preliminary limit is reset to \code{NaN}. The Rice distribution \code{\link{quarice}} is one such example that internally already traps an \eqn{F = 1} by returning \eqn{x(F{=}1) = }\code{NA}.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lmom2par}} }
\examples{
lmr <- lmoms(c(33, 37, 41, 54, 78, 91, 100, 120, 124))
supdist(lmom2par(lmr, type="gov" )) # Lower = 27.41782, Upper = 133.01470
supdist(lmom2par(lmr, type="gev" )) # Lower = -Inf,     Upper = 264.4127

supdist(lmom2par(lmr, type="wak" ))               # Lower = 16.43722, Upper = NaN
supdist(lmom2par(lmr, type="wak" ), trapNaN=TRUE) # Lower = 16.43722, Upper = 152.75126
#$support  16.43722  152.75126
#$fexpons        NA  -16
#$finite       TRUE  TRUE
\dontrun{
para <- vec2par(c(0.69, 0.625), type="kmu") # very flat tails and narrow peak!
supdist(para, delexp=1   )$support # [1] 0        NaN
supdist(para, delexp=0.5 )$support # [1] 0.000000 3.030334
supdist(para, delexp=0.05)$support # [1] 0.000000 3.155655
# This distribution appears to have a limit at PI and the delexp=0.5
}
}
\keyword{utility (distribution)}
\keyword{distribution (utility)}
\concept{distribution bounds}
\concept{distribution bound}
\concept{bounds}
\concept{distribution support}
\concept{distribution limits}
\concept{distribution limit}
\concept{tail limits}
\concept{tail limit}

