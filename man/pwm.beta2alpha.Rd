\name{pwm.beta2alpha}
\alias{pwm.beta2alpha}
\alias{pwm.alpha2beta}
\title{Conversion of Beta to Alpha Probability-Weighted Moments (PWMs) or Alpha to Beta PWMs }
\description{
Conversion of \dQuote{beta} (the well known ones) to \dQuote{alpha} probability-weighted moments (PWMs) by \code{\link{pwm.beta2alpha}} or alpha to beta PWMs by \code{\link{pwm.alpha2beta}}. The relations between the \eqn{\alpha} and \eqn{\beta} PWMs are
\deqn{\alpha_r = \sum^r_{k=0} (-1)^k {r \choose k} \beta_k\mbox{,}}
and
\deqn{\beta_r = \sum^r_{k=0} (-1)^k {r \choose k} \alpha_k\mbox{.}}
Lastly, note that the \eqn{\beta} are almost exclusively used in the literature. Because each is a linear combination of the other, they are equivalent in meaning but not numerically.
}
\usage{
pwm.beta2alpha(pwm)

pwm.alpha2beta(pwm)
}
\arguments{
  \item{pwm}{A vector of alpha or beta probability-weighted moments depending on which related function is called.}
}
\value{
  If \eqn{\beta_r \rightarrow \alpha_r} (\code{\link{pwm.beta2alpha}}), a vector of the \eqn{\alpha_r}. Note that convention is the have a \eqn{\alpha_0}, but this is placed in the first index \code{i=1} vector. Alternatively, if \eqn{\alpha_r \rightarrow \beta_r} (\code{\link{pwm.alpha2beta}}), a vector of the \eqn{\beta_r}.
}
\references{
# NEED
}
\author{W.H. Asquith}
\seealso{\code{\link{pwm}}, \code{\link{pwm2lmom}} }
\examples{
X <- rnorm(100)
pwm(X)$betas
pwm.beta2alpha(pwm(X)$betas)
pwm.alpha2beta(pwm.beta2alpha(pwm(X)$betas))
}
\keyword{probability-weighted moment (theoretical)}
