\name{tlmrnor}
\alias{tlmrnor}
\title{ Compute Select TL-moment ratios of the Normal Distribution }
\description{
This function computes select TL-moment ratios of the Normal distribution for defaults of \eqn{\mu = 0} and \eqn{\sigma = 1}. This function can be useful for plotting the trajectory of the distribution on TL-moment ratio diagrams of \eqn{\tau^{(t_1,t_2)}_2}, \eqn{\tau^{(t_1,t_2)}_3}, \eqn{\tau^{(t_1,t_2)}_4}, \eqn{\tau^{(t_1,t_2)}_5}, and \eqn{\tau^{(t_1,t_2)}_6}. In reality, \eqn{\tau^{(t_1,t_2)}_2} is  dependent on the values for \eqn{\mu} and \eqn{\sigma}.
}
\usage{
tlmrnor(trim=NULL, leftrim=NULL, rightrim=NULL, mu=0, sigma=1)
}
\arguments{
  \item{trim}{Level of symmetrical trimming to use in the computations.
Although \code{NULL} in the argument list, the default is 0---the usual L-moment ratios are returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample.}
  \item{rightrim}{Level of trimming of the right-tail of the sample.}
  \item{mu}{Location parameter (mean) of the distribution.}
  \item{sigma}{Scale parameter (standard deviation) of the distribution.}
}
\value{
  An \R \code{list} is returned.

  \item{tau2}{A vector of the \eqn{\tau^{(t_1,t_2)}_2} values.}
  \item{tau3}{A vector of the \eqn{\tau^{(t_1,t_2)}_3} values.}
  \item{tau4}{A vector of the \eqn{\tau^{(t_1,t_2)}_4} values.}
  \item{tau5}{A vector of the \eqn{\tau^{(t_1,t_2)}_5} values.}
  \item{tau6}{A vector of the \eqn{\tau^{(t_1,t_2)}_6} values.}
}
\note{
The function uses numerical integration of the quantile function of the distribution through the \code{\link{theoTLmoms}} function.
}
\author{ W.H. Asquith }
\seealso{\code{\link{quanor}}, \code{\link{theoTLmoms}}}
\examples{
\dontrun{
tlmrnor(leftrim=2, rightrim=1)
tlmrnor(leftrim=2, rightrim=1, mu=100, sigma=1000) # another slow example
}
}
\keyword{L-moment (trimmed distribution)}
\keyword{Distribution: Normal}
