\name{tlmrgum}
\alias{tlmrgum}
\title{ Compute Select TL-moment ratios of the Gumbel Distribution }
\description{
This function computes select TL-moment ratios of the Gumbel distribution for defaults of \eqn{\xi = 0} and \eqn{\alpha = 1}. This function can be useful for plotting the trajectory of the distribution on TL-moment ratio diagrams of \eqn{\tau^{(t_1,t_2)}_2}, \eqn{\tau^{(t_1,t_2)}_3}, \eqn{\tau^{(t_1,t_2)}_4}, \eqn{\tau^{(t_1,t_2)}_5}, and \eqn{\tau^{(t_1,t_2)}_6}. In reality, \eqn{\tau^{(t_1,t_2)}_2} is  dependent on the values for \eqn{\xi} and \eqn{\alpha}.
}
\usage{
tlmrgum(trim=NULL, leftrim=NULL, rightrim=NULL, xi=0, alpha=1)
}
\arguments{
  \item{trim}{Level of symmetrical trimming to use in the computations.
Although \code{NULL} in the argument list, the default is 0---the usual L-moment ratios are returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample.}
  \item{rightrim}{Level of trimming of the right-tail of the sample.}
  \item{xi}{Location parameter of the distribution.}
  \item{alpha}{Scale parameter of the distribution.}
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
\seealso{\code{\link{quagum}}, \code{\link{theoTLmoms}}}
\examples{
\dontrun{
tlmrgum(trim=2)
tlmrgum(trim=2, xi=2) # another slow example
}
}
\keyword{L-moment (trimmed distribution)}
\keyword{Distribution: Gumbel}
