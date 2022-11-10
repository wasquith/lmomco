\name{parTLgpa}
\alias{parTLgpa}
\title{Estimate the Parameters of the Generalized Pareto Distribution using Trimmed L-moments}
\description{
This function estimates the parameters of the Generalized Pareto distribution given
the the trimmed L-moments (TL-moments) for \eqn{t=1} of the data in TL-moment object with a trim level of unity (\code{trim=1}). The parameters are computed as
\deqn{\kappa = \frac{10-45\tau^{(1)}_3}{9\tau^{(1)}_3+10} \mbox{,}}
\deqn{\alpha = \frac{1}{6}\lambda^{(1)}_2(\kappa+2)(\kappa+3)(\kappa+4) \mbox{, and}}
\deqn{\xi = \lambda^{(1)}_1 - \frac{\alpha(\kappa+5)}{(\kappa+2)(\kappa+3)} \mbox{.}}
}
\usage{
parTLgpa(lmom, ...)
}
\arguments{
  \item{lmom}{A TL-moment object created by \code{\link{TLmoms}}.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{gpa}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parTLgpa}.}
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{W.H. Asquith}
\seealso{\code{\link{TLmoms}}, \code{\link{lmomTLgpa}},
\code{\link{cdfgpa}}, \code{\link{pdfgpa}}, \code{\link{quagpa}}
}
\examples{
TL <- TLmoms(rnorm(20),trim=1)
parTLgpa(TL)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Generalized Pareto}
