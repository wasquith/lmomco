\name{parcau}
\alias{parcau}
\title{Estimate the Parameters of the Cauchy Distribution}
\description{
This function estimates the parameters of the Cauchy distribution from the trimmed L-moments (TL-moments) having trim level 1. The relations between distribution parameters and the TL-moments (trim=1) are seen under \code{\link{lmomcau}}.
}
\usage{
parcau(lmom, ...)
}
\arguments{
  \item{lmom}{A TL-moment object from \code{\link{TLmoms}} with \code{trim=1}.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{cau}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parcau}.}
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{ W.H. Asquith}
\seealso{\code{\link{TLmoms}}, \code{\link{lmomcau}}, \code{\link{cdfcau}}, \code{\link{pdfcau}}, \code{\link{quacau}}
}
\examples{
X1 <- rcauchy(20)
parcau(TLmoms(X1,trim=1))
}
\keyword{distribution (parameters)}
\keyword{Distribution: Cauchy}
