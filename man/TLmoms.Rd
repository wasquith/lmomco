\name{TLmoms}
\alias{TLmoms}
\title{ The Sample Trimmed L-moments and L-moment Ratios }
\description{
Compute the sample trimmed L-moments (TL-moments) for a vector. The level of symmetrical trimming is specified. The mathematical expression for a TL-moment is seen under \code{\link{TLmom}}.  The \code{\link{TLmoms}} function loops across that expression and the \code{\link{TLmom}} function for each \code{nmom}=\eqn{r} set in the argument list.
}
\usage{
TLmoms(x, nmom, trim=NULL, leftrim=NULL, rightrim=NULL, vecit=FALSE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{trim}{Level of symmetrical trimming to use in the computations. Although \code{NULL} is in the argument list, the default is 0---the usual L-moment is returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
  \item{rightrim}{Level of trimming of the right-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
    \item{vecit}{A logical to return the first two \eqn{\lambda_i \in 1,2} and then the \eqn{\tau_i \in 3,\cdots} where the length of the returned vector is controlled by the \code{nmom} argument. This argument will store the trims in the attributes of the returned vector, but caution is advised if \code{\link{vec2par}} were to be used on the vector because that function does not consult the trimming.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the TL-moments. First element is
\eqn{\hat{\lambda}^{(t_1,t_2)}_1}, second element is \eqn{\hat{\lambda}^{(t_1,t_2)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\hat{\tau}^{(t_1,t_2)}}, third element is \eqn{\hat{\tau}^{(t_1,t_2)}_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which will equal \code{trim} if symmetrical trimming was used.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which will equal \code{trim} if symmetrical trimming was used.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{TLmoms}.}
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational statistics and data analysis, v. 43, pp. 299-314.
}
\author{ W.H. Asquith}

\seealso{\code{\link{TLmom}}, \code{\link{lmoms}}, and \code{\link{lmorph}}}
\examples{
X1 <- rcauchy(30)
TL <- TLmoms(X1,nmom=6,trim=1)

# This trimming with remove the 1 and the two 4s. All values passed on to the TLmom()
# function then are equal and number of L-moments is too big as well. TLmom() returns
# NaN but these are intercepted and systematically changed to NAs.
TLmoms(c(1,2,2,2,4,4), leftrim=1, rightrim=2, nmom=6)$lambdas
# [1]  2  0  0 NA NA NA

# Example of zero skewness (Berry Boessenkool)
TLmoms(c(3.2, 4.4, 4.8, 2.6, 3.6))
}
\keyword{L-moment (sample)}
\keyword{L-moment (trimmed sample)}
