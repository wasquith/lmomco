\name{TLmom}
\alias{TLmom}
\title{A Sample Trimmed L-moment }
\description{
A sample trimmed L-moment (TL-moment) is computed for a vector. The \eqn{r \ge 1} order of the L-moment is specified as well as the level of symmetrical trimming.  A trimmed TL-moment
\eqn{\hat{\lambda}^{(t_1,t_2)}_r} is
\deqn{ \hat{\lambda}^{(t_1,t_2)}_r = \frac{1}{r}\sum^{n-t_2}_{i=t_1+1}
           \left[ \frac{\sum\limits^{r-1}_{k=0}{ (-1)^k {r-1 \choose k}
                                                 {i-1 \choose r+t_1-1-k}
                                                 {n-i \choose t_2+k}
                                         }}{{n \choose r+t_1+t_2}}
           \right] x_{i:n} \mbox{,}}
where \eqn{t_a} represents the trimming level of the \eqn{t_2}-largest or \eqn{t_1}-smallest values, \eqn{r} represents the order of the L-moment, \eqn{n} represents the sample size, and \eqn{x_{i:n}} represents the \eqn{i}th sample order statistic (\eqn{x_{1:n} \le x_{2:n} \le \dots \le x_{n:n}}).
}
\usage{
TLmom(x, order, trim=NULL, leftrim=NULL, rightrim=NULL, sortdata=TRUE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{order}{L-moment order to use in the computations.
Default is 1 (the mean).}
  \item{trim}{Level of symmetrical trimming to use in the computations. Although \code{NULL} is in the argument list, the default is 0---the usual L-moment is returned.}
  \item{leftrim}{Level of trimming of the left-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
  \item{rightrim}{Level of trimming of the right-tail of the sample, which should be left to \code{NULL} if no or symmetrical trimming is used.}
  \item{sortdata}{A logical switch on whether the data should be sorted. The default is TRUE.}
}
\value{
  An \R \code{list} is returned.

  \item{lambda}{The TL-moment of order=\code{order}, \eqn{\hat{\lambda}^{(t_1,t_2)}_r} where \eqn{r} is the moment order, \eqn{t_1} is left-tail trimming, and \eqn{t_2} is right-tail trimming.}
  \item{order}{L-moment order computed. Default is 1 (the mean).}
  \item{trim}{Level of symmetrical trimming used in the computation.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which will equal \code{trim} if symmetrical trimming was used.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which will equal \code{trim} if symmetrical trimming was used.}
}
\note{
The presence of the \code{sortdata} switch can be dangerous. L-moment computation requires that the data be sorted into the \dQuote{order statistics}. Thus the default behavior of \code{sortdata=TRUE} is required when the function is called on its own. In practice, this function would almost certainly not be used on its own because multiple trimmed L-moments would be needed. Multiple trimmed L-moments are best computed by \code{\link{TLmoms}}, which calls \code{\link{TLmom}} multiple times. The function \code{\link{TLmoms}} takes over the sort operation on the data and passes \code{sortdata=FALSE} to \code{\link{TLmom}} for efficiency. (The point of this discussion is that CPU time is not wasted sorting the data more than once.)
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics and Data Analysis, v. 43, pp. 299--314.
}
\author{W.H. Asquith}
\seealso{\code{\link{TLmoms}}}
\examples{
X1 <- rcauchy(30)
TL <- TLmom(X1,order=2,trim=1)
}
\keyword{L-moment (sample)}
\keyword{L-moment (trimmed sample)}
