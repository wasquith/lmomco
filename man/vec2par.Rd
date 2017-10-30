\name{vec2par}
\alias{vec2par}
\title{Convert a Vector of Parameters to a Parameter Object of a Distribution}
\description{
This function converts a vector of parameters to a parameter object of a distribution. The type of distribution is specified in the
argument list:
\code{aep4}, \code{cau}, \code{exp}, \code{gam}, \code{gep}, \code{gev},
\code{glo}, \code{gno}, \code{gpa}, \code{gum}, \code{kap},
\code{kur}, \code{lap}, \code{lmrq}, \code{ln3}, \code{nor},
\code{pe3}, \code{ray}, \code{revgum}, \code{rice}, \code{st3},
\code{texp}, \code{wak}, and \code{wei}. These abbreviations and only these are used in routing logic within \pkg{lmomco}. There is no provision for fuzzy matching. However, if the distribution type is not identified, then the function issues a warning, but goes ahead and creates the parameter list and of course can not check for the validity of the parameters. If one has a need to determine on-the-fly the number of parameters in a distribution as supported in \pkg{lmomco}, then see the \code{\link{dist.list}} function.
}
\usage{
vec2par(vec, type, nowarn=FALSE, paracheck=TRUE, ...)
}
\arguments{
  \item{vec}{A vector of parameter values for the distribution specified by type.}
  \item{type}{Three character distribution type (for example, \code{type='gev'}).}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity. Overriding of this check might be extremely important and needed for use of the distribution quantile function in the context of TL-moments with nonzero trimming.}
  \item{...}{Additional arguments for the \code{\link{are.par.valid}} call that is made internally.}
}
\value{
  An \R \code{list} is returned. This list should contain at least the following items, but some distributions such as the \code{revgum} have extra.

  \item{type}{The type of distribution in three character format.}
  \item{para}{The parameters of the distribution.}
  \item{source}{Attribute specifying source of the parameters---\dQuote{vec2par}.}
}
\details{
If the distribution is a Reverse Gumbel (\code{type=revgum}) or Generalized Pareto (\code{type=gpa}), which are 2-parameter or 3-parameter distributions, the third or fourth value in the vector is the \eqn{\zeta} of the distribution. \eqn{\zeta} represents the fraction of the sample that is noncensored, or number of observed (noncensored) values divided by the sample size. The \eqn{\zeta} represents censoring on the right, that is there are unknown observations above a threshold or the largest observed sample. Consultation of \code{\link{parrevgum}} or \code{\link{pargpaRC}} should elucidate the censoring discussion.
}
\note{
If the \code{type} is not amongst the official list given above, then the \code{type} given is loaded into the \code{type} element of the returned list and an other element \code{isuser = TRUE} is also added. There is no \code{isuser} created if the distribution is supported by \pkg{lmomco}. This is an attempt to given some level of flexibility so that others can create their own distributions or conduct research on derivative code from \pkg{lmomco}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmom2par}}, \code{\link{par2vec}}}
\examples{
para <- vec2par(c(12,123,0.5),'gev')
Q <- quagev(0.5,para)

my.custom <- vec2par(c(2,2), type='myowndist') # Think about making your own
}
\keyword{utility (vector)}
