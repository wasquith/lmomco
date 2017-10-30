\name{rlmomco}
\alias{rlmomco}
\title{Random Variates of a Distribution}
\description{
This function generates random variates for the specified distribution in the parameter object argument. See documentation about the parameter object is seen in \code{\link{lmom2par}} or \code{\link{vec2par}}. The prepended \code{r} in the function name is to parallel the built-in distribution syntax of \R but of course reflects the \pkg{lmomco} name in the function. An assumption is made that the user knows that they are providing appropirate (that is valid) distribution parameters. This is evident by the \preformatted{paracheck = FALSE} argument passed to the \code{\link{par2qua}} function.
}
\usage{
rlmomco(n, para)
}
\arguments{
  \item{n}{Number of samples to generate}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
}
\value{
  Vector of quantile values.
}
\note{
The action of this function in \R idiom is \code{par2qua(runif(n), para)} for the distribution parameters \code{para}, the \R function \code{runif} is the Uniform distribution, and \code{n} being the simulation size. 
}
\author{ W.H. Asquith}
\seealso{\code{\link{dlmomco}}, \code{\link{plmomco}}, \code{\link{qlmomco}}, \code{\link{slmomco}} }
\examples{
lmr      <- lmoms(rnorm(20)) # generate 20 standard normal variates
para     <- parnor(lmr) # estimate parameters of the normal
simulate <- rlmomco(20,para) # simulate 20 samples using lmomco package

lmr  <- vec2lmom(c(1000,500,.3)) # first three lmoments are known
para <- lmom2par(lmr,type="gev") # est. parameters of GEV distribution
Q    <- rlmomco(45,para) # simulate 45 samples
PP   <- pp(Q)            # compute the plotting positions
plot(PP,sort(Q))         # plot the data up
}
\keyword{random variate}
\keyword{quantile function}
\concept{random variate}
\keyword{The lmomco functions}
\keyword{The lmomco function mimics of R nomenclature}
