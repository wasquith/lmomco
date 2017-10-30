\name{lmomsRCmark}
\alias{lmomsRCmark}
\title{ Sample L-moments Moments for Right-Tail Censoring by a Marking Variable }
\description{
Compute the sample L-moments for right-tail censored data set in which censored data values are identified by a marking variable. Extension of left-tail censoring can be made using \code{\link{fliplmoms}} and the example therein.
}
\usage{
lmomsRCmark(x, rcmark=NULL, nmom=5, flip=NA, flipfactor=1.1)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{rcmark}{The right-tail censoring (upper) marking variable for unknown threshold: 0 is uncensored, 1 is censored.}
  \item{nmom}{Number of L-moments to return.}
  \item{flip}{Do the data require flipping so that left-censored data can be processed as such. If the flip is a logical and \code{TRUE}, then \code{flipfactor} \eqn{\times} \eqn{\mathrm{max}(x)} (the maximum of \code{x}) is used. If the \code{flip} is a numeric, then it is used as the flip.}
  \item{flipfactor}{The value that is greater than 1, which is multiplied on the maximum of \code{x} to determine the flip, if the \code{flip} is not otherwise provided.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\hat{\lambda}^{(0,0)}_1}, second element is \eqn{\hat{\lambda}^{(0,0)}_2}, and so on.  \emph{The returned mean is NOT unflipped.}}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\hat{\tau}^{(0,0)}}, third element is \eqn{\hat{\tau}^{(0,0)}_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal \code{NULL} if asymmetrical trimming was used. This is not currently implemented as no one has done the derivations.}
  \item{leftrim}{Level of left-tail trimming used in the computation. This is not currently implemented as no one has done the derivations.}
  \item{rightrim}{Level of right-tail trimming used in the computation. This is not currently implemented as no one has done the derivations.}
  \item{n}{The complete sample size.}
  \item{n.cen}{The number of right-censored data values.}
  \item{flip}{The flip used in the computations for support of left-tail censoring.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomsRCmark}.}
}
\references{
Wang, Dongliang, Hutson, A.D., Miecznikowski, J.C., 2010, L-moment estimation for parametric survival models given censored data: Statistical Methodology, v. 7, no. 6, pp. 655--667.

Helsel, D.R., 2005, Nondetects and data analysis---Statistics for censored environmental data: Hoboken, New Jersey, John Wiley, 250 p.
}
\author{W.H. Asquith}
\seealso{ \code{\link{lmomRCmark}}, \code{\link{fliplmoms}} }
\examples{
# Efron, B., 1988, Logistic regression, survival analysis, and the
# Kaplan-Meier curve: Journal of the American Statistical Association,
# v. 83, no. 402, pp.414--425
# Survival time measured in days for 51 patients with a marking
# variable in the "time,mark" ensemble. If marking variable is 1,
# then the time is right-censored by an unknown censoring threshold.
Efron <-
c(7,0,  34,0,  42,0,  63,0,  64,0,  74,1,  83,0,  84,0,  91,0,
108,0,  112,0,  129,0,  133,0,  133,0,  139,0,  140,0,  140,0,
146,0,  149,0,  154,0,  157,0,  160,0,  160,0,  165,0,  173,0,
176,0,  185,1,  218,0,  225,0,  241,0,  248,0,  273,0,  277,0,
279,1,  297,0,  319,1,  405,0,  417,0,  420,0,  440,0,  523,1,
523,0,  583,0,  594,0,  1101,0,  1116,1,  1146,0,  1226,1,
1349,1,  1412,1, 1417,1);

# Break up the ensembles into to vectors
ix <- seq(1,length(Efron),by=2)
T  <- Efron[ix]
Efron.data <- T;
Efron.rcmark <- Efron[(ix+1)]

lmr.RC <- lmomsRCmark(Efron.data, rcmark=Efron.rcmark)
lmr.ub <- lmoms(Efron.data)
lmr.noRC <- lmomsRCmark(Efron.data)
PP <- pp(Efron.data)
plot(PP, Efron.data, col=(Efron.rcmark+1), ylab="DATA")
lines(PP, qlmomco(PP, lmom2par(lmr.noRC, type="kap")), lwd=3, col=8)
lines(PP, qlmomco(PP, lmom2par(lmr.ub, type="kap")))
lines(PP, qlmomco(PP, lmom2par(lmr.RC, type="kap")), lwd=2, col=2)
legend(0,1000,c("uncensored L-moments by indicator (Kappa distribution)",
                "unbiased L-moments (Kappa)",
           "right-censored L-moments by indicator (Kappa distribution)"),
                lwd=c(3,1,2), col=c(8,1,2))

########
ZF <- 5 # discharge of undetection of streamflow
Q <- c(rep(ZF,8), 116, 34, 56, 78, 909, 12, 56, 45, 560, 300, 2500)
Qc <- Q == ZF; Qc <- as.numeric(Qc)
lmr     <- lmoms(Q)
lmr.cen <- lmomsRCmark(Q, rcmark=Qc, flip=TRUE)
flip <- lmr.cen$flip
fit  <- pargev(lmr);                     fit.cen <- pargev(lmr.cen)
F <- seq(0.001, 0.999, by=0.001)
Qfit     <-        qlmomco(    F, fit    )
Qfit.cen <- flip - qlmomco(1 - F, fit.cen) # remember to reverse qdf
plot(pp(Q),sort(Q), log="y", xlab="NONEXCEED PROB.", ylab="QUANTILE")
lines(F, Qfit);   lines(F, Qfit.cen,col=2)
}
\keyword{L-moment (sample)}
\keyword{data censoring}
