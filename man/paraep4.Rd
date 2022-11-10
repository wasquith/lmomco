\name{paraep4}
\alias{paraep4}
\title{Estimate the Parameters of the 4-Parameter Asymmetric Exponential Power Distribution}
\description{
This function estimates the parameters of the 4-parameter Asymmetric Exponential Power distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relation between distribution parameters and L-moments is seen under \code{\link{lmomaep4}}. Relatively straightforward, but difficult to numerically achieve, optimization is needed to extract the parameters from the L-moments.

Delicado and Goria (2008) argue for numerical methods to use the following objective function
\deqn{\epsilon(\alpha, \kappa, h) = \log(1 + \sum_{r=2}^4 (\hat\lambda_r - \lambda_r)^2)\mbox{,}}
and subsequently solve directly for \eqn{\xi}. This objective function was chosen by Delicado and Goria because the solution surface can become quite flat for away from the minimum.  The author of \pkg{lmomco} agrees with the findings of those authors from limited exploratory analysis and the development of the algorithms used here under the rubic of the \dQuote{DG} method. This exploration resulted in an alternative algorithm using tabulated initial guesses described below. An evident drawback of the Delicado-Goria algorithm, is that precision in \eqn{\alpha} is may be lost according to the observation that this parameter can be analytically computed given \eqn{\lambda_2}, \eqn{\kappa}, and \eqn{h}.

It is established practice in L-moment theory of four (and similarly three) parameter distributions to see expressions for \eqn{\tau_3} and \eqn{\tau_4} used for numerical optimization to obtain the two higher parameters (\eqn{\alpha} and \eqn{h}) first and then see analytical expressions directly compute the two lower parameters (\eqn{\xi} and \eqn{\alpha}). The author made various exploratory studies by optimizing on \eqn{\tau_3} and \eqn{\tau_4} through a least squares objective function. Such a practice seems to perform acceptably when compared to that recommended by Delicado and Goria (2008) when the initial guesses for the parameters are drawn from pretabulation of the relation between \eqn{\{\alpha, h\}} and \eqn{\{\tau_3, \tau_4\}}.

Another optimization, referred to here as the \dQuote{A} (Asquith) method, is available for parameter estimation using the following objective function
\deqn{\epsilon(\kappa, h) = \sqrt{(\hat\tau_3 - \tau_3)^2 + (\hat\tau_4 - \tau_4)^2}\mbox{,}}
and subsequently solve directly for \eqn{\alpha} and then \eqn{\xi}. The \dQuote{A} method appears to perform better in \eqn{\kappa} and \eqn{h} estimation and quite a bit better in \eqn{\alpha} and and \eqn{\xi} as seemingly expected because these last two are analytically computed (Asquith, 2014). The objective function of the \dQuote{A} method defaults to use of the \eqn{\sqrt{x}} but this can be removed by setting \code{sqrt.t3t4=FALSE}.

The initial guesses for the \eqn{\kappa} and \eqn{h} parameters derives from a hashed environment in in file \cr \file{sysdata.rda} (\env{.lmomcohash$AEPkh2lmrTable}) in which the \eqn{\{\kappa, h\}} pair having the minimum \eqn{\epsilon(\kappa, h)} in which \eqn{\tau_3} and \eqn{\tau_4} derive from the table as well.  The file \file{SysDataBuilder.R} provides additional technical details on how the \code{AEPkh2lmrTable} was generated. The table represents a systematic double-loop sweep through \code{\link{lmomaep4}} for
\deqn{\kappa \mapsto \{-3 \le \log(\kappa) \le 3, \Delta\log(\kappa)=0.05\}\mbox{,}}
and
\deqn{h \mapsto \{-3 \le \log(h) \le 3, \Delta\log(h)=0.05\}\mbox{.}}

The function will not return parameters if the following lower (estimated) bounds of \eqn{\tau_4} are not met: \cr \eqn{\tau_4 \ge 0.77555(|\tau_3|) - 3.3355(|\tau_3|)^2 + 14.196(|\tau_3|)^3 - 29.909(|\tau_3|)^4 + 37.214(|\tau_3|)^5 - 24.741(|\tau_3|)^6 + 6.7998(|\tau_3|)^7}. For this polynomial, the residual standard error is RSE = 0.0003125 and the maximum absolute error for \eqn{\tau_3{:}[0,1] < 0.0015}. The actual coefficients in \code{\link{paraep4}} have additional significant figures. However, the argument \code{snap.tau4}, if set, will set \eqn{\tau_4} equal to the prediction from the polynomial. This value of \eqn{\tau_4} should be close enough numerically to the boundary because the optimization is made using a log-transformation to ensure that \eqn{\alpha}, \eqn{\kappa}, and \eqn{h} remain in the positive domain---though the argument \code{nudge.tau4} is provided to offset \eqn{\tau_4} upward just incase of optimization problems.
}
\usage{
paraep4(lmom, checklmom=TRUE, method=c("A", "DG", "ADG"),
        sqrt.t3t4=TRUE, eps=1e-4, checkbounds=TRUE, kapapproved=TRUE,
        snap.tau4=FALSE, nudge.tau4=0,
        A.guess=NULL, K.guess=NULL, H.guess=NULL, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the L-moments be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{method}{Which method for parameter estimation should be used. The \dQuote{A} or \dQuote{DG} methods. The \dQuote{ADG} method will run both methods and retains the salient optimization results of each but the official parameters in \code{para} are those from the \dQuote{A} method. Lastly, all minimization is based on the \code{optim} function using the Nelder--Mead method and default arguments.}
  \item{sqrt.t3t4}{If true and the method is \dQuote{A}, then the square root of the sum of square errors in \eqn{\tau_3} and \eqn{\tau_4} are used instead of sum of square differences alone.}
 \item{eps}{A small term or threshold for which the square root of the sum of square errors in \eqn{\tau_3} and \eqn{\tau_4} is compared to to judge \dQuote{good enough} for the alogrithm to set the \code{ifail} on return in addition to convergence flags coming from the \code{optim} function. Note that \code{eps} is only used if the \dQuote{A} or \dQuote{ADG} methods are triggered because the other method uses the scale parameter which in reality could be quite large relative to the other two shape parameters, and a reasonable default for such a secondary error threshold check would be ambiguous.}
  \item{checkbounds}{Should the lower bounds of \eqn{\tau_4} be verified and if sample \eqn{\hat\tau_3} and \eqn{\hat\tau_4} are outside of these bounds, then \code{NA} are returned for the solutions.}
  \item{kapapproved}{Should the Kappa distribution be fit by \code{\link{parkap}} if \eqn{\hat\tau_4} is below the lower bounds of \eqn{\tau_4}? This fitting is only possible if \code{checkbounds} is true. The Kappa and AEP4 overlap partially. The AEP4 extends \eqn{\tau_4} above Generalized Logistic and Kappa extends \eqn{\tau_4} below the lower bounds of \eqn{\tau_4} for AEP4 and extends all the way to the theoretical limits as used within \code{\link{are.lmom.valid}}.}
  \item{snap.tau4}{A logical to \dQuote{snap} the \eqn{\tau_4} upwards to the lower boundary if the given \eqn{\tau_4} is lower than the boundary described in the polynomial.}
  \item{nudge.tau4}{An offset to the snapping of \eqn{\tau_4} intended to move \eqn{\tau_4} just above the lower bounds in case of optimization problems. (The absolute value of the nudge is made internally to ensure only upward adjustment by an addition operation.)}
  \item{A.guess}{A user specified guess of the \eqn{\alpha} parameter to provide to the optimization of any of the methods. This argument just superceeds the simple initial guess of \eqn{\alpha = 1}.}
  \item{K.guess}{A user specified guess of the \eqn{\kappa} parameter to supercede that derived from the \env{.lmomcohash$AEPkh2lmrTable} in file \file{sysdata.rda}.}
  \item{H.guess}{A user specified guess of the \eqn{h} parameter to supercede that derived from the \env{.lmomcohash$AEPkh2lmrTable} in file \file{sysdata.rda}.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{aep4}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{paraep4}.}
  \item{method}{The method as specified by the \code{method}.}
  \item{ifail}{A numeric failure code.}
  \item{ifailtext}{A text message for the failure code.}
  \item{L234}{Optional and dependent on method \dQuote{DG} or \dQuote{ADG}. Another \R \code{list} containing the optimization details by the \dQuote{DG} method along with the estimated parameters in \code{para_L234}. The \dQuote{_234} is to signify that optimization is made using \eqn{\lambda_2}, \eqn{\lambda_3}, and \eqn{\lambda_4}.  The parameter values in \code{para} are those only when the \dQuote{DG} method is used.}
  \item{T34}{Optional and dependent on method \dQuote{A} or \dQuote{ADG}. Another \R \code{list} containing the optimization details by the \dQuote{A} method along with the estimated parameters in \code{para_T34}. The \dQuote{_T34} is to signify that opimization is being conducted using \eqn{\tau_3} and \eqn{\tau_4} only. The parameter values in \code{para} are those by the \dQuote{A} method.}

The values for \code{ifail} or produced by three mechanisms. First, the convergence number emanating from the \code{optim} function itself. Second, the integer 1 is used when the failure is attributable to the \code{optim} function. Third, the interger 2 is a general attempt to have a singular failure by sometype of \code{eps} outside of \code{optim}. Fourth, the integer 3 is used to show that the parameters fail against a parameter validity check in \code{\link{are.paraep4.valid}}. And fifth, the integer 4 is used to show that the sample L-moments are below the lower bounds of the \eqn{\tau_4} polynomial shown here.

Additional and self explanatory elements on the returned list will be present if the Kappa distribution was fit instead.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.

Delicado, P., and Goria, M.N., 2008, A small sample comparison of maximum likelihood,
moments and L-moments methods for the asymmetric exponential power distribution:
Computational Statistics and Data Analysis, v. 52, no. 3, pp. 1661--1673.
}
\seealso{\code{\link{lmomaep4}}, \code{\link{cdfaep4}}, \code{\link{pdfaep4}}, \code{\link{quaaep4}}, \code{\link{quaaep4kapmix}}
}
\examples{
# As a general rule AEP4 optimization can be CPU intensive

\dontrun{
lmr <- vec2lmom(c(305, 263, 0.815, 0.631))
plotlmrdia(lmrdia()); points(lmr$ratios[3], lmr$ratios[4], pch=16, cex=3)
PAR <- paraep4(lmr, snap.tau4=TRUE) # will just miss the default eps
FF <- nonexceeds(sig6=TRUE)
plot(FF, quaaep4(FF, PAR), type="l", log="y")
lmomaep4(PAR) # 305, 263, 0.8150952, 0.6602706 (compare to those in lmr) }

\dontrun{
PAR <- list(para=c(100, 1000, 1.7, 1.4), type="aep4")
lmr <- lmomaep4(PAR)
aep4 <- paraep4(lmr, method="ADG")
print(aep4) # }

\dontrun{
PARdg  <- paraep4(lmr, method="DG")
PARasq <- paraep4(lmr, method="A")
print(PARdg)
print(PARasq)
F <- c(0.001, 0.005, seq(0.01,0.99, by=0.01), 0.995, 0.999)
qF <- qnorm(F)
ylim <- range( quaaep4(F, PAR), quaaep4(F, PARdg), quaaep4(F, PARasq) )
plot(qF, quaaep4(F, PARdg), type="n", ylim=ylim,
     xlab="STANDARD NORMAL VARIATE", ylab="QUANTILE")
lines(qF, quaaep4(F, PAR), col=8, lwd=10) # the true curve
lines(qF, quaaep4(F, PARdg),  col=2, lwd=3)
lines(qF, quaaep4(F, PARasq), col=3, lwd=2, lty=2)
# See how the red curve deviates, Delicado and Goria failed
# and the ifail attribute in PARdg is TRUE. Note for lmomco 2.3.1+
# that after movement to log-exp transform to the parameters during
# optimization that this "error" as described does not appear to occur.

print(PAR$para)
print(PARdg$para)
print(PARasq$para)

ePAR1dg <- abs((PAR$para[1] - PARdg$para[1])/PAR$para[1])
ePAR2dg <- abs((PAR$para[2] - PARdg$para[2])/PAR$para[2])
ePAR3dg <- abs((PAR$para[3] - PARdg$para[3])/PAR$para[3])
ePAR4dg <- abs((PAR$para[4] - PARdg$para[4])/PAR$para[4])

ePAR1asq <- abs((PAR$para[1] - PARasq$para[1])/PAR$para[1])
ePAR2asq <- abs((PAR$para[2] - PARasq$para[2])/PAR$para[2])
ePAR3asq <- abs((PAR$para[3] - PARasq$para[3])/PAR$para[3])
ePAR4asq <- abs((PAR$para[4] - PARasq$para[4])/PAR$para[4])

MADdg  <- mean(ePAR1dg,  ePAR2dg,  ePAR3dg,  ePAR4dg)
MADasq <- mean(ePAR1asq, ePAR2asq, ePAR3asq, ePAR4asq)

# We see that the Asquith method performs better for the example
# parameters in PAR and inspection of the graphic will show that
# the Delicado and Goria solution is obviously off. (See Note above)
print(MADdg)
print(MADasq)

# Repeat the above with this change in parameter to
# PAR <- list(para=c(100, 1000, .7, 1.4), type="aep4")
# and the user will see that all three methods converged on the
# correct values. }
}
\keyword{distribution (parameters)}
\keyword{Distribution: Asymmetric Exponential Power}
