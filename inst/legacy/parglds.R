"parTLgld" <-
function(lmom,result='best',verbose=FALSE,extract=0,initkh=NULL) {
    if(is.null(lmom$trim)) {
      warning("An L-moment object and not a TL-moment object appears to have been passed")
      return()
    }
    if(length(lmom$trim) == 1 && lmom$trim != 1) {
      warning("Attribute of TL-moments is not trim=1--can not complete parameter estimation")
      return()
    }

    LM1 <- lmom$lambdas[1]
    LM2 <- lmom$lambdas[2]

    T3 <- lmom$ratios[3]
    T4 <- lmom$ratios[4]
    T5 <- lmom$ratios[5]

    # Four parameter distributions do not normally need the
    # fifth L-moment, but for the GLD as implemented here--we do
    # This error message is to help with this fact.
    if(is.null(T5)) {
      warning("The fifth L-moment ratio TAU5 is undefined")
      return()
    }

    # THE TRIMMED L-MOMENTS OF GLD
    estTLla1 <- function(La2,La3,La4) {
      La1 <- LM1 - 6*La2*(1/((La3+3)*(La3+2)) - 1/((La4+3)*(La4+2)))
    }

    estTLla2 <- function(LM2,La3,La4) {
      D1 <- (La3+2)*(La3+3)*(La3+4)
      D2 <- (La4+2)*(La4+3)*(La4+4)
      return((LM2/6)/( La3/D1 + La4/D2 ) )
    }

    estTLtau3 <- function(K,H) {
      D1 <- (K+5)*(K+4)*(K+3)*(K+2)
      D2 <- (H+5)*(H+4)*(H+3)*(H+2)
      G  <- K*(H+4)*(H+3)*(H+2)+H*(K+4)*(K+3)*(K+2)
      t3 <- (10/9)*(K*(K-1)*D2 - H*(H-1)*D1)/((K+5)*(H+5)*G)
      return(t3)
    }
    estTLtau4 <- function(K,H) {
      D1 <- (K+6)*(K+5)*(K+4)*(K+3)*(K+2)
      D2 <- (H+6)*(H+5)*(H+4)*(H+3)*(H+2)
      G  <- K*(H+4)*(H+3)*(H+2)+H*(K+4)*(K+3)*(K+2)
      t4 <- (5/4)*(K*(K-1)*(K-2)*D2 + H*(H-1)*(H-2)*D1)/((K+6)*(H+6)*(K+5)*(H+5)*G)
      return(t4)
    }
    estTLtau5 <- function(K,H) {
      D1 <- (K+7)*(K+6)*(K+5)*(K+4)*(K+3)*(K+2)
      D2 <- (H+7)*(H+6)*(H+5)*(H+4)*(H+3)*(H+2)
      G  <- K*(H+4)*(H+3)*(H+2)+H*(K+4)*(K+3)*(K+2)
      t5 <- (7/5)*(K*(K-1)*(K-2)*(K-3)*D2 - H*(H-1)*(H-2)*(H-3)*D1)/
                                         ((K+7)*(H+7)*(K+6)*(H+6)*(K+5)*(H+5)*G)
      return(t5)
    }
    # END OF TRIMMED L-MOMENTS OF GLD

    # ORDINARY TAU3, TAU4, AND TAU5 OF GLD GIVEN PARAMETERS
    esttau3 <- function(La3,La4) {
      N1 <- La3*(La3-1)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-1)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t3 <- (N1 - N2)/(D1*D2)
      return(t3)
    }
    esttau4 <- function(La3,La4) {
      N1 <- La3*(La3-2)*(La3-1)*(La4+4)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-2)*(La4-1)*(La3+4)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+4)*(La4+4)*(La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t4 <- (N1 + N2)/(D1*D2)
      return(t4)
    }
    esttau5 <- function(La3,La4) {
      N1 <- La3*(La3-3)*(La3-2)*(La3-1)*(La4+5)*(La4+4)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-3)*(La4-2)*(La4-1)*(La3+5)*(La3+4)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+5)*(La4+5)*(La3+4)*(La4+4)*(La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t5 <- (N1 - N2)/(D1*D2)
      return(t5)
    }
    # END ORDINARY L-MOMENTS


    # Define the objective function
    fn <- function(x) {
      La3 <- x[1]
      La4 <- x[2]
      t3  <- estTLtau3(La3,La4) # TRIMMED LMOMENTS USED FOR SOLUTION, BUT CONVERSION
      t4  <- estTLtau4(La3,La4) # TO ORDINARY USED TO TEST LMOMENT VALIDITY
      ss  <- ((T3-t3)^2 + (T4-t4)^2)
      return(ss)
    }

    # This function is a scaled down version of are.pargld.valid
    # to avoid the unnecessary overhead of computing the first
    # L-moment and building the standard parameter object.
    validgld <- function(La2,La3,La4) {
      if(is.na(La2)) return(FALSE)
      if(La2 == -Inf || La2 == Inf) return(FALSE)
      #if(verbose == TRUE) cat(c("validgld-Checking Theorem 1.3.3: ",La2, La3, La4,"\n"))
      # Test that second L-moment is suitable
      for(F in seq(0,1,by=0.00001)) {
        tmp <- La2*(La3*F^(La3-1) + La4*(1-F)^(La4-1))
        #cat(c(La2,La3,La4,tmp))
        if(tmp < 0) return(FALSE)
      }

      #if(verbose == TRUE) cat(c("validgld-Checking by region: ",La2, La3, La4,"\n"))
      # ratios define the curved lines in figure1.3-1 of K&D
      ratio6 <- -1/La3
      ratio5 <-  1/La4
      # See Theorem 1.3.33 of Karian and Dudewicz and figure1.3-1
      if(La3 <= -1 && La4 >=  1) {         # REGION 1
        return(FALSE) # ordinary L-moments do not exist in REGION 1
        return(TRUE)
      }
      else if(La3 >=  1 && La4 <= -1) {    # REGION 2
        return(FALSE) # ordinary L-moments do not exist in REGION 2
        return(TRUE)
      }
      else if(La3 >= 0 && La4 >= 0) {      # REGION 3
        return(TRUE)
      }
      else if(La3 <= 0 && La4 <= 0) {      # REGION 4
        return(TRUE)
      }
      else if((La4 < ratio6 && La4 >= -1) && La3 >= 1) {  # REGION 6
        return(TRUE)
      }
      else if(La4 > ratio5 && (La3 >= -1 && La3 <= 0)) { # REGION 5
        return(TRUE)
      }
      return(FALSE)
    }

    # Are the L-moments valid?  Not fully certain that this check
    # is needed if the GLD proves to be valid through the validgld
    # function above. Early testing of pargld() before the suitabiliy
    # of the second L-moment for GLD showed that the validlmom() test
    # was needed.  Consider this for safety at any rate.
    validlmom <- function(sol,attempt) {
      LM3 <- sol$par[1]
      LM4 <- sol$par[2]
      # TL-moment (t=1) based GLD only valid for > -2 parameters--Hosking email
      if(LM3 <= -2 || LM4 <= -2) return(FALSE)
      T3 <- esttau3(LM3,LM4)
      T4 <- esttau4(LM3,LM4)
      T5 <- esttau5(LM3,LM4)
      #if(verbose == TRUE) cat(c("validlmom(region,Tau3,Tau4,Tau5)?:\n",attempt,T3,T4,T5))
      if(abs(T3) > 1)  return(FALSE)
      if(T4 < (0.25*(5*T3^2 - 1)) || T4 > 1) return(FALSE)
      if(abs(T5) > 1)  return(FALSE)
      return(TRUE)
    }

   if(is.null(initkh)) {
     g <- 13
     M <- matrix(nrow = g, ncol = 2)
     M[1,]  <- c(-2.5,2.5)   # REGION 1
     M[2,]  <- c(2.5,-2.5)   # REGION 2
     M[3,]  <- c(2.5,2.5)    # REGION 3
     M[4,]  <- c(-.5,-.5)    # REGION 4
     M[5,]  <- c(-1.5,-1.5)  # REGION 4
     M[6,]  <- c(-2.5,-2.5)  # REGION 4
     M[7,]  <- c(-3.5,-3.5)  # REGION 4
     M[8,]  <- c(-4.5,-4.5)  # REGION 4
     M[9,]  <- c(-5.5,-5.5)  # REGION 4
     M[10,] <- c(-6.5,-6.5)  # REGION 4
     M[11,] <- c(-7.5,-7.5)  # REGION 4
     M[12,] <- c(-.5,2.5)    # REGION 5
     M[13,] <- c(2.5,-.5)    # REGION 6
   }
   else {
     M <- matrix(nrow = 1, ncol = 2)
     M[1,] <- initkh
   }
   each_count            <- 0
   each_attempt          <- vector(mode = 'numeric', length = 1)
   each_initialK         <- vector(mode = 'numeric', length = 1)
   each_initialH         <- vector(mode = 'numeric', length = 1)
   each_error            <- vector(mode = 'numeric', length = 1)
   each_interpretederror <- vector(mode = 'numeric', length = 1)
   each_xi               <- vector(mode = 'numeric', length = 1)
   each_alpha            <- vector(mode = 'numeric', length = 1)
   each_kappa            <- vector(mode = 'numeric', length = 1)
   each_h                <- vector(mode = 'numeric', length = 1)
   each_t5diff           <- vector(mode = 'numeric', length = 1)

   WIDTH <- getOption("width")

   if(verbose == TRUE) {
     cat("SUMMARY OF INCREMENTAL OPTIMIZATIONS\n")
     cat("  Q(F) = X + A*(F^K + (1-F)^H)\n")
     cat(c(rep("-",WIDTH),"\n"),sep="")
     cat("Attempt      X        A        K        H            TLtau5_diff     SumSqError    Diagnostics\n")
     cat(c(rep("-",WIDTH),"\n"),sep="")
   }
   for(i in seq(1,nrow(M))) {
     # Test for NaN from the KEK guess
     if(is.nan(M[i,1]) || is.nan(M[i,2])) next

     r       <- optim(M[i,],fn) # THE MAGIC IS HERE!!!!!!!!!!
     e       <- r$value # extract error
     K       <- r$par[1] # extract the two solutions
     H       <- r$par[2]
     TL5diff <- abs(T5 - estTLtau5(K,H)) # compute difference
     A       <- estTLla2(LM2,K,H)
     X       <- estTLla1(A,K,H)

     # BEGIN SECTION FOR DETAILED OUTPUT
     if(verbose == TRUE) {
       cat(c(" ",i,sprintf("%6.6f",X),",",
                   sprintf("%6.6f",A),",",
                   sprintf("%6.6f",K),",",
                   sprintf("%6.6f",H),"   ",
                   sprintf("%6.6f",TL5diff),"  ",
                   sprintf("%8.10f",e)))
     }
     if(validlmom(r,attempt=i) == FALSE) {
       if(verbose == TRUE) cat(c("       inval:tau_r\n"))
       next
     }
     if(verbose == TRUE) cat(c("       val:tau_r--"))
     if(validgld(A,K,H) == FALSE) {
       if(verbose == TRUE) cat(c("inval:GLD\n"))
       next
     }
     if(verbose == TRUE) cat(c("val:GLD\n"))
     # END SECTION FOR DETAILED OUTPUT

     each_count                <- each_count + 1
     each_attempt[each_count]  <- i
     each_initialK[each_count] <- M[i,1]
     each_initialH[each_count] <- M[i,2]
     each_xi[each_count]       <- X
     each_alpha[each_count]    <- A
     each_kappa[each_count]    <- K
     each_h[each_count]        <- H
     each_t5diff[each_count]   <- TL5diff
     each_error[each_count]    <- e
   }
   if(verbose == TRUE) cat(c(rep("-",WIDTH),"\n"),sep="")

   EACH <- data.frame(attempt    = each_attempt,
                      x          = each_xi,
                      a          = each_alpha,
                      k          = each_kappa,
                      h          = each_h,
                      absDeltau5 = each_t5diff,
                      error      = each_error,
                      initial_k  = each_initialK,
                      initial_h  = each_initialH)
   EACH <- EACH[order(EACH$error),]

   # Preparing final best guess . . .
   para <- vector(mode="numeric", length=4)
   names(para) <- c("xi","alpha","kappa","h")
   para[1]  <- EACH$x[1]
   para[2]  <- EACH$a[1]
   para[3]  <- EACH$k[1]
   para[4]  <- EACH$h[1]
   tau5diff <- EACH$absDeltau5[1]
   error    <- EACH$error[1]

   if(result == 'best') {
     return(list(type       = 'gld',
                 para       = para,
                 error      = error,
                 absDelTau5 = tau5diff,
                 source     = "parTLgld"))
   }
   else if(result == 'dataframe') {
     if(extract > 0) {
       return(list(type       = 'gld',
                   para       = c(EACH[extract,]$x,
                                  EACH[extract,]$a,
                                  EACH[extract,]$k,
                                  EACH[extract,]$h),
                   error      = EACH[extract,]$error,
                   absDelTau5 = EACH[extract,]$absDelTau5,
                   source     = "parTLgld"))
     }
     else {
       return(EACH)
     }
   }
   else {
     warning("result argument is not 'best' or 'dataframe'")
   }
}
"pargld" <-
function(lmom, result='best', verbose=FALSE, extract=0, initkh=NULL) {

    if(length(lmom$source) == 1 && lmom$source == "TLmoms") {
      if(lmom$trim != 0) {
        warning("Attribute of TL-moments is not trim=0--can not complete parameter estimation")
        return()
      }
    }

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    LM1 <- lmom$L1
    LM2 <- lmom$L2

    T3 <- lmom$TAU3
    T4 <- lmom$TAU4
    T5 <- lmom$TAU5

    # Four parameter distributions do not normally need the
    # fifth L-moment, but for the GLD as implemented here--we do
    # This error message is to help with this fact.
    if(is.null(T5)) {
      warning("The fifth L-moment ratio TAU5 is undefined")
      return()
    }

    estla1 <- function(La2,La3,La4) {
      La1 <- LM1 - La2*(1/(La3+1) - 1/(La4+1))
    }

    estla2 <- function(LM2,La3,La4) {
      return(LM2/(La3/((La3+1)*(La3+2)) + La4/((La4+1)*(La4+2))))
    }

    esttau3 <- function(La3,La4) {
      N1 <- La3*(La3-1)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-1)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t3 <- (N1 - N2)/(D1*D2)
      return(t3)
    }
    esttau4 <- function(La3,La4) {
      N1 <- La3*(La3-2)*(La3-1)*(La4+4)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-2)*(La4-1)*(La3+4)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+4)*(La4+4)*(La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t4 <- (N1 + N2)/(D1*D2)
      return(t4)
    }
    esttau5 <- function(La3,La4) {
      N1 <- La3*(La3-3)*(La3-2)*(La3-1)*(La4+5)*(La4+4)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-3)*(La4-2)*(La4-1)*(La3+5)*(La3+4)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+5)*(La4+5)*(La3+4)*(La4+4)*(La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t5 <- (N1 - N2)/(D1*D2)
      return(t5)
    }

    # Define the objective function
    fn <- function(x) {
      La3 <- x[1]
      La4 <- x[2]
      t3  <- esttau3(La3,La4)
      t4  <- esttau4(La3,La4)
      ss  <- ((T3-t3)^2 + (T4-t4)^2)
      return(ss)
    }

    # This function is a scaled down version of are.pargld.valid
    # to avoid the unnecessary overhead of computing the first
    # L-moment and building the standard parameter object.
    validgld <- function(La2,La3,La4) {
      if(is.na(La2)) return(FALSE)
      if(La2 == -Inf || La2 == Inf) return(FALSE)
      #if(verbose == TRUE) cat(c("validgld-Checking Theorem 1.3.3: ",La2, La3, La4,"\n"))
      # Test that second L-moment is suitable
      for(F in seq(0,1,by=0.00001)) {
        tmp <- La2*(La3*F^(La3-1) + La4*(1-F)^(La4-1))
        #cat(c(La2,La3,La4,tmp))
        if(tmp < 0) return(FALSE)
      }

      #if(verbose == TRUE) cat(c("validgld-Checking by region: ",La2, La3, La4,"\n"))
      # ratios define the curved lines in figure1.3-1 of K&D
      ratio6 <- -1/La3
      ratio5 <-  1/La4
      # See Theorem 1.3.33 of Karian and Dudewicz and figure1.3-1
      if(La3 <= -1 && La4 >=  1) {         # REGION 1
        return(FALSE) # ordinary L-moments do not exist in REGION 1
        return(TRUE)
      }
      else if(La3 >=  1 && La4 <= -1) {    # REGION 2
        return(FALSE) # ordinary L-moments do not exist in REGION 2
        return(TRUE)
      }
      else if(La3 >= 0 && La4 >= 0) {      # REGION 3
        return(TRUE)
      }
      else if(La3 <= 0 && La4 <= 0) {      # REGION 4
        return(TRUE)
      }
      else if((La4 < ratio6 && La4 >= -1) && La3 >= 1) {  # REGION 6
        return(TRUE)
      }
      else if(La4 > ratio5 && (La3 >= -1 && La3 <= 0)) { # REGION 5
        return(TRUE)
      }
      return(FALSE)
    }

    # Are the L-moments valid?  Not fully certain that this check
    # is needed if the GLD proves to be valid through the validgld
    # function above. Early testing of pargld() before the suitabiliy
    # of the second L-moment for GLD showed that the validlmom() test
    # was needed.  Consider this for safety at any rate.
    validlmom <- function(sol,attempt) {
      LM3 <- sol$par[1]
      LM4 <- sol$par[2]
      # L-moment based GLD only valid for > -1 parameters--Hosking email
      if(LM3 <= -1 || LM4 <= -1) return(FALSE)
      T3 <- esttau3(LM3,LM4)
      T4 <- esttau4(LM3,LM4)
      T5 <- esttau5(LM3,LM4)
      #if(verbose == TRUE) cat(c("validlmom(region,Tau3,Tau4,Tau5)?:\n",attempt,T3,T4,T5))
      if(abs(T3) > 1)  return(FALSE)
      if(T4 < (0.25*(5*T3^2 - 1)) || T4 > 1) return(FALSE)
      if(abs(T5) > 1)  return(FALSE)
      return(TRUE)
    }

   if(is.null(initkh)) {
     g <- 14
     M <- matrix(nrow = g, ncol = 2)
     M[1,]  <- c(-2.5,2.5)   # REGION 1
     M[2,]  <- c(2.5,-2.5)   # REGION 2
     M[3,]  <- c(2.5,2.5)    # REGION 3
     M[4,]  <- c(-.5,-.5)    # REGION 4
     M[5,]  <- c(-1.5,-1.5)  # REGION 4
     M[6,]  <- c(-2.5,-2.5)  # REGION 4
     M[7,]  <- c(-3.5,-3.5)  # REGION 4
     M[8,]  <- c(-4.5,-4.5)  # REGION 4
     M[9,]  <- c(-5.5,-5.5)  # REGION 4
     M[10,] <- c(-6.5,-6.5)  # REGION 4
     M[11,] <- c(-7.5,-7.5)  # REGION 4
     M[12,] <- c(-.5,2.5)    # REGION 5
     M[13,] <- c(2.5,-.5)    # REGION 6

     # See equation 25 of Karvanen, Eriksson, and Koivunen (2002)
     # Adaptive Score Functions for Maximum Likelihood ICA
     # Journal of VLSI Signal Processing, vol. 32, pp. 83--92.
     KEKguess <- (3+7*T4+(1+98*T4+T4^2)^0.5)/(2*(1-T4))
     M[14,] <- c(KEKguess,KEKguess)
   }
   else {
     M <- matrix(nrow = 1, ncol = 2)
     M[1,] <- initkh
   }
   each_count            <- 0
   each_attempt          <- vector(mode = 'numeric', length = 1)
   each_initialK         <- vector(mode = 'numeric', length = 1)
   each_initialH         <- vector(mode = 'numeric', length = 1)
   each_error            <- vector(mode = 'numeric', length = 1)
   each_interpretederror <- vector(mode = 'numeric', length = 1)
   each_xi               <- vector(mode = 'numeric', length = 1)
   each_alpha            <- vector(mode = 'numeric', length = 1)
   each_kappa            <- vector(mode = 'numeric', length = 1)
   each_h                <- vector(mode = 'numeric', length = 1)
   each_t5diff           <- vector(mode = 'numeric', length = 1)

   WIDTH <- getOption("width")

   if(verbose == TRUE) {
     cat("SUMMARY OF INCREMENTAL OPTIMIZATIONS\n")
     cat("  Q(F) = X + A*(F^K + (1-F)^H)\n")
     cat(c(rep("-",WIDTH),"\n"),sep="")
     cat("Attempt      X        A        K        H            tau5_diff     SumSqError    Diagnostics\n")
     cat(c(rep("-",WIDTH),"\n"),sep="")
   }
   for(i in seq(1,nrow(M))) {
     # Test for NaN from the KEK guess
     if(is.nan(M[i,1]) || is.nan(M[i,2])) next

     r      <- optim(M[i,],fn) # THE MAGIC IS HERE!!!!!!!!!!
     e      <- r$value # extract error
     K      <- r$par[1] # extract the two solutions
     H      <- r$par[2]
     T5diff <- abs(T5 - esttau5(K,H)) # compute difference
     A      <- estla2(LM2,K,H)
     X      <- estla1(A,K,H)

     # BEGIN SECTION FOR DETAILED OUTPUT
     if(verbose == TRUE) {
        cat(c(" ",i,"  ",
           sprintf("%6.6f",X),", ",
           sprintf("%6.6f",A),", ",
           sprintf("%6.6f",K),", ",
           sprintf("%6.6f",H),"   ",
           sprintf("%6.6f",T5diff),"  ",
           sprintf("%8.10f",e)),sep="")
     if(validlmom(r,attempt=i) == FALSE) {
       if(verbose == TRUE) cat(c("       inval:tau_r\n"))
       next
     }
     if(verbose == TRUE) cat(c("       val:tau_r--"))
     if(validgld(A,K,H) == FALSE) {
       if(verbose == TRUE) cat(c("inval:GLD\n"))
       next
     }
     if(verbose == TRUE) cat(c("val:GLD\n"))
     # END SECTION FOR DETAILED OUTPUT

     each_count                <- each_count + 1
     each_attempt[each_count]  <- i
     each_initialK[each_count] <- M[i,1]
     each_initialH[each_count] <- M[i,2]
     each_xi[each_count]       <- X
     each_alpha[each_count]    <- A
     each_kappa[each_count]    <- K
     each_h[each_count]        <- H
     each_t5diff[each_count]   <- T5diff
     each_error[each_count]    <- e
   }
   if(verbose == TRUE) cat(c(rep("-",WIDTH),"\n"),sep="")

   EACH <- data.frame(attempt    = each_attempt,
                      x          = each_xi,
                      a          = each_alpha,
                      k          = each_kappa,
                      h          = each_h,
                      absDelTau5 = each_t5diff,
                      error      = each_error,
                      initial_k  = each_initialK,
                      initial_h  = each_initialH)
   EACH <- EACH[order(EACH$error),]

   # Preparing final best guess . . .
   para <- vector(mode="numeric", length=4)
   names(para) <- c("xi","alpha","kappa","h")
   para[1]  <- EACH$x[1]
   para[2]  <- EACH$a[1]
   para[3]  <- EACH$k[1]
   para[4]  <- EACH$h[1]
   tau5diff <- EACH$absDelTau5[1]
   error    <- EACH$error[1]

   if(result == 'best') {
     return(list(type       = 'gld',
                 para       = para,
                 error      = error,
                 absDelTau5 = tau5diff,
                 source     = "pargld"))
   }
   else if(result == 'dataframe') {
     if(extract > 0) {
       return(list(type       = 'gld',
                   para       = c(EACH[extract,]$x,
                                  EACH[extract,]$a,
                                  EACH[extract,]$k,
                                  EACH[extract,]$h),
                   error      = EACH[extract,]$error,
                   absDelTau5 = EACH[extract,]$absDelTau5,
                   source     = "pargld"))
     }
     else {
       return(EACH)
     }
   }
   else {
     warning("result argument is not 'best' or 'dataframe'")
   }
}
