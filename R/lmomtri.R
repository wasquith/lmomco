"lmomtri" <-
function(para, paracheck=TRUE) {
    # DEPRECATED MAY 2026 nmom <- match.arg(nmom)
    # DEPRECATED MAY 2026 nmom <- as.numeric(nmom)

    if(paracheck) {
       if(! are.partri.valid(para)) return()
    }
    attributes(para$para) <- NULL

    MIN  <- para$para[1]
    MODE <- para$para[2]
    MAX  <- para$para[3]
    A <- (MAX-MIN); B <- (MODE-MIN); C <- (MAX-MODE)
    ABS_LEFT_OVER_RIGHT <- abs(B/C) # as part of May 2026 changes

    L1 <- (MIN + MODE + MAX)/3
    L2 <- ((MIN-MODE)^2/(MAX-MIN) - (MIN+MODE) + 2*MAX)/15

    tmp1 <- 2*MIN*(B/A)^3 + (12/7)*A^(-3)*B^4 + 2*MAX-2*MAX*(B/A)^3 - (12/7)*A^(-3)*C^4
    tmp2 <- (12/5)*2*A^(-2)*C^3 - (12/3)*A^(-1)*C^2
    tmp3 <- -3*((1/A^2) * (MIN*B^2 + (4/5)*B^3 ) + MAX - (1/A^2) * (MAX*B^2 - (4/5)*C^3 + (4/3)*A*C^2))
    L3   <- tmp1 + tmp2 + tmp3 + L1
    LCV  <- L2/L1
    TAU3 <- L3/L2

    "trirats" <- function(x) {
      myPlst <- list(t3 = c( 2.09473699830127e-07, -0.258427192321678, 0.000813523008402076,
                             0.267735800912562, 0.083539386009353, 0.103792393415911,
                             1.1235299039319, 2.50893936592746, 2.95985604098822,
                             2.23014491179716, 1.13782291985567, 0.393669866501206,
                             0.0874340488117483, 0.0102291668349848, 0, -0.000126355387612369,
                             0, 2.14008028158768e-06, 0, -2.18232884781946e-08),
                     t4 = c( 0.0901349547176836, -9.21544372657868e-05, -0.123029444766339,
                            -0.0267561923267532, 0.0327240614465926, -0.656790348277229,
                            -1.85192915044293, -2.41342071512778, -1.90097991627064,
                            -0.969596294880085, -0.317484445715324, -0.0603555992143965,
                            -0.00416931903790304, 0.000442849305814251, 0, -2.44917241833515e-05,
                             0, 6.74776329981755e-07, 0, -9.11225652915777e-09),
                     t5 = c( 2.55962420237385e-06, 0.0047273193737089, 0.00673917574549263,
                            -0.0375002891059526, 0.339063191566633, 1.21158859349206,
                             1.62612385165448, 0.970492776955565, -0.0568087391623817,
                            -0.547508600661963, -0.448099986683432, -0.196524432565279,
                            -0.0505100126115241, -0.00652735793993698, 0, 9.16188021590211e-05,
                             0, -1.68350058057341e-06, 0, 1.82509292864212e-08),
                     t6 = c( 0.0164331340923569, -0.000336312205790996, 0.0355477156765414,
                             0.00262341574829825, -0.00535564820450852, 1.03704335704812,
                             3.68197811052433, 6.24245398891586, 6.51916134071964,
                             4.56843605968637, 2.20913127358098, 0.730060999480297,
                             0.155406018803995, 0.0174567726997667, 0, -0.000199592834921634,
                             0, 3.14937242578995e-06, 0, -3.01770503062676e-08))
      xp <- x; xp[x > 1] <- 1/x[x > 1]; xp[xp < -3] <- -3; xp <- log10(xp)
      strs <- c("t3", "t4", "t5", "t6"); zzz <- list(left_over_right=x)
      for(str in strs) {
        a  <- get(str, myPlst)
        mx <- stats::poly(xp, degree=length(a)-1, raw=TRUE)
        b  <- a[ 1]; m <- a[-1]
        ox <- matrix(nrow=nrow(mx), ncol=ncol(mx))
        for(i in seq_len(ncol(mx))) ox[,i] <- m[i]*mx[,i]
        yp <- rowSums(ox) + b
        if(str %in% c("t3", "t5")) yp[x > 1] <- -yp[x > 1]
        zzz[[str]] <- yp
      }
      return(zzz)
    }
    zzz <- trirats(ABS_LEFT_OVER_RIGHT)

    z <- list(lambdas=c(L1, L2, L3,   zzz[["t4"]]*L2, zzz[["t5"]]*L2, zzz[["t6"]]*L2),
              ratios=c(NA, LCV, TAU3, zzz[["t4"]],    zzz[["t5"]],    zzz[["t6"]]),
              trim=0, rightrim=0, leftrim=0, tau3poly=zzz[["t3"]],
              source="lmomtri")
    return(z)
    # ------------------ abandoned May 25, 2026, in favor of the polynomial approximations above
    E11  <- L1
    E22  <- L2 + E11
    E33  <- (L3 + 3*E22 - E11)/2
    nmom <- NULL # part of deprecation but to keep R CRAN checks passing
    if(nmom == 3) {
       z <- list(lambdas=c(L1, L2, L3, NA, NA),
                 ratios=c(NA, LCV, TAU3, NA, NA),
                 trim=0, rightrim=0, leftrim=0, E33err=NA,
                 source="lmomtri")
       return(z)
    } else {
       E33p <- E44 <- E55 <- NULL
       try(E33p <- expect.max.ostat(3, para=para, qua=quatri), silent=TRUE)
       try(E44  <- expect.max.ostat(4, para=para, qua=quatri), silent=TRUE)
       try(E55  <- expect.max.ostat(5, para=para, qua=quatri), silent=TRUE)
       if(is.null(E33p)) {
          warning("E33p could not be computed by expect.max.ostat() numerical integration ",
                  "so E33err will be NA")
          E33p <- NA
       }
       if(is.null(E44)) {
          warning("E44 could not be computed by expect.max.ostat() numerical integration ",
                  "so Tau4 will be NA")
          E44 <- NA
       }
       if(is.null(E55)) {
          warning("E55 could not be computed by expect.max.ostat() numerical integration ",
                  "so Tau5 will be NA")
          E55 <- NA
       }
       L4 <-           5*E44 - 10*E33 +  6*E22 - 1*E11
       L5 <- 14*E55 - 35*E44 + 30*E33 - 10*E22 + 1*E11
       TAU4 <- L4/L2; TAU5 <- L5/L2

       E33percent.error <- 100*((E33 - E33p)/E33)

       z <- list(lambdas=c(L1, L2, L3, L4, L5),
                 ratios=c(NA, LCV, TAU3, TAU4, TAU5),
                 trim=0, rightrim=0, leftrim=0, E33err=E33percent.error,
                 source="lmomtri")
       # DISABLED return(z)
       return("disabled")
    }
}

