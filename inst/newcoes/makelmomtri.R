df <- NULL
for(i in seq_len(4000)) {
  l <- -runif(1, min=0, max=1000); r <- runif(1, min=0, max=1000)
  lmr <- theoLmoms(vec2par(c(l,0, r), type="tri"), nmom=6, subdivisions=200)
  t3 <- lmr$ratios[3]; t4 <- lmr$ratios[4]
  t5 <- lmr$ratios[5]; t6 <- lmr$ratios[6]
  df <- rbind(df, data.frame(l=l, r=r, left_over_right=abs(l/r), t3=t3, t4=t4, t5=t5, t6=t6))
}
df <- df[order(df$t3),]

left_over_right <- c(10^seq(-3, -1, by=0.01), seq(0.1, 10, by=0.005), 10^seq(1, 3, by=0.01))
left_over_right <- sort( unique(left_over_right) )
t3x <- approx(log10(df$left_over_right), df$t3, xout=log10(left_over_right), rule=2)$y
plot(left_over_right, t3x, log="x", type="l")

t4x <- approx(log10(df$left_over_right), df$t4, xout=log10(left_over_right), rule=2)$y
plot(left_over_right, t4x, log="x", type="l")

t5x <- approx(log10(df$left_over_right), df$t5, xout=log10(left_over_right), rule=2)$y
plot(left_over_right, t5x, log="x", type="l")

t6x <- approx(log10(df$left_over_right), df$t6, xout=log10(left_over_right), rule=2)$y
plot(left_over_right, t6x, log="x", type="l")

dfo <- data.frame(left_over_right, t3=t3x, t4=t4x, t5=t5x, t6=t6x)

t3s <- unique(c(seq(min(dfo$t3), max(dfo$t3), by=0.001), max(dfo$t3)))
t4s <- approx(dfo$t3, dfo$t4, xout=t3s, rule=2)$y

plotlmrdia(lmrdia(), xlim=c(-0.2, .2), ylim=c(0,0.2))
lines(dfo$t3, dfo$t4, lwd=8, col="grey80")
lines(t3s, t4s, col="red1")
t4fromt3 <- lm(t4s~I(t3s^2)+I(t3s^4)+I(t3s^6)+I(t3s^8))
lines(t3s, predict(t4fromt3))
print(t4fromt3$coefficients, 16)
#print(t4fromt3$coefficients, 8)
#           (Intercept)               I(t3s^2)               I(t3s^4)               I(t3s^6)
# 9.012821605998006e-02 -1.789333565397665e+00 -1.448666575025029e+01  6.190828635274557e+02
#              I(t3s^8)
#-2.965804718201972e+04
print((abs(min(dfo$t3)) + max(dfo$t3))/2, 16) # 0.1428571650715543
numbers::contfrac(0.1428571650715543,  tol=1E-7) # 1 / 7
print(max(dfo$t4), 16) # 0.09013546774995014
numbers::contfrac(0.09013546774995014, tol=1E-7) # 233 / 2585
print(min(dfo$t4), 16) # 0.04761902517520374
numbers::contfrac(0.04761902517520374, tol=1E-7) # 1 / 21






wnt <- dfo$left_over_right <= 1
x <- log10(dfo$left_over_right[wnt])
strs <- c("t3", "t4", "t5", "t6")
strs <- strs[1:4]

Penv <- new.env()
for(str in strs) {
  message("---- ", toupper(str), " ------------ Polynomial notation [m]f(x) ------------ ")
  y <- dfo[,str] # response variable
  y <- y[wnt]
  ERR <- +Inf; KEY <- ""; A <- NA # characteristics of the optimal solution for str
  for(m in 1:20) {
      key <- paste0("[", m, "]"); message(key, appendLF=FALSE)
      lmmod <- lm(y~poly(x, degree=m, raw=TRUE))
      a <- coefficients(lmmod)
      yp <- predict(lmmod)
      err <- sum(1*(y - yp)^2, na.rm=TRUE) # accumulate a square error
      if(err < ERR) { # if we have found a new minimum error, let us preserve those results
          ERR <- err; KEY <- key
          a[is.na(a)] <- 0
          assign(str, a, envir=Penv)
          message("*", appendLF=FALSE)
          txt <- paste0(key, " err=", round(err, digits=12))
          plot( 10^x, y, log="x", pch=21, col=grey(0.7), bg="white", xlab="|l/r|", ylab=str)
          lines(10^x, yp, lwd=3, col="forestgreen"); mtext(txt, line=1); #stop()
      }
    message("")
  }
}


Plst <- as.list(Penv)
txt <- paste0("myPlst <- list(\n",  paste(
 sapply(seq_len(length(Plst)), function(k) paste0("         ", names(Plst)[k],
                                     paste0(" = c(", paste(get(names(Plst)[k], Plst),
                collapse=", ")), ")")),
                collapse=",\n"), ")")
message(txt)
eval(str2expression(txt))

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
    xp <- x; xp[x > 1] <- 1/x[x > 1]
    xp[xp < -3] <- -3; xp <- log10(xp)
    strs <- c("t3", "t4", "t5", "t6")
    zzz <- list(left_over_right=x)
    for(str in strs) {
      a  <- get(str, myPlst)
      mx <- poly(xp, degree=length(a)-1, raw=TRUE)
      b  <- a[ 1]; m <- a[-1]
      ox <- matrix(nrow=nrow(mx), ncol=ncol(mx))
      for(i in seq_len(ncol(mx))) ox[,i] <- m[i]*mx[,i]
      yp <- rowSums(ox) + b
      if(str %in% c("t3", "t5")) yp[x > 1] <- -yp[x > 1]
      zzz[[str]] <- yp
    }
    return(zzz)
  }


zzz <- trirats(dfo$left_over_right)

plot(dfo$t3, dfo$t4, log="", col=grey(0.8))
lines(zzz[["t3"]], zzz[["t4"]])
