setwd(dirname(rstudioapi::getSourceEditorContext()$path))

####################################################################################################
#
#  ad88888ba  8b        d8  88b           d88     8888888ba,    88   ad8888ba  8888888888  ad8888ba
# d8"     "8b  Y8,    ,8P   888b         d888     88     `"8b   88  d8"    "8b     88     d8"    "8b
# Y8,           Y8,  ,8P    88`8b       d8'88     88       `8b  88  Y8,            88     Y8,
# `Y8aaaaa,      "8aa8"     88 `8b     d8' 88     88        88  88  `Y8aaaa,       88     `Y8aaaa,
#   `"""""8b,     `88'      88  `8b   d8'  88     88        88  88     `"""8b,     88       `"""8b,
#         `8b      88       88   `8b d8'   88     88        8P  88         `8b     88            `8b
# Y8a     a8P      88       88    `888'    88     88     .a8P   88  Y8a    a8P     88     Y8a    a8P
#  "Y88888P"       88       88     `8'     88     8888888Y"'    88   "Y8888P"      88      "Y8888P"
#
####################################################################################################

library(lmomco)

# L1, first L-moment, arithmetic mean
# L2, second L-moment, L-variation
# T3, third L-moment ratio (L3/L2), L-skew
# T4, third L-moment ratio (L4/L2), L-kurtosis
# T5, third L-moment ratio (L5/L2), L-cinco (yes, some in literature, Mel Schaefer)
# T6, third L-moment ratio (L6/L2), named (L-seis ?)

linety <- c(4, 1, 1, 1)
colors <- c("red", "forestgreen", "black",  "blue")
dtypes <- c("aep4",  "gld",       "st3",    "pdq4")

T4s <- sort( c(seq(-0.24, 0.99, by=0.01), 0.999, 0, -0.005, +0.005, -0.249) ) # strategic Tau4 values
T4s <- round( T4s, digits=4) # trim decimals beyond those we likely will be caring about to tidy tables
T4s <- unique(T4s) # insurance against any double counting
L1 <- 0      # arithmetic mean of zero w/o loss of generality
L2 <- 1      # a unit of L-variation w/o loss of generality for location-scale families
T6df <- NULL # a data frame to construct for each Tau4 and each distribution
for(T4 in T4s) {
  message(T4, appendLF=FALSE)
  ggg <- c(L1, L2, 0, T4, 0)
  lmr <- lmomco::vec2lmom(ggg, checklmom=FALSE) # structure needed by lmomco package
  for(dtype in dtypes) { # for each of the distribution types
    if(dtype == "pdq4" & T4 > 0.845) next # beyond apparent computational limits in R for the math of PDQ4
    para <- NULL # initialization for the try() logic
    try(para <- lmomco::lmom2par(lmr, type=dtype, checklmom=FALSE))
    if(  is.null(para))       next # some type of "failure" or out of bounds for given distribution
    if(any(is.na(para$para))) next # some type of "failure" or out of boudns for given distribution
    if(dtype == "aep4" & para$type == "kap") next # default has AEP4 fitting Kappa if below the line.
    message(paste0("-", dtype), appendLF=FALSE)
    if(dtype == "gld") message(paste0("(", round(c(para$para[3], para$para[4]), digits=5), ")"),
                                appendLF=FALSE)
    blmr <- lmomco::par2lmom(para) # lmomco conversion of parameters back to L-moments
    if(length(blmr$ratios) == 4) lmr$ratios <- c(lmr$ratios, NA, NA) # not all returns have vector to
    if(length(blmr$ratios) == 5) lmr$ratios <- c(lmr$ratios, NA) # 5 or 6 L-moments, so pad here
    tlmr <- lmomco::theoTLmoms(para, nmom=6) # theoretical, num. integration of quantile function
    #print(c(dtype, blmr$ratios[c(4,6)], tlmr$ratios[c(4,6)]))
    df <- data.frame(tau4=T4, type=dtype, special_parameter=NA,
                     tau4bak=blmr$ratios[4], tau6bak=blmr$ratios[6],
                     tau4int=tlmr$ratios[4], tau6int=tlmr$ratios[6])
    if(dtype == "gld") {
      df$special_parameter <- mean(c(para$para[3], para$para[4]))
    } else if(dtype == "pdq4" | dtype == "st3") {
      df$special_parameter <- mean(  para$para[3]  )
    }
    T6df <- rbind(T6df, df) # data assembly
  }
  message("")
}

# Now for some minor post-assembly manipulation of the data frame of the Tau6 values
T6df$tau6 <- T6df$tau6bak # store into a canonical column
T6df$tau6[is.na(T6df$tau6)] <- T6df$tau6int[is.na(T6df$tau6)] # back file if lmomco does not have
# native parameter ---> to 6 L-moment return, so must use the numerical integration versions
T6df <- T6df[ , c("type", "special_parameter", "tau4", "tau6", "tau4bak", "tau6bak", "tau4int", "tau6int")] # re-order
T6df$tau6err <- abs(T6df$tau6bak - T6df$tau6int) # compute an error for later study, if needed
T6df <- T6df[order(T6df$tau4), ]  # sort the table by L-kurtosis
T6df$type <- as.factor(T6df$type) # perhaps for some future operations, convert to factor
T6df$col  <- NA # initialize a color and then set them
T6df$lty  <- NA # initialize a color and then set them
for(dtype in dtypes) T6df$col[T6df$type == dtype] <- colors[dtypes == dtype]
for(dtype in dtypes) T6df$lty[T6df$type == dtype] <- linety[dtypes == dtype]



####################################################################################################
#
#   888888888888     ,d8         888888888888  ad8888ba,     88888888ba,         db    888888888888
#        88        ,d888              88      8P'    "Y8     88      `"8b       d88b        88
#        88      ,d8" 88              88     d8              88        `8b     d8'`8b       88
#        88    ,d8"   88              88     88,dd888bb,     88         88    d8'  `8b      88
#        88  ,d8"     88    aaaaaaaa  88     88P'    `8b     88         88   d8YaaaaY8b     88
#        88  8888888888888  """"""""  88     88       d8     88         8P  d8""""""""8b    88
#        88           88              88     88a     a8P     88      .a8P  d8'        `8b   88
#        88           88              88      "Y88888P"      88888888Y"'  d8'          `8b  88
#
####################################################################################################
# There might exist a need to have one canonical parent-like data object containing the relations
# between Tau4 and Tau6. So, let us begin assembling one with important distributions, cross-check
# against other sources when possible and otherwise tinker such that proper behavior is seen near
# Tau4 = 0 (uniform distribution) and (or) for cases of exactly L2 = 1.
# ***** NORMAL DISTRIBUTION *****
NOR <- lmomco::theoTLmoms(vec2par(c(0, 1/sqrt(pi)), type="nor"), nmom=6, verbose=TRUE)
NOR$ratios[4] <- lmomco::par2lmom(vec2par(c(0, 1/sqrt(pi)), type="nor"))$ratios[4]
TAU46dists <- list(nor=c(NOR$ratios[4], NOR$ratios[6]))

# ***** TUKEY LAMBDA DISTRIBUTION ***** (This is not a constant L2 distribution, see lmomco::GLD)
# https://en.wikipedia.org/wiki/Tukey_lambda_distribution
Lambdas <- sort(c(0.001, -0.001, seq(-1, 1, by=0.005))) # see getting "close to zero" in the vector
TukeyL2 <- (2/Lambdas) * (- 1/(1+Lambdas) +  2/(2+Lambdas))
TukeyL4 <- (2/Lambdas) * (- 1/(1+Lambdas) + 12/(2+Lambdas) -  30/(3+Lambdas) +  20/(4+Lambdas))
TukeyL6 <- (2/Lambdas) * (- 1/(1+Lambdas) + 30/(2+Lambdas) - 210/(3+Lambdas) + 560/(4+Lambdas) -
                                                             630/(5+Lambdas) + 252/(6+Lambdas))
TukeyT4 <- TukeyL4 / TukeyL2
TukeyT6 <- TukeyL6 / TukeyL2


# Symmetrical L-moment Ratio Tables from Richard Vogel (March 25, 2025)
# We will see that the Stable is "new" per se to the lmomco infrastructure, the Student t with
# df matches the Student 3t from lmomco, the Tukey matches Wikipedia entry, and the Power
# Exponential is the AEP4 from lmomco with Tau3 = 0.
Stable  <- read.table("StableDistribution.txt", sep="|", header=TRUE)
Student <- read.table("StudentT.txt",           sep="|", header=TRUE)
Tukey   <- read.table("SymTukeyLambda.txt",     sep="|", header=TRUE)
PowExp  <- read.table("PowerExponential.txt",   sep="|", header=TRUE)

plot(                                      T6df$tau6[T6df$type == "aep4"],
     approx(PowExp$tau4, PowExp$tau6, xout=T6df$tau4[T6df$type == "aep4"])$y)
abline(0, 1) # equal value line
TAU46dists$pwrexp <- data.frame(alpha=PowExp$alpha, tau4=PowExp$tau4, tau6=PowExp$tau6)
TAU46dists$aep4   <- data.frame(tau4=T6df$tau4[T6df$type == "aep4"], tau6=T6df$tau6[T6df$type == "aep4"])


plot(                                        T6df$tau6[T6df$type == "st3"],
     approx(Student$tau4, Student$tau6, xout=T6df$tau4[T6df$type == "st3"])$y)
abline(0, 1) # equal value line
TAU46dists$st2 <- data.frame(df=Student$df, tau4=Student$tau4, tau6=Student$tau6)
TAU46dists$st3 <- data.frame(tau4=Student$tau4, tau6=Student$tau6)


plot(                                    TukeyT6,
     approx(Tukey$tau4, Tukey$tau6, xout=TukeyT4)$y)
abline(0, 1) # equal value line
TAU46dists$tukeylam <- data.frame(lambda=Lambdas, lambda2=TukeyL2, tau4=TukeyT4, tau6=TukeyT6)
TAU46dists$tukeylam$tau4[1] <- TAU46dists$tukeylam$tau6[1] <- 1
TAU46dists$tukeylam$lambda2[TAU46dists$tukeylam$lambda == 0] <-
                  approx(   TAU46dists$tukeylam$lambda, y=TAU46dists$tukeylam$L2,   xout=0)$y
TAU46dists$tukeylam$lambda2[TAU46dists$tukeylam$lambda == 0] <- 1
TAU46dists$tukeylam$tau4[   TAU46dists$tukeylam$lambda == 0] <-
                  approx(   TAU46dists$tukeylam$lambda, y=TAU46dists$tukeylam$tau4, xout=0)$y
TAU46dists$tukeylam$tau6[   TAU46dists$tukeylam$lambda == 0] <-
                  approx(   TAU46dists$tukeylam$lambda, y=TAU46dists$tukeylam$tau6, xout=0)$y


TAU46dists$symstable <- Stable
j <- TAU46dists$symstable[0,]
j[1,] <- c(NA, 1, NA, NA)
TAU46dists$symstable <- rbind(TAU46dists$symstable, j)
TAU46dists$symstable <-   TAU46dists$symstable[order(TAU46dists$symstable$lambda2), ]
TAU46dists$symstable$tau4[TAU46dists$symstable$lambda2 == 1] <- # this is precise unity
                  approx( TAU46dists$symstable$lambda2, y=TAU46dists$symstable$tau4, xout=1)$y
TAU46dists$symstable$tau6[TAU46dists$symstable$lambda2 == 1] <- # this is precise unity
                  approx( TAU46dists$symstable$lambda2, y=TAU46dists$symstable$tau6, xout=1)$y


TAU46dists$gld  <- data.frame(khparas=T6df$special_parameter[T6df$type == "gld" ],
                              tau4=T6df$tau4[T6df$type == "gld" ],
                              tau6=T6df$tau6[T6df$type == "gld" ])
TAU46dists$pdq4 <- data.frame(kappa=T6df$special_parameter[T6df$type == "pdq4" ],
                              tau4=T6df$tau4[T6df$type == "pdq4"],
                              tau6=T6df$tau6[T6df$type == "pdq4"])


# We have now completed a list() object containing lookup tables of the Tau4 and Tau6 of
# symmetrical distributions for reference and creation (if needed externally) of T4-T6 diagrams.
# Note, the ST3 in lmomco has internally polynomial-based versions of the T4 and T6, see lmomst3().
####################################################################################################



####################################################################################################
#
#   888888888888     ,d8         888888888888  ad8888ba,        88888888ba,    88         db
#        88        ,d888              88      8P'    "Y8        88      `"8b   88        d88b
#        88      ,d8" 88              88     d8                 88        `8b  88       d8'`8b
#        88    ,d8"   88              88     88,dd888bb,        88         88  88      d8'  `8b
#        88  ,d8"     88    aaaaaaaa  88     88P'    `8b        88         88  88     d8YaaaaY8b
#        88  8888888888888  """"""""  88     88       d8        88         8P  88    d8""""""""8b
#        88           88              88     88a     a8P        88      .a8P   88   d8'        `8b
#        88           88              88      "Y88888P"         88888888Y"'    88  d8'          `8b
#
####################################################################################################
apptxt <- c(" L2 = 1 and Tau3 set to zero (PowerExponential)",
            " L2 = 1 and Tau3 set to zero",
            " L2 = 1 and Tau3 defined as zero",
            " L2 = 1 and Tau3 defined as zero")
par(lend=2, ljoin=1, mgp=c(2.3, 0.8, 0))
plot(c(-0.25, 1), c(0, 1), type="n", las=1, xlim=c(0,0.4), ylim=c(0, 0.4),
     xlab="L-kurtosis (Tau4), dimensionless", ylab="Tau6 (6th L-moment ratio), dimensionless")
#points(2.502947e-01, 1.128522e-01, pch=16, cex=3, col=grey(0.5)) # MRVA residuals, 155960 monthlies
# MRVA skew is -3.605823e-02
axis(3, axTicks(1), labels=FALSE, mgp=c(2, 0.7, 0), lwd=0, lwd.ticks=1)
axis(4, axTicks(2), labels=FALSE, mgp=c(2, 0.7, 0), lwd=0, lwd.ticks=1)


for(dtype in dtypes) {
  tmp <- T6df[T6df$type == dtype, ]
  lines(tmp$tau4, tmp$tau6, col=tmp$col, lwd=1.5, lty=tmp$lty)
}


lines(TAU46dists$symstable$tau4, TAU46dists$symstable$tau6, lwd=1.5, lty=5, col="tan2"      )
lines(TAU46dists$tukeylam$tau4,  TAU46dists$tukeylam$tau6,  lwd=1.5, lty=2, col="indianred4")

wnt <- TAU46dists$symstable$lambda2 == 1 # this precise matching is what we wanted perfect 1 on the table
points(TAU46dists$symstable$tau4[wnt], TAU46dists$symstable$tau6[wnt], cex=1.4, pch=15, col="tan2"      )
wnt <- TAU46dists$tukeylam$lambda2  == 1 # this precise matching is what we wanted perfect 1 on the table
points(TAU46dists$tukeylam$tau4[ wnt], TAU46dists$tukeylam$tau6[ wnt], cex=1.4, pch=16, col="indianred4")

points(TAU46dists$nor[1], TAU46dists$nor[2], pch=10, lwd=1.3, cex=1.8, col=grey(0.3)) # NORMAL DIST


txt <- paste0(toupper(dtypes), apptxt, " (lmomco package)")
txt <- c(txt, "Symmetric stable distribution (L2 != 1), symbol at L2 = 1")
txt <- c(txt, "Tukey lambda distribution (L2 != 1), symbol at L2 = 1")
txt <- c(txt, "Normal distribution")
legend("topleft", txt, bty="o", cex=0.8, lwd=1.5, seg.len=3.5, inset=0.01, box.lty=0,
       lty=   c(rep( 1,  length(colors)),  5,   2, NA), y.intersp = 1.1,
       pt.cex=c(rep( NA, length(colors)),  1.4, 1.4, 1.8),
       pt.lwd=c(rep( NA, length(colors)),  1, 1, 1.3),
       col=   c(               colors, "tan2", "indianred4", grey(0.3)),
       pch=   c(rep(NA, length(colors)), 15,  16,   10))
mtext("L-moment Ratio Diagram of Symmetrical Distributions (L-skew = 0)", line=1, font=2)


####################################################################################################
#
#    88888888ba,    88888888888  888b      88   ad88888ba   88  888888888888  8b        d8
#    88      `"8b   88           8888b     88  d8"     "8b  88       88        Y8,    ,8P
#    88        `8b  88           88 `8b    88  Y8,          88       88         Y8,  ,8P
#    88         88  88aaaaa      88  `8b   88  `Y8aaaaa,    88       88          "8aa8"
#    88         88  88"""""      88   `8b  88    `"""""8b,  88       88           `88'
#    88         8P  88           88    `8b 88          `8b  88       88            88
#    88      .a8P   88           88     `8888  Y8a     a8P  88       88            88
#    88888888Y"'    88888888888  88      `888   "Y88888P"   88       88            88
#
####################################################################################################
FF <- seq(0.01, 0.99, by=0.01) # nonexceedance probabilities
T4 <- 0.3 # NOR, GLD, AEP4, PDQ4, ST3
ggg <- c(L1, L2, 0, T4, 0) # vector of L-moments, Lmoments package needs this vector with L4 but
                           # because L2 = 1, we can use T4, but be careful is L2 != 1
lmr <- lmomco::vec2lmom(ggg, checklmom=FALSE) # structure needed by lmomco package

para <- lmom2par(lmr, type="gld"); gldt6 <- lmomco::theoTLmoms(para, nmom=6)$ratios[6]
plot(qlmomco(FF, para), dlmomco(qlmomco(FF, para), para), type="l", las=1,
     ylim=c(0, 1), lwd=3, col="forestgreen", ylab="Probability density",
     xlab="Symmetrical quantile (mean zero and L2 = 1) from 0.01 to 0.99 nonexceedance probability",)
axis(3, axTicks(1), labels=FALSE, mgp=c(2, 0.7, 0), lwd=0, lwd.ticks=1)
axis(4, axTicks(2), labels=FALSE, mgp=c(2, 0.7, 0), lwd=0, lwd.ticks=1)

para   <- lmomco::lmom2par(lmr, type="aep4")
aep4t6 <- lmomco::theoTLmoms(para, nmom=6)$ratios[6]
lines(        qlmomco(FF, para),
      dlmomco(qlmomco(FF, para), para), lty=4, col="red")

para   <- lmomco::lmom2par(lmr, type="pdq4")
pdq4t6 <- lmomco::theoTLmoms(para, nmom=6)$ratios[6]
lines(        qlmomco(FF, para),
      dlmomco(qlmomco(FF, para), para), col="blue")

para   <- lmomco::lmom2par(lmr, type="st3" )
st3t6  <- lmomco::theoTLmoms(para, nmom=6)$ratios[6]
lines(        qlmomco(FF, para),
      dlmomco(qlmomco(FF, para), para), col="black")

myt4s <- sprintf("%0.3f", round(c(aep4t6, gldt6, pdq4t6, st3t6), digits=3))

apptxt <- c(" L2 = 1 and Tau3 set to zero (PowerExponential)",
            " L2 = 1 and Tau3 set to zero",
            " L2 = 1 and Tau3 defined as zero",
            " L2 = 1 and Tau3 defined as zero")
txt <- paste0(toupper(dtypes), apptxt, " (lmomco package) (Tau6 = ", myt4s[1:4], ")")
legend("topleft", txt, bty="o", cex=0.8, lwd=c(1, 3, 1, 1), seg.len=3.5, inset=0.01, box.lty=0,
       lty=   c(4, rep( 1, length(colors)-1)), y.intersp = 1.1,
       pt.cex=c(rep( 1, length(colors))),
       col=   c(               colors),
       pch=   NA)
mtext(paste0("Symmetrical (Tau3 = 0) Probability Densities for Tau4 = ", T4,
             " for L-scale (L2) = 1"), line=1, font=2)



stop()





lmr <- vec2lmom(c(L1, 10*L2, 0, 0.6, 0), checklmom=FALSE)
slmr <- lmoms(rlmomco(1E5, lmom2par(lmr, type="gld")), nmom=8)
points(slmr$ratios[4], slmr$ratios[6])

lmr <- vec2lmom(c(L1, 10*L2, 0, 0, 0), checklmom=FALSE)
slmr <- lmoms(rlmomco(1E5, lmom2par(lmr, type="pdq4")), nmom=6)
points(slmr$ratios[4], slmr$ratios[6])


lmr <- vec2lmom(c(L1, 10*L2, 0, 0.3, 0), checklmom=FALSE)
slmr <- lmoms(rlmomco(1E5, lmom2par(lmr, type="st3")), nmom=6)
points(slmr$ratios[4], slmr$ratios[6])



lmr <- vec2lmom(c(L1, 10*L2, 0, 0.3, 0), checklmom=FALSE)
slmr <- lmoms(rlmomco(1E5, lmom2par(lmr, type="aep4")), nmom=6)
points(slmr$ratios[4], slmr$ratios[6])

# Exponential power distribution
# https://math.wm.edu/~leemis/chart/UDRF/PDFs/Exponentialpower.pdf
quapowexp <- function(f, para=NULL) {
  A <- 1/para$para[1]
  K <- 1/para$para[2]
  return( (A*log(1-log(1-f)))^K )
}
plotlmrdia(lmrdia())
for(i in 1:1000) {
  A <- runif(1, min=0, max=30); K <- runif(1, min=0, max=30)
  tlmr <- theoLmoms(para=list(type="powexp", para=c(1/A, 1/K)), quafunc=quapowexp)
  points(tlmr$ratios[3], tlmr$ratios[4], pch=16, cex=0.5)
}



T4s <- sort( c(seq(-0.24, 0.99, by=0.01), 0.999, 0, -0.005, +0.005) ) # strategic Tau4 values
T6s <- seq(0.1, 0.6, by=0.1)
L1 <- 0     # arithmetic mean of zero w/o loss of generality
L2 <- 1     # a unit of L-variation w/o loss of generality for location-scale families
T46df <- NULL # a data frame to construct for each Tau4 and each distribution
for(T4 in T4s) {
  for(T6 in T6s) {
    message(T4, appendLF=FALSE)
    ggg <- c(L1, L2, 0, T4, 0, T6)
    juha_np6 <- Lmoments::lmom2normpoly6(ggg)
    if(! any(is.na(juha_np6))) {
      juha_tlmr6 <- lmomco::theoTLmoms(juha_np6, quafunc=Lmoments::normpoly_inv, nmom=6)$ratios
      if(! is.na(juha_tlmr6[6])) {
        df <- data.frame(tau4=T4, tau6=T6, type="normpoly6",
                         tau4int=juha_tlmr6[4], tau6int=juha_tlmr6[6])
        T46df <- rbind(T46df, df) # data assembly
      }
    }
  }
}
