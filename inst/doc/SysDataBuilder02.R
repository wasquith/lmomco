setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Run SysDataBuilder01.R first when remaking the sysdata.rda.

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
load(file.choose()) # go find the sysdata.rda

# L1, first L-moment, arithmetic mean
# L2, second L-moment, L-variation
# T3, third L-moment ratio (L3/L2), L-skew
# T4, third L-moment ratio (L4/L2), L-kurtosis
# T5, third L-moment ratio (L5/L2), L-cinco (yes, some in literature, Mel Schaefer)
# T6, third L-moment ratio (L6/L2), named (L-seis ?)

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
                     lam2bak=blmr$lambdas[2], lam2int=tlmr$lambdas[2],
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
T6df <- T6df[ , c("type", "lam2bak", "lam2int", "special_parameter",
                  "tau4", "tau6", "tau4bak", "tau6bak", "tau4int", "tau6int")] # re-order
T6df$tau6err <- abs(T6df$tau6bak - T6df$tau6int) # compute an error for later study, if needed
T6df <- T6df[order(T6df$tau4), ]  # sort the table by L-kurtosis


# There might exist a need to have one canonical parent-like data object containing the relations
# between Tau4 and Tau6. So, let us begin assembling one with important distributions, cross-check
# against other sources when possible and otherwise tinker such that proper behavior is seen near
# Tau4 = 0 (uniform distribution) and (or) for cases of exactly L2 = 1.
# ***** NORMAL DISTRIBUTION *****
NOR <- lmomco::theoTLmoms(lmomco::vec2par(c(0, 1/sqrt(pi)), type="nor"), nmom=6, verbose=TRUE)
NOR$ratios[4] <- lmomco::par2lmom(lmomco::vec2par(c(0, 1/sqrt(pi)), type="nor"))$ratios[4]
tau46list <- list(nor=data.frame(tau4=NOR$ratios[4], tau6=NOR$ratios[6]))

# ***** TUKEY LAMBDA DISTRIBUTION ***** (This is not a constant L2 distribution, see lmomco::GLD)
# https://en.wikipedia.org/wiki/Tukey_lambda_distribution
Lambdas <- sort(c(0.001, -0.001, seq(-1, 1, by=0.005))) # see getting "close to zero" in the vector
TukeyL2 <- (2/Lambdas) * (- 1/(1+Lambdas) +  2/(2+Lambdas))
TukeyL4 <- (2/Lambdas) * (- 1/(1+Lambdas) + 12/(2+Lambdas) -  30/(3+Lambdas) +  20/(4+Lambdas))
TukeyL6 <- (2/Lambdas) * (- 1/(1+Lambdas) + 30/(2+Lambdas) - 210/(3+Lambdas) + 560/(4+Lambdas) -
                                                             630/(5+Lambdas) + 252/(6+Lambdas))
TukeyT4 <- TukeyL4 / TukeyL2
TukeyT6 <- TukeyL6 / TukeyL2

# Symmetrical L-moment Ratio Tables from Richard Vogel (March 25, 2024)
# We will see that the Stable is "new" per se to the lmomco infrastructure, the Student t with
# df matches the Student 3t from lmomco, the Tukey matches Wikipedia entry, and the Power
# Exponential is the AEP4 from lmomco with Tau3 = 0.
Stable  <- read.table("t4t6/StableDistribution.txt",   sep="|", header=TRUE)
Student <- read.table("t4t6/StudentT.txt",             sep="|", header=TRUE)
Tukey   <- read.table("t4t6/SymTukeyLambda.txt",       sep="|", header=TRUE)
PowExp  <- read.table("t4t6/PowerExponential.txt",     sep="|", header=TRUE)

plot(                                      T6df$tau6[T6df$type == "aep4"],
     approx(PowExp$tau4, PowExp$tau6, xout=T6df$tau4[T6df$type == "aep4"])$y)
abline(0, 1) # equal value line
tau46list$pwrexp <- data.frame(alpha=PowExp$alpha, tau4=PowExp$tau4, tau6=PowExp$tau6)
tau46list$aep4   <- data.frame(tau4=T6df$tau4[T6df$type == "aep4"], tau6=T6df$tau6[T6df$type == "aep4"])


plot(                                        T6df$tau6[T6df$type == "st3"],
     approx(Student$tau4, Student$tau6, xout=T6df$tau4[T6df$type == "st3"])$y)
abline(0, 1) # equal value line
tau46list$st2 <- data.frame(df=Student$df, tau4=Student$tau4, tau6=Student$tau6)
tau46list$st3 <- data.frame(tau4=Student$tau4, tau6=Student$tau6)


plot(                                    TukeyT6,
     approx(Tukey$tau4, Tukey$tau6, xout=TukeyT4)$y)
abline(0, 1) # equal value line
tau46list$tukeylam <- data.frame(lambda=Lambdas, lambda2=TukeyL2, tau4=TukeyT4, tau6=TukeyT6)
tau46list$tukeylam$tau4[1] <- tau46list$tukeylam$tau6[1] <- 1
tau46list$tukeylam$lambda2[tau46list$tukeylam$lambda == 0] <-
                  approx(   tau46list$tukeylam$lambda, y=tau46list$tukeylam$L2,   xout=0)$y
tau46list$tukeylam$lambda2[tau46list$tukeylam$lambda == 0] <- 1
tau46list$tukeylam$tau4[   tau46list$tukeylam$lambda == 0] <-
                  approx(   tau46list$tukeylam$lambda, y=tau46list$tukeylam$tau4, xout=0)$y
tau46list$tukeylam$tau6[   tau46list$tukeylam$lambda == 0] <-
                  approx(   tau46list$tukeylam$lambda, y=tau46list$tukeylam$tau6, xout=0)$y


tau46list$symstable <- Stable
j <- tau46list$symstable[0,]
j[1,] <- c(NA, 1, NA, NA)
tau46list$symstable <- rbind(tau46list$symstable, j)
tau46list$symstable <-   tau46list$symstable[order(tau46list$symstable$lambda2), ]
tau46list$symstable$tau4[tau46list$symstable$lambda2 == 1] <- # this is precise unity
                  approx( tau46list$symstable$lambda2, y=tau46list$symstable$tau4, xout=1)$y
tau46list$symstable$tau6[tau46list$symstable$lambda2 == 1] <- # this is precise unity
                  approx( tau46list$symstable$lambda2, y=tau46list$symstable$tau6, xout=1)$y

tau46list$gld_byt6tukeylam <- tau46list$tukeylam

tau46list$gld_byt5opt <- data.frame(khparas=T6df$special_parameter[T6df$type == "gld" ],
                                       tau4=T6df$tau4[T6df$type == "gld" ],
                                       tau6=T6df$tau6[T6df$type == "gld" ])
tau46list$pdq4 <- data.frame(kappa=T6df$special_parameter[T6df$type == "pdq4" ],
                              tau4=T6df$tau4[T6df$type == "pdq4"],
                              tau6=T6df$tau6[T6df$type == "pdq4"])


# We have now completed a list() object containing lookup tables of the Tau4 and Tau6 of
# symmetrical distributions for reference and creation (if needed externally) of T4-T6 diagrams.
# Note, the ST3 in lmomco has internally polynomial-based versions of the T4 and T6, see lmomst3().
####################################################################################################


assign("tau46list", tau46list, .lmomcohash)

save(.lmomcohash, file="sysdata.rda")

