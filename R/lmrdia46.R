"lmrdia46" <- function() {
   tau46list <- get("tau46list", envir=.lmomcohash)

   cau <- matrix(nrow = 1, ncol = 2)
   cau[1,] <- c(0.34280842, 0.20274358) # see lmomcau.Rd :: Examples
   tau46list$cau <- cau

   sla <- matrix(nrow = 1, ncol = 2)
   sla[1,] <- c(0.30420472, 0.18900723) # see lmomsla.Rd :: Examples
   tau46list$sla <- sla

   # Gamma Difference Distribution
   tau4 <- c(0.1226017195408909, 0.1227, 0.123, 0.125, seq(0.13, 1, by=0.01))
   "gddsymt46f" <- function(t4) { # print(coefficients(LM))
     coe <- c( -0.0969112,    2.1743687, -12.8878580,  47.8931168, -108.0871549,
              156.9200440, -139.5599813,  69.3492358, -14.7052424)
     ix <- seq_len(length(coe))-1
     sapply(t4, function(t) sum(coe[ix+1]*t^ix))
   } # This function is from Note section in lmomgdd.R
   gdd <- matrix(nrow=length(tau4), ncol=2)
   gdd[,1] <- tau4; gdd[,2] <- gddsymt46f(tau4)
   colnames(gdd) <- c("tau4", "tau6")
   tau46list$symgdd <- gdd
   return(tau46list)
}
