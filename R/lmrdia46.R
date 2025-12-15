"lmrdia46" <- function(usrtrim=FALSE) {
   tau46list <- get("tau46list", envir=.lmomcohash)

   #cau <- matrix(nrow = 1, ncol = 2)
   #cau[1,] <- c(0.34280842, 0.20274358) # see lmomcau.Rd :: Examples
   cau <- data.frame(tau4=0.34280842, tau6=0.20274358)
   tau46list$cau <- cau

   #sla <- matrix(nrow = 1, ncol = 2)
   #sla[1,] <- c(0.30420472, 0.18900723) # see lmomsla.Rd :: Examples
   sla <- data.frame(tau4=0.30420472, tau6=0.18900723)
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
   tau46list$symgdd <- as.data.frame( gdd )

   if(usrtrim) {
     t4lim <- par()$usr[1:2]; t6lim <- par()$usr[3:4]
     for(d in names(tau46list)) {
       t4t6  <- tau46list[[d]]
       kt4   <- ncol(t4t6)-1; kt6 <- ncol(t4t6)
       t4s   <- seq(t4lim[1], t4lim[2], by=diff(t4lim)/1000)
       if(nrow(t4t6) < 2) {
         at4t6 <- as.matrix( data.frame(tau4=t4t6[,kt4], tau6=t4t6[,kt6]) )
       } else {
         df <- data.frame(tau4=t4s)
         df$tau6 <- approx(t4t6[,kt4], t4t6[,kt6], xout=t4s)$y
         if(ncol(t4t6) > 2) {
           for(nm in names(t4t6)[1:(kt4-1)]) {
             df[,nm] <- approx(t4t6[,kt4], t4t6[,nm], xout=t4s)$y
           }
         }
         at4t6 <- df[ , names(t4t6)]
       }
       at4t6 <- as.matrix( at4t6 )
       # Error in `[<-.data.frame`(`*tmp*`, !wnt, 2, value = NA) :
       #  missing values are not allowed in subscripted assignments of data frames
       wnt <- t4lim[1] <= at4t6[,kt4] & at4t6[,kt4] <= t4lim[2]
       at4t6[! wnt, kt4] <- NA
       wnt <- t6lim[1] <= at4t6[,kt6] & at4t6[,kt6] <= t6lim[2]
       at4t6[! wnt, kt6] <- NA
       tau46list[[d]] <- at4t6
     }
   }

   return(tau46list)
}
