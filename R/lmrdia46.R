"lmrdia46" <- function() {
   tau46list <- get("tau46list", envir=.lmomcohash)

   cau <- matrix(nrow = 1, ncol = 2)
   cau[1,] <- c(0.34280842, 0.20274358) # see lmomcau.Rd :: Examples
   tau46list$cau <- cau

   sla <- matrix(nrow = 1, ncol = 2)
   sla[1,] <- c(0.30420472, 0.18900723) # see lmomsla.Rd :: Examples
   tau46list$sla <- sla

   return(tau46list)
}
