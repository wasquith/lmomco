"cdfst3" <- function(x, para=NULL, paracheck=TRUE) {
   if(paracheck) {
      if(! are.parst3.valid(para)) return()
   }
   U <- para$para[1]
   A <- para$para[2]
   N <- para$para[3]

   SMALL.NU <- 1.001  # arrived from manual experiments involving theoLmoms() testing
   LARGE.NU <- 10^5.5 # arrived from manual experiments involving theoLmoms() testing

   if(N < SMALL.NU) N <- SMALL.NU
   if(N > LARGE.NU) N <- LARGE.NU

   f <- pt((x-U)/A, N)
   names(f) <- NULL
   return(f)
}

