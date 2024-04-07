"quast3" <-
function(f, para, paracheck=TRUE) {
   if(! check.fs(f)) return()
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

   x <- U + A*qt(f, N)
   names(x) <- NULL
   return(x)

   stop("Should not be here in execution")
}

