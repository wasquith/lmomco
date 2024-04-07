"pdfst3" <-
function(x, para, paracheck=TRUE) {
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

   f <- dt((x-U)/A, N)/A
   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}
