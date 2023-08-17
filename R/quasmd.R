"quasmd" <-
function(f, para, paracheck=TRUE) {
   if(! check.fs(f)) return()
   if(paracheck == TRUE) {
      if(! are.parsmd.valid(para)) return()
   }
   U <- para$para[1] # location U >= 0
   A <- para$para[2] # scale A > 0
   B <- para$para[3] # shape B > 0
   Q <- para$para[4] # shape Q > 0
   # 0 <= x <= +Inf

   x <- U + A * ( (1-f)^(-1/Q) - 1 )^(1/B)
   names(x) <- NULL
   return(x)
}

