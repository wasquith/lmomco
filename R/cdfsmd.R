"cdfsmd" <-
function(x, para) {
   if(! are.parsmd.valid(para)) return()
   U <- para$para[1] # location U >= 0
   A <- para$para[2] # scale A > 0
   B <- para$para[3] # shape B > 0
   Q <- para$para[4] # shape Q > 0
   # 0 <= x <= +Inf

   f <- 1 - ( 1 + ( (x - U) / A)^B )^(-Q)

   return(f)
}

