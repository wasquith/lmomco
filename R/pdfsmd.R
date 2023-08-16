"pdfsmd" <-
function(x, para) {
   if(! are.parsmd.valid(para)) return()
   U <- para$para[1] # location U >= 0
   A <- para$para[2] # scale A > 0
   B <- para$para[3] # shape B > 0
   Q <- para$para[4] # shape Q > 0
   # 0 <= x <= +Inf

   #a <- B*Q*x^(B-1)
   #b <- A^B * (1 + ((x-U)/A)^B)^(Q+1)
   #f <- a / b
   suppressWarnings( a <-   log(B) + log(Q) + (B-1)*log(x - U) )
   #In log(x - U) : NaNs produced
   suppressWarnings( b <- B*log(A) + log((1 + ((x-U)/A)^B)^(Q+1)) )
   #In log((1 + ((x-U)/A)^B)^(Q+1)) : NaNs produced
   f <- exp( a - b )

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0
   return(f)
}
