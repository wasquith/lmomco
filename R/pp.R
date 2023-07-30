"pp" <-
function(x, A=NULL, B=NULL, a=0, sort=TRUE, ties.method="first", ...) {

   if(! is.null(a)) {
      if(a < 0 | a > 0.50) {
         warning("Plotting position parameter a is invalid, not in [0,0.5]")
         return()
      }
      A <- -a
      B <- 1 - 2*a
   } else {
     if(is.null(A)) {
       warnings("Plotting position parameter A is NULL")
       return(NULL)
     }
     if(is.null(B)) {
       warnings("Plotting position parameter B is NULL")
       return(NULL)
     }
     if(A <= -1) {
       warnings("Plotting position parameter A <= -1, invalid")
       return(NULL)
     }
     if(B <= -1) {
       warnings("Plotting position parameter B <= -1, invalid")
       return(NULL)
     }
     if(A >= B) {
       warnings("Plotting position parameters B < A, invalid")
       return(NULL)
     }
   }

   denom <- length(x) + B
   ranks <- rank(x, ties.method=ties.method)

   if(sort) {
      return( (sort(ranks) + A) / denom)
   } else {
      return(      (ranks  + A) / denom)
   }
}

# plot(c(-1,1), c(-1,1), type="n")
# for(A in seq(-1,1,by=.02)) {
#   for(B in seq(-1,1, by=.02)) {
#     the.pp <- pp(1:10, A=A, B=B, a=NULL)
#     if(is.null(the.pp)) next
#     suppressWarnings(col <- as.numeric(check.fs(the.pp)))
#     points(A,B, col=col, pch=16)
#   }
# }
