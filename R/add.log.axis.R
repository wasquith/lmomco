"add.log.axis" <-
function(make.labs=FALSE, logs=c(2, 3, 4, 5, 6, 7, 8, 9), side=1,
         two.sided=FALSE, label=NULL, x=NULL, col.ticks=1, ...) {
   if(! is.null(x)) { # aligning or snapping limits in untransformed x
     x <- log10(x[! is.na(x) & x > 0])
     y <- abs(as.integer(min(x))) + 1
     lim <- 10^(range(x+y)) # complicating sign reversals at log10(1)
     # so lets offset based on orders of magnitude below zero and offset later
     needs <- log10(1:10)
     wants <- log10(lim) - as.integer(log10(lim))
     left  <- ifelse(wants[1] == 0, 0, rev(needs[needs <= wants[1]])[1])
     right <- ifelse(wants[2] == 0, 0,     needs[needs >= wants[2]][1])
     return(10^((as.integer(log10(lim)) + c(left,right))-y))
   }
   ifelse(side == 2 | side == 4, interval <- par()$usr[3:4],
                                 interval <- par()$usr[1:2])
   interval <- 10^interval

   other.side <- switch(as.character(side), "1"=3, "2"=4, "3"=1, "4"=1)
   interval <- sort(interval)
   min <- interval[1];    max <- interval[2]
   log.beg  <- as.integer(log10(min)) - 1
   log.end  <- as.integer(log10(max)) + 1
   the.logs <- vector(mode="numeric")
   for(cycle in log.beg:log.end) {
      the.logs.in.cycle <- logs*10^(cycle)
      the.logs <- c(the.logs, the.logs.in.cycle)
   }
   if(make.labs) {
      axis(side, at=the.logs, labels=the.logs, tcl=0, col.ticks=col.ticks,  col=NA, ...)
      mtext(label, line=2, side=side)
   } else {
      axis(side, at=the.logs, labels=NA, col.ticks=col.ticks,  col=NA, ...)
      if(two.sided) axis(other.side, at=the.logs,  labels=NA, col.ticks=col.ticks,  col=NA, ...)
   }
}

