"headrick.sheng.lalpha" <- function(x, bycovFF=FALSE, a=0.5, digits=8, ...) {

   pitems <- n <- NA
   names(pitems) <- "Number of items"
   names(n) <- "Sample size"
   zz <- list(alpha=NA, pitems=pitems, n=n, text="", source="headrick.sheng.alpha")

   if(! inherits(x, "matrix") && ! inherits(x, "data.frame")) {
     warning("'x' must be either a data.frame or a matrix, returning NULL")
     return(zz)
   }

   pitems <- ncol(x); names(pitems) <- "Number of items"
   zz$pitems <- pitems

   if(pitems < 2) {
     zz$text <- "less than two columns in x"
     return(zz)
   }
   ii <- ""
   for(i in seq_len(pitems)) {
     if(length(unique(x[,i])) == 1) {
       ii <- ifelse(ii == "", paste0(i, ""), paste0(ii, ",", i))
     }
     if(ii != "") {
       zz$text <- paste0("no variation in column(s) '", ii, "' degenerate condition")
     }
   }

   if(bycovFF) {
     x <- as.data.frame(x)
     n <- nrow(x); names(n) <- "Sample size"
     zz$n <- n
     lco2  <- matrix(nrow=pitems, ncol=pitems)
     for(i in seq_len(pitems)) {
       x[, i+pitems] <- pp(x[,i], sort=FALSE, a=a, ...)
       if(length(unique(x[,i])) == 1) {
         lco2[i, i] <- 0
       } else {
         lco2[i, i] <- lmoms(x[,i], nmom=2)$lambdas[2]
       }
     }
     x <- as.matrix(x)
     for(i in seq_len(pitems)) {
       for(j in  seq_len(pitems)) {
         if(i == j) next
         lco2[i,j] <- 2*stats::cov(x[,i], x[, j+pitems])
       }
     }
   } else {
     if(is.matrix(x)) {
       lco2 <- x
       zz$n <- "sample size is unknown, 2nd L-comoments were given"
     } else {
       lco2 <- Lcomoment.matrix(x, k=2)$matrix
       n <- nrow(x); names(n) <- "Sample size"
       zz$n <- n
     }
   }

   #print(lco2)

   dd <- dim(lco2)
   if(dd[1] != dd[2]) {
     warning("'lco2' not a square matrix")
     return(zz)
   }
   d <- dd[1]
   A <- sum(diag(as.matrix(lco2)))
   B <- sum(sapply(1:d, function(i) {
        sum(sapply(1:d, function(j) {
            ifelse(i == j, 0, lco2[i,j])
        } ))
   } ))
   alpha <- (d / (d-1)) * (1 - A / (A+B))
   if(alpha > 1) {
     alpha <- NA
     zz$text  <- paste0(zz$text, " and\n", "L-alpha computed as > 1, ",
                        "pathological issues in 'lco2', small sample size?")
   }
   names(alpha) <- "L-alpha"
   zz$alpha <- round(alpha, digits=digits)
   return(zz)
}


"lalpha" <- function(x, bycovFF=FALSE, a=0.5, digits=8, ...) {
   return(headrick.sheng.lalpha(x, bycovFF=bycovFF, a=a, digits=digits, ...))
}

