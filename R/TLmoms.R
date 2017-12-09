"TLmoms" <-
function(x, nmom=5, trim=NULL, leftrim=NULL, rightrim=NULL, vecit=FALSE) {
  if(nmom < 1) {
    warning("Number of L-moments is less than 1")
    return()
  }
  if(! is.null(trim) && trim < 0) {
    warning("Trimming value is less than 0")
    return()
  }
  if(! is.null(leftrim) && leftrim < 0) {
    warning("Left rimming value is less than 0")
    return()
  }
  if(! is.null(rightrim) && rightrim < 0) {
    warning("Right trimming value is less than 0")
    return()
  }
  if(is.null(trim) && is.null(leftrim) && is.null(rightrim)) {
    trim <- 0
  }

  # This traps for the entire sample being the same but this does leak if whole
  # sample is not but after trimming is!
  if(length(unique(x)) == 1) stop("all values are equal--TLmoments can not be computed")

  x <- sort(x)
  n <- length(x)

  if(nmom > n) {
    stop("More TLmoments requested by parameter 'nmom' than data points available in 'x'")
  }

  L <- R <- rep(NA, nmom)
  for(r in seq(1,nmom)) {
    lambda <- TLmom(x, trim=trim, leftrim=leftrim,
                    rightrim=rightrim, order=r, sortdata=FALSE)
    lr <- lambda$lambda
    L[r] <- ifelse(is.nan(lr), NA, lr)
  }

  if(nmom >= 2) R[2] <- L[2]/L[1]

  if(nmom >= 3) for(r in seq(3,nmom)) R[r] <- ifelse(!is.finite(L[r]/L[2]), NA, L[r]/L[2])

  z <- list(lambdas = L, ratios = R,
            trim=trim, leftrim=leftrim,
            rightrim=rightrim, source="TLmoms")
  if(! vecit) return(z)
  if(nmom == 1) {
     z <- z$lambdas[1]
  } else if(nmom == 2) {
     z <- c(z$lambdas[1], z$lambdas[2])
  } else {
     z <- c(z$lambdas[1], z$lambdas[2], z$ratios[3:nmom])
  }
  attr(z, which="trim")     <- trim
  attr(z, which="rightrim") <- rightrim
  attr(z, which="leftrim")  <- leftrim
  attr(z, which="source")   <- source
  return(z)
}
