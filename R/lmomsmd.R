"lmomsmd" <-
function(para) {
  if(! are.parsmd.valid(para)) {
    warning("Parameters are invalid")
    return()
  }

  #print(par)
  names(para$para) <- NULL
  U <- para$para[1] # location U >= 0
  A <- para$para[2] # scale A > 0
  B <- para$para[3] # shape B > 0
  Q <- para$para[4] # shape Q > 0
  # 0 <= x <= +Inf

  IB <- 1/B
  BR <- exp( log(A) + lgamma(1+IB) )
  t1 <- exp( lgamma(1*Q - IB) - lgamma(1*Q)  )
  t2 <- exp( lgamma(2*Q - IB) - lgamma(2*Q)  )
  t3 <- exp( lgamma(3*Q - IB) - lgamma(3*Q)  )
  t4 <- exp( lgamma(4*Q - IB) - lgamma(4*Q)  )
  t5 <- exp( lgamma(5*Q - IB) - lgamma(5*Q)  )
  t6 <- exp( lgamma(6*Q - IB) - lgamma(6*Q)  )

  L1 <- BR * t1 + U
  L2 <- BR * (1*t1 -  1*t2)
  L3 <- BR * (1*t1 -  3*t2 +  2*t3)
  L4 <- BR * (1*t1 -  6*t2 + 10*t3 -   5*t4)
  L5 <- BR * (1*t1 - 10*t2 + 30*t3 -  35*t4 +  14*t5)
  L6 <- BR * (1*t1 - 15*t2 + 70*t3 - 140*t4 + 126*t5 - 42*t6)

  # 2 #  -1    0    0    0   0
  # 3 #  -3   +2    0    0   0
  # 4 #  -6  +10   -5    0   0
  # 5 # -10  +30  -35  +14   0
  # 6 # -15  +70 -140 +126 -42

  z <- list(lambdas=c(L1, L2,    L3,    L4,    L5,    L6),
            ratios= c(NA, L2/L1, L3/L2, L4/L2, L5/L2, L6/L2),
            trim=0, leftrim=0, rightrim=0,
            source="lmomsmd")
  return(z)
}
