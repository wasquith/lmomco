setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Change to the directory containing the top-level of the package
# and then source this file to source all of the R files.
# This utility is strictly for the core developer to test things
# without ever having to actually install the package!
for(f in list.files("./R", pattern="\\.R")) source(paste0("./R/", f))


  lmrsmd <- lmomsmd(vec2par(c(0, 1000, 10, 100), type="smd"))
  para <- parsmd( lmrsmd, emplims=FALSE )



MU <- 1
t3s <- sort(unique(c(-0.18, -0.178, -0.176, -0.174, -0.172,
                 seq(-0.17, +1, by=0.005), 0.999)))
t4s <- sort(unique(seq(0.1, +1, by=0.0005) ))
eT3s <- eT4s <- tT3s <- tT4s <- Bs <- Qs <- Iters <- NULL
plot(c(-0.25,1), c(0.1,1), type="n")
for(t3 in t3s) {
  for(t4 in t4s) {
    if(abs(t3) == 1) next
    if(abs(t4) == 1) next
    lmr <- vec2lmom(c(1000, 10, t3, t4), lscale=TRUE, checklmom=FALSE)
    if(! are.lmom.valid(lmr)) next
    para <- parsmd(lmr, emplims=FALSE)
    lmrsmd <- lmomsmd(para)
    if(! are.lmom.valid(lmrsmd)) next
    #if(! is.finite(lmrsmd$lambdas[1])) next
    #errt1a <- abs( MU - lmrsmd$lambdas[1] )
    #errt2a <- abs( t2 - lmrsmd$ratios[2]  )
    #errt3a <- abs( t3 - lmrsmd$ratios[3]  )
    #if(! is.finite(lmrsmd$ratios[ 4]))      next
    #if(        abs(lmrsmd$ratios[ 4]) >= 1) next
    #if(errt1a > 0.0001) next
    #if(errt2a > 0.0001) next
    #if(errt3a > 0.0001) next
    if(para$para[3]/para$para[4] > 1E+14) next
    if(para$para[3]/para$para[4] < 1E-14) next
    if(para$para[4]/para$para[3] > 1E+12) next
    if(para$para[4]/para$para[3] < 1E-12) next
    #if(para$iter != 1) {
      points(lmrsmd$ratios[3], lmrsmd$ratios[4], cex=0.4, pch=16, col=para$iter)
    #}
    if(para$ifail) {
      points(lmrsmd$ratios[3], lmrsmd$ratios[4], cex=0.2, pch=16, col="seagreen")
    } else {
      Bs <- c(Bs, para$para[3]); Qs <- c(Qs, para$para[4])
      tT3s <- c(tT3s, t3); tT4s <- c(tT4s, t4); Iters <- c(Iters, para$iter)
      eT3s <- c(eT3s, lmrsmd$ratios[3]); eT4s <- c(eT4s, lmrsmd$ratios[4])
    }
  }
}
lmrdia <- lmomco::lmrdia()
lines(lmrdia$gpa[,1], lmrdia$gpa[,2], col="blue")


SMD <- data.frame(TAU3=tT3s, TAU4=tT4s, B=Bs, Q=Qs, Iters=Iters)
save(SMD, file="SMD.RData")

plot(tT3s, tT4s, cex=0.4, pch=16, col=Iters)

tmp <- SMD[SMD$Iters == 1,]
x <- c(-0.188546284, -0.072683616, -0.047221860, -0.026665580,  0.006271186)
y <- c(0.1500420, 0.1158318, 0.1103633, 0.1078197, 0.1034958)
tmp <- tmp[tmp$TAU4 > approx(x, y, xout=tmp$TAU3, rule=2)$y,]
x <- c(-0.12033681, -0.12033681, -0.11029226, -0.03063668,
       0.1499298, -0.1679737, -0.1679737, -0.1679737, -0.1682314, -0.1700358,
       -0.1700358, -0.1677159, -0.1700358, -0.1679737, -0.1617872, -0.1195131)
y <- c(0.1782749, 0.1720434, 0.1611063, 0.1681009,
       0.1538342, 0.1538342, 0.1548486, 0.1560658, 0.1570802, 0.1538342,
       0.1461252, 0.1461252, 0.1449080, 0.1438936, 0.1438936, 0.163572)
for(i in 1:length(x)) {
  tmp <- tmp[sqrt((x[i]-tmp$TAU3)^2+(y[i]-tmp$TAU4)^2) > 0.001, ]
}
plot(tmp$TAU3, tmp$TAU4, cex=0.2, pch=16)
SMDmax <- aggregate(tmp, by=list(tmp$TAU3), max)
SMDmin <- aggregate(tmp, by=list(tmp$TAU3), min)
#lines(SMDmax$Group.1, SMDmax$TAU4, col="red")
#lines(SMDmin$Group.1, SMDmin$TAU4, col="red")

UPR <- lm(SMDmax$TAU4~SMDmax$Group.1+I(SMDmax$Group.1^2)+I(SMDmax$Group.1^3)+
                    I(SMDmax$Group.1^4)+I(SMDmax$Group.1^5)+I(SMDmax$Group.1^6)+
                    I(SMDmax$Group.1^7)+I(SMDmax$Group.1^8))
h <- coefficients(UPR)
names(h) <- NULL
message("UPR:")
print(h, digits=8)
lines(SMDmax$Group.1, fitted.values(UPR), col="red", lwd=2)

LWR <- lm(SMDmin$TAU4~SMDmin$Group.1+I(SMDmin$Group.1^2)+I(SMDmin$Group.1^3)+
                    I(SMDmin$Group.1^4)+I(SMDmin$Group.1^5)+I(SMDmin$Group.1^6)+
                    I(SMDmin$Group.1^7)+I(SMDmin$Group.1^8))
h <- coefficients(LWR)
names(h) <- NULL
message("LWR:")
print(h, digits=8)
lines(SMDmin$Group.1, fitted.values(LWR), col="red", lwd=2)

plot(tmp$TAU3, tmp$TAU4, cex=0.4, pch=16, xlim=c(-.2,0), ylim=c(0.1, 0.2))
lines(SMDmax$Group.1, fitted.values(UPR), col="red", lwd=2)
lines(SMDmin$Group.1, fitted.values(LWR), col="red", lwd=2)

plot(tmp$TAU3, tmp$TAU4, cex=0.4, pch=16, xlim=c(.9,1), ylim=c(0.9, 1))
lines(SMDmax$Group.1, fitted.values(UPR), col="red", lwd=2)
lines(SMDmin$Group.1, fitted.values(LWR), col="red", lwd=2)


stop()
