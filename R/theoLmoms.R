"theoLmoms" <-
function(para, nmom=5, minF=0, maxF=1, quafunc=NULL,
               nsim=50000, fold=5,
               silent=TRUE, verbose=FALSE, ...) {
  if(nmom < 1) {
    warning("Number of L-moments requested is less than 1")
    return()
  }
  zz <- theoTLmoms(para, nmom=nmom, trim=0, minF=minF, maxF=maxF, quafunc=quafunc,
                         nsim=nsim, fold=fold, silent=silent, verbose=verbose, ...)
  zz$source="theoLmoms"
  return(zz)
}
