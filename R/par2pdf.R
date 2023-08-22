"par2pdf" <-
function(x,para,...) {
    type <- para$type
    if(type == 'aep4') {
      return(pdfaep4(x,para))
    }
    else if(type == 'cau') {
      return(pdfcau(x,para))
    }
    else if(type == 'emu') {
      return(pdfemu(x,para))
    }
    else if(type == 'exp') {
      return(pdfexp(x,para))
    }
    else if(type == 'gam') {
      return(pdfgam(x,para))
    }
    else if(type == 'gep') {
      return(pdfgep(x,para))
    }
    else if(type == 'gev') {
      return(pdfgev(x,para))
    }
    else if(type == 'gld') {
      return(pdfgld(x,para,...))
    }
    else if(type == 'glo') {
      return(pdfglo(x,para))
    }
    else if(type == 'gno') {
      return(pdfgno(x,para))
    }
    else if(type == 'gov') {
      return(pdfgov(x,para))
    }
    else if(type == 'gpa') {
      return(pdfgpa(x,para))
    }
    else if(type == 'gum') {
      return(pdfgum(x,para))
    }
    else if(type == 'kap') {
      return(pdfkap(x,para))
    }
    else if(type == 'kmu') {
      return(pdfkmu(x,para))
    }
    else if(type == 'kur') {
      return(pdfkur(x,para))
    }
    else if(type == 'lap') {
      return(pdflap(x,para))
    }
    else if(type == 'lmrq') {
      return(pdflmrq(x,para))
    }
    else if(type == 'ln3') {
      return(pdfln3(x,para))
    }
    else if(type == 'nor') {
      return(pdfnor(x,para))
    }
    else if(type == 'pe3') {
      return(pdfpe3(x,para))
    }
    else if(type == 'pdq3') {
      return(pdfpdq3(x,para,...))
    }
    else if(type == 'pdq4') {
      return(pdfpdq4(x,para,...))
    }
    else if(type == 'ray') {
      return(pdfray(x,para))
    }
    else if(type == 'revgum') {
      return(pdfrevgum(x,para))
    }
    else if(type == 'rice') {
      return(pdfrice(x,para))
    }
    else if(type == 'sla') {
      return(pdfsla(x,para))
    }
    else if(type == 'smd') {
      return(pdfsmd(x,para))
    }
    else if(type == 'st3') {
      return(pdfst3(x,para))
    }
    else if(type == 'texp') {
      return(pdftexp(x,para))
    }
    else if(type == 'tri') {
      return(pdftri(x,para))
    }
    else if(type == 'wak') {
      return(pdfwak(x,para))
    }
    else if(type == 'wei') {
      return(pdfwei(x,para))
    }
    else {
      stop("Did not find a valid distribution type")
    }
}

