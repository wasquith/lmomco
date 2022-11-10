\name{dist.list}
\alias{dist.list}
\title{List of Distribution Names}
\description{
Return a list of the three character syntax identifying distributions supported within the \pkg{lmomco} package. The distributions are
\code{aep4}, \code{cau}, \code{emu}, \code{exp}, \code{gam}, \code{gep},
\code{gev}, \code{gld}, \code{glo}, \code{gno}, \code{gov},
\code{gpa}, \code{gum}, \code{kap}, \code{kmu}, \code{kur},
\code{lap}, \code{lmrq}, \code{ln3}, \code{nor}, \code{pdq3}, \code{pdq4},
\code{pe3}, \code{ray}, \code{revgum}, \code{rice}, \code{sla}, \code{st3},
\code{texp}, \code{tri}, \code{wak}, and \code{wei}. These abbreviations and only these are used in routing logic within \pkg{lmomco}. There is no provision for fuzzy matching. The full distributions names are available in \code{\link{prettydist}}.
}
\usage{
dist.list(type=NULL)
}
\arguments{
\item{type}{If \code{type} is not \code{NULL} and is one of the abbreviations shown above, then the number of parameters of that distribution are returned or a warning message is issued. This subtle feature might be useful for developers.}
}
\value{
A vector of distribution identifiers as listed above or the number of parameters for a given distribution type.
}
\author{W.H. Asquith}
\seealso{\code{\link{prettydist}}
}
\examples{
dist.list("gpa")

\dontrun{
# Build an L-moment object
LM <- vec2lmom(c(10000, 1500, 0.3, 0.1, 0.04))
lm2 <- lmorph(LM)  # convert to vectored format
lm1 <- lmorph(lm2) # and back to named format
dist <- dist.list()
# Demonstrate that lmom2par internally converts to needed L-moment object
for(i in 1:length(dist)) {
  # Skip Cauchy and Slash (need TL-moments).
  # Skip AEP4, Kumaraswamy, LMRQ, Student t (3-parameter), Truncated Exponential
  # are skipped because each is inapplicable to the given L-moments.
  # The Eta-Mu and Kappa-Mu are skipped for speed.
  if(dist[i] == 'aep4' | dist[i] == 'cau' | dist[i] == 'emu'  | dist[i] == 'gep' |
     dist[i] == 'kmu'  | dist[i] == 'kur' | dist[i] == 'lmrq' | dist[i] == 'tri' |
     dist[i] == 'sla'  | dist[i] == 'st3' | dist[i] == 'texp') next
  message(dist[i], " parameters : ",
          paste(round(lmom2par(lm1, type=dist[i])$para, digits=4), collapse=", "))
  message(dist[i], " parameters : ",
          paste(round(lmom2par(lm2, type=dist[i])$para, digits=4), collapse=", "))
} # }
}
\keyword{distribution (utility)}
