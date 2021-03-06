% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SigSetMethods.R
\docType{methods}
\name{oobR<-}
\alias{oobR<-}
\alias{oobR<-,SigSet-method}
\title{oobR replacement generic}
\usage{
oobR(x) <- value

\S4method{oobR}{SigSet}(x) <- value
}
\arguments{
\item{x}{object of \code{SigSet}}

\item{value}{new value}
}
\value{
a new \code{SigSet}
}
\description{
oobR replacement generic

Replace oobR slot of SigSet class
}
\examples{
sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
df <- oobR(sset)
df[1,1] <- 10
oobR(sset) <- df 
}
