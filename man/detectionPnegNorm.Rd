% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detection.R
\name{detectionPnegNorm}
\alias{detectionPnegNorm}
\title{Detection P-value based on normal fitting the negative controls}
\usage{
detectionPnegNorm(sset)
}
\arguments{
\item{sset}{a \code{SigSet}}
}
\value{
detection p-value
}
\description{
The function takes a \code{SigSet} as input, computes detection p-value
using negative control probes parametrized in a normal distribution and
returns a new \code{SigSet} with an updated pval slot.
}
\details{
Background of Grn and Red are estimated separately from negative control
probes-parameterized normal distribution. p-value is taken from the
minimum of the p-value of the two alleles (color depends on probe design).
}
\examples{
sset <- makeExampleSeSAMeDataSet()
sset <- detectionPnegNorm(sset)

}
