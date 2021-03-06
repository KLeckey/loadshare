% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HypoExp.R
\name{phypoexp2}
\alias{phypoexp2}
\title{cdf of the hypoexponential distribution}
\usage{
phypoexp2(x, rates)
}
\arguments{
\item{x}{the quantile at which the cdf is evaluated}

\item{rates}{a vector of rates used in the hypoexponential distribution}
}
\value{
the p-value at \code{x}
}
\description{
An alternative calculation for the cdf of the hypoexponential distribution 
based on a matrix representation. A faster but numerically more instable
implementation for the cdf is given in \code{sdprisk::phypoexp}
}
\seealso{
\code{\link{qhypoexp2}, \link{dhypoexp2}}
}
