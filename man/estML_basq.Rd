% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLEstimation.R
\name{estML_basq}
\alias{estML_basq}
\title{Maximum likelihood estimation in the Basquin model}
\usage{
estML_basq(data, returnMLEonly = TRUE, na.rm = TRUE)
}
\arguments{
\item{data}{The experiment data given in a list with vectors \code{s,w,Iv,I}
representing the parameters in each experiment.}

\item{returnMLEonly}{a boolean variable. If \code{TRUE}, only the MLE
estimator for \code{theta} is returned. Otherwise, other infomation such
as the log-likelihood value at the MLE is returned as well.}
}
\value{
The MLE for \code{theta}
}
\description{
A more efficient implementation for \code{\link{estML}} in the case
where \code{linkfun=h_basq}. Note that \code{thetaLower} and \code{thetaUpper}
are both chosen to be the default values (i.e. no restrictions to theta).
}
\seealso{
\code{\link{estML}}
}
