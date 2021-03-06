% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConSets.R
\name{ConSet2d_wald}
\alias{ConSet2d_wald}
\title{Computation of confidence sets based on Wald's test}
\usage{
ConSet2d_wald(data, alpha, BasquinModel, linkfun, gradlinkfun, precision,
  thetahat)
}
\arguments{
\item{data, alpha, BasquinModel, linkfun, gradlinkfun}{parameters chosen in \code{\link{ConSet}}}

\item{precision}{the precision in the representation of the ellipse. This value
corresponds to the \code{mesh_width} in \code{\link{ConSet}}}

\item{thetahat}{the maximum likelihood estimation (MLE) for the parameter
theta. The MLE can be found via \code{\link{estML}}.}
}
\value{
a matrix with three rows and a variable number of columns.
Details can be found in \code{\link{ConSet}}
}
\description{
An alternative implementation for confidence sets avoiding grid search.
This faster implementation makes use of the fact that the confidence set
based on Wald's test is an ellipse.
}
\seealso{
\code{\link{ConSet}}
}
