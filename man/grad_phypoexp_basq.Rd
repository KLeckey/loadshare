% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasquinModel.R
\name{grad_phypoexp_basq}
\alias{grad_phypoexp_basq}
\title{Gradient of the cdf of the hypoexponential distribution}
\usage{
grad_phypoexp_basq(theta, q, I0, Ic, I, s)
}
\arguments{
\item{theta}{value of theta where the gradient is evaluated}

\item{q}{quantile at which the cdf is evaluated}

\item{I0, Ic}{start and end in terms of broken wires for the rates}

\item{I, s}{parameter choice (max. number of wires and initial stress range)
in the unterlying experiment.}
}
\value{
gradient of the cdf with respect to \code{theta}
}
\description{
Gradient with respect to \code{theta} of the cdf of the 
hypoexponential distribution with link function \code{\link{h_basq}}.
}
