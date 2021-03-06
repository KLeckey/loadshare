% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasquinModel.R
\name{grad_basq}
\alias{grad_basq}
\title{Gradient of the link function (Basquin Model)}
\usage{
grad_basq(theta, x)
}
\arguments{
\item{theta}{a vector of length two used as the parameter in the link function}

\item{x}{a vector of real numbers used as arguments in the link function.}
}
\value{
a vector \code{y} in which \code{y[c(i,length(x)+i)]} represents the 
the gradient of \code{\link{h_basq}} with arguments \code{theta} and \code{x[i]}.
}
\description{
The gradient of the function \code{\link{h_basq}} with respect to \code{theta}.
}
\seealso{
\code{\link{h_basq}}
}
