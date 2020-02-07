
#' cdf of the hypoexponential distribution
#' 
#' An alternative calculation for the cdf of the hypoexponential distribution 
#' based on a matrix representation. A faster but numerically more instable
#' implementation for the cdf is given in \code{sdprisk::phypoexp}
#'  
#' @param x the quantile at which the cdf is evaluated
#' @param rates a vector of rates used in the hypoexponential distribution
#'
#' @return the p-value at \code{x}
#' 
#' @seealso \code{\link{qhypoexp2}, \link{dhypoexp2}}
#' 
#' @export  

phypoexp2<-function(x,rates){
  if(length(rates)==1){
    return(pexp(x, rate=rates))
  }
  l<-length(rates)
  D<-matrix(rep(0,l^2),nrow=l)
  D[col(D) - row(D) == 1] <- rates[1:(l-1)]
  D<-D+diag(-rates)
  eD<-Matrix::expm(D*x)
  return(1-sum(eD[1,]))
}

#' Quantile function of the hypoexponential distribution
#' 
#' An alternative calculation for the quantile function of the hypoexponential 
#' distribution based on a matrix representation. A faster but numerically 
#' more instable implementation for the cdf is given in \code{sdprisk::qhypoexp}
#'  
#' @param p the p-value at which the function is evaluated
#' @param rates a vector of rates used in the hypoexponential distribution
#' @param searchInterval an initial search interval for \code{\link{uniroot}}
#'
#' @return the quantile at \code{p}
#' 
#' @seealso \code{\link{phypoexp2}, \link{dhypoexp2}}
#' 
#' @export  

qhypoexp2<-function(p,rates, searchInterval=c(0,1e+10)){
  if(length(rates)==1){
    return(qexp(p, rate=rates))
  }
  F<-function(y){
    phypoexp2(y,rates)-p
  }
  result<-uniroot(f=F,interval=searchInterval,extendInt = "upX",tol=.Machine$double.eps^0.5)$root
  return(result)
}

#' Density function of the hypoexponential distribution
#' 
#' An alternative calculation for the density function of the hypoexponential 
#' distribution based on a matrix representation. A faster but numerically 
#' more instable implementation for the cdf is given in \code{sdprisk::dhypoexp}
#'  
#' @param x the position at which the density is evaluated
#' @param rates a vector of rates used in the hypoexponential distribution
#'
#' @return the value of the density at \code{x}
#' 
#' @seealso \code{\link{phypoexp2}, \link{qhypoexp2}}
#' 
#' @export

dhypoexp2<-function(x,rates){
  if(length(rates)==1){
    return(dexp(x, rate=rates))
  }
  l<-length(rates)
  D<-matrix(rep(0,l^2),nrow=l)
  D[col(D) - row(D) == 1] <- rates[1:(l-1)]
  D<-D+diag(-rates)
  eD<-D%*%Matrix::expm(D*x)
  return(-sum(eD[1,]))
}

#' Function to check for numeric instabilities
#' 
#' A boolean function that checks whether the fast implementation in
#' \code{sdprisk::phypoexp} has severe numerical instabilities. To this end,
#' the function evaluates the sum of all coefficients in the cdf, which should 
#' add up to one. If the sum deviates from one by more than the prespecified 
#' \code{precision} value, the function returns \code{FALSE} and 
#' \code{sdprisk::phypoexp} should not be used. Note that even in the case
#' where the return value is \code{TRUE}, other kinds of numeric instabilites 
#' can occur. In situations with numeric instabilities the more stable 
#' implementation in \code{\link{phypoexp2}} should be used.
#'  
#' @param rates a vector of rates used in the hypoexponential distribution
#'
#' @return a boolean value to indicate whether \code{sdprisk::phypoexp} seems 
#' numerically stable
#' 
#' @seealso \code{\link{phypoexp2}}
#' 
#' @export

checkHypoexp <- function(rates, precision=0.001){
  coeff <- function(i){
    return( prod( 1/ (1- rates[i]/rates[-i] )))
  }
  temp <- sapply(1:length(rates),FUN=coeff)
  return(abs(sum(temp)-1) <= precision)
}

#' Prediction intervals for the hypoexponential distribution
#' 
#' A function that returns a \code{(1-alpha)}-prediction interval for a 
#' hypoexponential distribution with given \code{rates}. The prediction 
#' interval is given as a vector containing the \code{alpha/2}-quantile
#' and the \code{(1-alpha/2)}-quantile of the distribution.
#' 
#' @param alpha the alpha-error in the prediction 
#' @param rates a vector of rates used in the hypoexponential distribution
#' @param use_matrix_rep a boolean variable to choose whether 
#' \code{\link{qhypoexp2}} or \code{sdprisk::qhypoexp} should be used to compute
#' quantiles. The variable is overwritten by an automatic choice if 
#' \code{auto_matrix_rep=TRUE}
#' @param auto_matrix_rep a boolean variable. If \code{TRUE}, the choice 
#' for \code{use_matrix_rep} is done automatically by checking for numeric 
#' instabilities via \code{\link{checkHypoexp}} 
#'
#' @return a \code{(1-alpha)}-prediction interval
#' 
#' @seealso \code{\link{qhypoexp2}}
#' 
#' @export

predHypoexp<-function(alpha,rates, use_matrix_rep=TRUE, auto_matrix_rep=TRUE){
  if(alpha<=0|alpha>=1){
    stop("alpha has to be between zero and one.")
  }
  if(length(rates)==1){
    return(c(qexp(alpha/2,rate=rates),qexp(1-alpha/2,rate=rates)))
  }
  if(auto_matrix_rep){
    #Check for numeric instabilities in qhypoexp
    use_matrix_rep <- !checkHypoexp(rates)
  }
  if(use_matrix_rep){
    return(c(qhypoexp2(alpha/2,rates), qhypoexp2(1-alpha/2,rates)))
  }
  else{
    return(c(sdprisk::qhypoexp(alpha/2,rates), sdprisk::qhypoexp(1-alpha/2,rates)))
  }
}