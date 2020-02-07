#' Maximum likelihood estimation
#'
#' A function to compute the MLE for \code{theta} in model given by \code{linkfun}.
#'
#' @param data The experiment data given in a list with vectors \code{s,w,Iv,I}
#' representing the parameters in each experiment.
#' @param BasquinModel A boolean variable to indicate whether the Basquin model
#' is used. If \code{TRUE}, \code{linkfun} is automatically set to
#' \code{\link{h_basq}} and the more efficient function \code{\link{estML_basq}}
#' is used to compute the MLE.
#' @param linkfun the link function used in the model.
#' @param optim_start the parameter used as \code{par} in \code{\link{optim}}
#' @param thetaLower,thetaUpper restrictions to \code{theta} in the MLE
#' computation. The first coordinate of the estimator is restricted to be between
#' \code{thetaLower[1]} and \code{thetaUpper[1]} whereas the second coordinate
#' is contained in the interval (\code{thetaLower[1], thetaUpper[1]})
#' @param iterateOptim a boolean variable. If \code{TRUE} then \code{optim} is
#' used repeatedly until no improvement is made in the maximization process.
#' @param returnMLEonly a boolean variable. If \code{TRUE}, only the MLE
#' estimator for \code{theta} is returned. Otherwise, other infomation such
#' as the log-likelihood value at the MLE is returned as well.
#'
#' @return The MLE for \code{theta}
#'
#' @seealso \code{\link{estML_basq}}
#'
#' @export


estML<-function(data, BasquinModel=TRUE, linkfun=NULL,
                optim_start=c(28,3), thetaLower=c(-Inf,Inf), thetaUpper=c(Inf,Inf),
                iterateOptim =TRUE, returnMLEonly=TRUE, na.rm=TRUE){

  DataSanityCheck(data)

  #Use the more efficient estML_basq in the Basquin model
  if(BasquinModel){
    return(estML_basq(data, returnMLEonly))
  }

  if(is.null(linkfun)){
    stop("linkfun needs to be specified unless the default Basquin model is used.")
  }

  #Produce vector of inputs x for the link function
  x_vec<-rep(data$s*data$I,data$Iv)/
    (data$I-getbreakvec(data$Iv))

  #Log-likelihood function:
  logML<-function(theta){
    ## Penalize theta outside the allowed range with a -Inf log-likelihood
    if(any(theta<thetaLower) || any(theta>thetaUpper))
      return(-Inf)
    lam<-linkfun(theta,x_vec)
    return(sum(log(lam)-lam*data$w, na.rm = TRUE))
  }
  maxML<-optim(par=optim_start,logML,control=list(fnscale=-1))
  #In order to check for a bad numerical optim result, we test for improved results when optim_start is changed
  if(iterateOptim){
    improving<-TRUE
    while(improving){
      temp<-optim(par=maxML$par,logML,control=list(fnscale=-1))
      if(temp$value<maxML$value+0.01){
        improving<-FALSE
      }
      else{
        maxML<-temp
      }
    }
  }
  if(returnMLEonly){
    return(maxML$par)
  }
  else{
    return(maxML)
  }
}

#' Maximum likelihood estimation in the Basquin model
#'
#' A more efficient implementation for \code{\link{estML}} in the case
#' where \code{linkfun=h_basq}. Note that \code{thetaLower} and \code{thetaUpper}
#' are both chosen to be the default values (i.e. no restrictions to theta).
#'
#' @param data The experiment data given in a list with vectors \code{s,w,Iv,I}
#' representing the parameters in each experiment.
#' @param returnMLEonly a boolean variable. If \code{TRUE}, only the MLE
#' estimator for \code{theta} is returned. Otherwise, other infomation such
#' as the log-likelihood value at the MLE is returned as well.
#'
#' @return The MLE for \code{theta}
#'
#' @seealso \code{\link{estML}}

estML_basq <- function(data, returnMLEonly=TRUE, na.rm=TRUE){

  totalWires <- data$I
  N <- sum(data$Iv)

  if(length(totalWires)==1)
    totalWires <- rep(totalWires, length(data$s))

  x_vec <- rep(data$s*totalWires, data$Iv)/
    (rep(totalWires,data$Iv) - getbreakvec(data$Iv))

  theta1_ML <- function(theta2){
    return( log(sum(data$w*x_vec^theta2, na.rm=na.rm)/N))
  }
  deriv_theta1_ML <- function(theta2){
    return( sum(data$w*x_vec^theta2*log(x_vec), na.rm=na.rm)/sum(data$w*x_vec^theta2,na.rm=na.rm))
  }
  deriv_logML_theta2 <- function(theta2){
    return( sum( log(x_vec), na.rm=na.rm) - deriv_theta1_ML(theta2)*N)
  }


  theta2_hat <- uniroot(f = deriv_logML_theta2,interval= c(0,10),
                          extendInt = "yes")$root
  theta1_hat <- theta1_ML(theta2_hat)

  if(returnMLEonly){
    return(c(theta1_hat,theta2_hat))
  }
  else{
    logML<-function(theta){
      lam<-h_basq(theta,x_vec)
      return(sum(log(lam))-sum(lam*data$w))
    }
    return(list(par=c(theta1_hat,theta2_hat), value=logML(c(theta1_hat,theta2_hat))))
  }
}


