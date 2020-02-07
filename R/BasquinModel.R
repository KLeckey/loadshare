
#' Link function based on the Basquin Model
#' 
#' Default choice for the link function in all predictions. 
#'  
#' @param theta a vector of length two used as the parameter in the link function
#' @param x a vector of real numbers used as arguments in the link function.
#'
#' @return a vector with the same length as x.
#' 
#' @seealso \code{\link{getRates}}
#' 
#' @export  

h_basq<-function(theta,x){
  if(length(theta)!=2){
    stop("Dimension of theta has to be 2 for the Basquin model.")
  }
  return(exp(-theta[1])*(x^(theta[2])))
}

#' Gradient of the link function (Basquin Model)
#' 
#' The gradient of the function \code{\link{h_basq}} with respect to \code{theta}. 
#'  
#' @param theta a vector of length two used as the parameter in the link function
#' @param x a vector of real numbers used as arguments in the link function.
#'
#' @return a vector \code{y} in which \code{y[c(i,length(x)+i)]} represents the 
#' the gradient of \code{\link{h_basq}} with arguments \code{theta} and \code{x[i]}.
#' 
#' @seealso \code{\link{h_basq}}

grad_basq<-function(theta,x){
  if(length(theta)!=2){
    stop("Dimension of theta has to be 2 for the Basquin model.")
  }
  y1<- -exp(-theta[1])*x^(theta[2])
  y2<- log(x)*exp(-theta[1])*x^(theta[2])
  return(c(y1,y2))
}

#' Gradient of the rate function (Basquin Model)
#' 
#' The gradient of the \code{i}-th coordinate of the function 
#' \code{\link{getRates}} with respect to 
#' \code{theta} when using link function \code{\link{h_basq}}. 
#'  
#' @param theta,I,s parameters used in \code{\link{getRates}}. The missing 
#' variable \code{Iv} can be any value which is at least \code{i}.
#' @param i coordinate(s) of \code{\link{getRates}} which is/are considered. 
#'
#' @return a vector (or a matrix if \code{length(i)>1}) stating the values 
#' of the gradient at position \code{theta}
#' 
#' @seealso \code{\link{getRates},\link{h_basq}}


grad_lambda_basq <- function(theta, I, s, i){
  x <- s*I/(I-Iv)
  if(length(x)>1){
    N <- length(x)
    res1 <- rep(-1,N)*h_basq(theta,x)
    res2 <- log(x)*h_basq(theta,x)
    return(c(rbind(res1,res2)))
  }else{
    return(c(-1,log(x))*h_basq(theta,x))
  }
}

#' Prediction of information matrix (Basquin Model)
#' 
#' The simplified formula to predict the information matrix when the Basquin 
#' model is used.
#'  
#' @param s vector of stress ranges in the experiments
#' @param I total number of wires in each experiment
#' @param Iv vector of total number of broken wires observed in each experiment
#'
#' @return a matrix estimating the Fisher information
#' 
#' @seealso \code{\link{predImat}}

predImat_basq <- function(s,Iv,I){
  if(length(I)==1){
    I <- rep(I,length(s))
  }
  x_vec <- rep(s*I,Iv)/
    (rep(I,Iv)-getbreakvec(Iv))
  log_vec <- log(x_vec)
  
  return(matrix(c(sum(Iv),-sum(log_vec),-sum(log_vec), sum(log_vec^2)),nrow=2))
}

#' Coefficients in the cdf of the hypoexponential distribution
#' 
#' Coefficients \code{a_i(theta)} appearing in the cdf of the 
#' hypoexponential distribution.
#'  
#' @param theta,I,s parameter choices for the model
#' @param i index of the coefficient that is returned by the function
#' @param I0,Ic start and end in terms of broken wires for the rates
#'
#' @return the coefficient a_i(theta) of the cdf of the hypoexponential 
#' distribution

coeff_basq <- function(theta,i, I0,Ic, I, s){
  if(Ic-I0==1){
    return(1)
  }
  ### After simplification: ###
  indices <- I0:(Ic-1)
  indices <- indices[which(indices!=i)]
  factors<- 1/(1-((I-indices)/(I-i))^(theta[2]))
  return(prod(factors))
}

#' Gradient of the coefficients in the cdf of the hypoexponential distribution
#' 
#' Gradient with respect to \code{theta} of the coefficient given in 
#' \code{\link{coeff_basq}} 
#'  
#' @param theta value of theta where the gradient is evaluated
#' @param i,I0,Ic,I,s parameter choice for the coefficient in \code{\link{coeff_basq}} 
#'
#' @return value of the gradient at position \code{theta}

grad_coeff_basq <- function(theta,i, I0,Ic, I, s){
  if(Ic-I0==1){
    return(c(0,0))
  }
  indices <- I0:(Ic-1)
  indices <- indices[which(indices!=i)]
  fraction <- (I-indices)/(I-i)
  summands <- log(fraction)*(fraction)^(theta[2])/(1-fraction^(theta[2]))
  return(c(0, coeff_basq(theta,i,I0, Ic, I, s)*sum(summands)))
}


#' Coefficient vector in the cdf of the hypoexponential distribution
#' 
#' Vector of coefficients computed via \code{\link{coeff_basq}}.
#'  
#' @param theta,I0,Ic,I,s parameter choices in \code{\link{coeff_basq}}
#'
#' @return vector of values given by \code{\link{coeff_basq}} for \code{i=I0:Ic}

coeff_vec_basq <- function(theta,I0,Ic,I,s){
  N <- Ic-I0
  if(N==1){
    return(1)
  }
  res <- rep(0,N)
  for(i in I0:(Ic-1)){
    res[i-I0+1] <- coeff_basq(theta,i, I0,Ic,I,s)
  }
  return(res)
}

#' Derivative of the coefficient vector with respect to \code{theta[2]}
#' 
#' Derivative of \code{\link{coeff_vec_basq}} with respect to the second 
#' coordinate of \code{theta} evaluated at \code{theta}.
#'  
#' @param theta,i,I0,Ic,I,s parameter choice in \code{\link{coeff_vec_basq}}
#'
#' @return Derivative of \code{\link{coeff_vec_basq}} with respect to \code{theta[2]}

grad2_coeff_vec_basq <- function(theta, I0, Ic, I, s){
  N <- Ic-I0
  if(N==1){
    return(0)
  }
  res<- rep(0,N)
  for(i in I0:(Ic-1)){
    res[i-I0+1] <-  grad_coeff_basq(theta,i,I0,Ic,I,s)[2]
  }
  return(res)
}


#' Gradient of the cdf of the hypoexponential distribution
#' 
#' Gradient with respect to \code{theta} of the cdf of the 
#' hypoexponential distribution with link function \code{\link{h_basq}}.
#'  
#' @param theta value of theta where the gradient is evaluated
#' @param q quantile at which the cdf is evaluated
#' @param I0,Ic start and end in terms of broken wires for the rates
#' @param I,s parameter choice (max. number of wires and initial stress range)
#' in the unterlying experiment.
#'
#' @return gradient of the cdf with respect to \code{theta}

grad_phypoexp_basq <- function(theta, q , I0, Ic, I, s){
  
  indices <- I0:(Ic-1)
  grad_a_2_vec <- grad2_coeff_vec_basq(theta,I0,Ic,I,s)
  a_vec <- coeff_vec_basq(theta, I0, Ic, I, s)
  x_vec <- I*s/(I-indices)
  lambda_vec <- h_basq(theta, x_vec)
  
  res1 <- -q*dhypoexp2(q,lambda_vec)#First entry of gradient after simplification
  
  summands <- a_vec*log(x_vec)*lambda_vec*exp(-q*lambda_vec)*q-
    grad_a_2_vec*exp(-q*lambda_vec)
  res2 <- sum(summands)
  
  return(c(res1,res2))
}

#' Gradient of the quantile function of the hypoexponential distribution
#' 
#' Gradient with respect to \code{theta} of the quantile function of the 
#' hypoexponential distribution with link function \code{\link{h_basq}}.
#' The calculation of the gradient is based on the implicit function theorem.
#'  
#' @param theta value of theta where the gradient is evaluated
#' @param alpha p-value at which the quantile function is evaluated
#' @param I0,Ic start and end in terms of broken wires for the rates
#' @param I,s parameter choice (max. number of wires and initial stress range)
#' in the unterlying experiment.
#' 
#' @return Gradient of the quantile function with respect to \code{theta}

grad_qhypoexp_basq <- function(theta, alpha, I0, Ic, I, s){
  x_vec <- I*s/(I-I0:(Ic-1))
  lambda_vec <- h_basq(theta, x_vec)
  quant <- qhypoexp2(alpha, lambda_vec) 
  
  result <- -grad_phypoexp_basq(theta, quant, I0,Ic,I, s)/dhypoexp2(quant,lambda_vec)
  
  return(result)
}

