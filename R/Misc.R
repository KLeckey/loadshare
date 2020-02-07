
#' Function to check whether data is given in the correct format
#' 
#' An auxilliary function that produces an error if \code{data} is not given 
#' in the correct form for subsequent functions. 
#'  
#' @param data the data used for other function
#'
#' @return \code{invisible(NULL)}. Function produces an error if the data
#' has a wrong format
#' 
#' @export

DataSanityCheck <- function(data=list(w=NULL, s=NULL, Iv=NULL, I=35)){
  if(length(data$w)!=sum(data$Iv))
    stop("Mismatching data: length(data$w)!=sum(data$Iv)")
  if(length(data$Iv)!=length(data$s))
    stop("Mismatching data: length(data$s)!=length(data$Iv)")
  if( (length(data$I)!=1) & (length(data$I)!=length(data$s)))
    stop("Mismatching data: length(data$I) has to either be 1 or 
         equal to length(data$Iv)")
  
  return(invisible(NULL))
}



#' Auxilliary function 
#' 
#' An auxilliary function often used to compute a vector of rates in the 
#' Basquin model. The function expands a vector (a_1,...,a_N) of integers to
#' a longer vector (0,...,a_1-1,....,0,....,a_N-1).  
#'  
#' @param vec_a a vector of integers
#'
#' @return another vector \code{b} with a total length of \code{sum(vec_a)}. 
#' The vector \code{b} is a concatenation of the vectors 
#' \code{0:(vec_a[1]-1)},...,\code{0:(vec_a[N]-1)} with \code{N=length(vec_a)}
#' 
#' @seealso \code{\link{getRates}}

getbreakvec<-function(vec_a){
  k<-length(vec_a)
  n<-sum(vec_a)
  return(0:(n-1)-rep(c(0,cumsum(vec_a[-k])),vec_a))
}

#' Function to compute rates
#' 
#' A function to compute the rates of a hypoexponential distribution for a given
#' link function \code{linkfun}.
#'  
#' @param theta,s,Iv,I model parameters used in \code{linkfun}
#' @param linkfun link function used for the rate. Default choice is 
#' \code{\link{h_basq}}
#' 
#' @return vector of rates for the hypoexponential distribution.
#' 
#' @seealso \code{\link{h_basq},\link{getbreakvec}}
#' 
#' @export

getRates <- function(theta, s, Iv, I, linkfun=h_basq){
  if(length(I)==1){
    I <- rep(I, length(s))
  }
  x_vec <- rep(s*I, Iv)/(rep(I, Iv)-getbreakvec(Iv))
  rates <- linkfun(theta, x_vec)
  return(rates)
}

#' Inverse matrix with Cramer's rule
#' 
#' Function to compute the inverse of a 2x2 matrix using Cramer's rule
#'  
#' @param A a 2x2 matrix
#' 
#' @return Inverse of \code{A}
#' 
#' @export

inv_cramer<- function(A){
  if(!is.matrix(A) | any(dim(A)!=2)){
    stop("Only 2x2 matrices permitted.")
  }
  detA <- A[1,1]*A[2,2]-A[1,2]*A[2,1]
  if(detA!=0){
    return( matrix(c(A[2,2],-A[2,1],-A[1,2],A[1,1])/detA,nrow=2))
  }else{
    warning("Matrix has no inverse.")
    return(NULL)
  }
}
