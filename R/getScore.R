
#' Coverage rates of prediction intervals
#' 
#' A computation of the coverage rate of prediction intervals given as a
#' matrix \code{predInts}. The values that need to be predicted are given as
#' a vector \code{jumptimes}
#'  
#' @param jumptimes a vector of length N storing the times that are to be 
#' predicted
#' @param predInts a N by 2 matrix in which \code{predInts[n,]} denotes the
#' prediction interval for \code{jumptimes[n]} 
#'
#' @return a number between zero and one that represents the relative number
#' of points in \code{jumptimes} that are covered by the corresponding 
#' prediction interval
#' 
#' @export  

getCovRate <- function(jumptimes, predInts){
  if(is.matrix(predInts)){
    cov <- sum( (predInts[,1]<= jumptimes)*(predInts[,2]>= jumptimes), na.rm = TRUE)
    return(cov/length(jumptimes))
  }
  else if(length(jumptimes)==1 && length(predInts)==2){
    cov <- (predInts[1]<= jumptimes)*(predInts[2]>= jumptimes)
    return(cov)
  }
  else{
    stop("predInts needs to be an n by 2 matrix if length(jumptimes)=n")
  }
}


#' Interval lengths of prediction intervals
#' 
#' A computation of the average interval length of prediction intervals given 
#' as a matrix \code{predInts}. 
#'  
#' @param predInts a N by 2 matrix in which \code{predInts[n,]} denotes the
#' prediction interval for \code{jumptimes[n]} 
#'
#' @return the average length (mean) of the intervals in \code{predInts}
#' 
#' @export  

getIntLength <- function(predInts){
  if(is.matrix(predInts)){
    return(sum(predInts[,2]-predInts[,1], na.rm = TRUE)/length(predInts[,1]))
  }
  else{
    return(predInts[2]-predInts[1])
  }
  
}

#' Interval score of prediction intervals
#' 
#' A computation of the interval score by Gneiting and Raftery (2007) of 
#' prediction intervals given as a
#' matrix \code{predInts}. The values that need to be predicted are given as
#' a vector \code{jumptimes}. 
#' 
#' @param jumptimes a vector of length N storing the times that are to be 
#' predicted
#' @param predInts a N by 2 matrix in which \code{predInts[n,]} denotes the
#' prediction interval for \code{jumptimes[n]} 
#' @param alpha the alpha-error chosen for the prediction intervals
#'
#' @return the average interval score of the prediction intervals
#' 
#' @export

getScore <- function(jumptimes, predInts, alpha){
  
  if(is.matrix(predInts)){
    temp <- sum((predInts[,1]-jumptimes)*(predInts[,1]>jumptimes), na.rm=TRUE)+
            sum((jumptimes - predInts[,2])*(predInts[,2]<jumptimes), na.rm=TRUE)+
            sum(jumptimes*is.na(predInts[,1])) #Added: Penality for empty prediction intervals
    temp <- temp/length(jumptimes)
    
    return(getIntLength(predInts)+2/alpha * temp)
  }
  else{
    temp <- (predInts[1]-jumptimes)*(predInts[1]>jumptimes)+
               (jumptimes - predInts[2])*(predInts[2]<jumptimes)
    
    return(getIntLength(predInts)+2/alpha * temp)    
  }
}