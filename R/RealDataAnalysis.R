#' Load experiment data
#'
#' A function to load the experiment results performed in the civil engineering
#' department of the TU Dortmund university. The data contains
#' the stress ranges \code{stress}, the times \code{jumptimes} when wires break
#' and the times \code{t_jm} between consecutive wire breaks in each experiment.
#' Each experiment is based on a concrete beam with I=35 wires.
#' The data in experiment SB06a was obtained after increasing the stress
#' range in SB06 since only one wire break was observed after a long running
#' time. We recommend to remove the data from SB06a since it is not comparable
#' to the other data.
#'
#' @param mergeToList a boolean variable. If \code{TRUE}, the function returns
#' the data in form of a list that suits the format required in
#' \code{\link{DataSanityCheck}}
#' @param excludeSB06a a boolean variable to decide whether the results
#' from experiment SB06a should be removed from the list created when
#' \code{mergeToList=TRUE}
#'
#' @return The experiment data in form a list
#'
#' @export

loadData <- function(mergeToList=TRUE, excludeSB06a=TRUE){

  if(!mergeToList){
    return(list(jumptimes=loadshare::jumptimes,stress=loadshare::stress,t_jm=loadshare::t_jm))
  }
  else{
    data_w<-unlist(loadshare::t_jm, use.names = FALSE)

    data_Iv<-lengths(loadshare::t_jm, use.names=FALSE)
    data_J<-length(data_Iv)

    s_full <- unlist(loadshare::stress, use.names = FALSE)
    s <- s_full[ 1+cumsum(c(0,data_Iv[-data_J]))]

    if(excludeSB06a){
      data_w<-data_w[-((sum(data_Iv[-data_J])+1):length(data_w))]
      s<-s[1:11]
      data_Iv <- data_Iv[-data_J]
      data_J <- data_J-1
    }
    return(list(w=data_w, s=s, Iv=data_Iv, I=35))
  }
}

#' Split data into observed and prediction for a leave-one-out analysis
#'
#' A function to that split \code{data} into a list \code{observations} of
#' observed data and suitable information (\code{pred_initial_value,pred_jumps})
#' on data that is to be predicted. The observed data contains
#' all experiment except for (parts of) experiment number \code{num_experiment}.
#' The jump times up to index \code{obs_in_experiment} of experiment number
#' \code{num_experiment} are also added to the observed data in
#' \code{observations}. The remaining jump times are considered to be
#' part of the data that is to be predicted. In particular,
#' \code{pred_initial_value} equals either 0 (if \code{obs_in_experiment=0})
#' or is given as jump number \code{obs_in_experiment} of experiment number
#' \code{num_experiment}.
#'
#'
#' @param data the entire data set as a list supported by
#' \code{\link{DataSanityCheck}}
#' @param num_experiment the number of the experiment that is to be predicted
#' in the leave-one-out analysis
#' @param obs_in_experiment the number of jumps that are observed
#' in experiment number \code{num_experiment} which do not need to be predicted
#'
#' @return a list which contains another list \code{observations} of observed
#' data in the leave-one-out analysis and the information \code{pred_initial_value}
#' and \code{pred_jumps} of the data that needs to be predicted.
#'
#' @seealso \code{\link{predNext_LOO}}
#'
#' @export

leaveOneOutData <- function(data, num_experiment, obs_in_experiment){

  I <- data$I

  if(length(I)==1){
    I <- rep(I, length(data$Iv))
  }

  inds_experiment <- (cumsum(c(0,data$Iv))[num_experiment]+1):
                     (cumsum(data$Iv)[num_experiment])

  inds_prediction <- inds_experiment[(obs_in_experiment+1):
                                       length(inds_experiment)]

  obs_w <- data$w[-inds_prediction]

  if(obs_in_experiment>0){
    obs_s <- data$s
    obs_Iv <- data$Iv
    obs_Iv[num_experiment] <- obs_in_experiment
    obs_I <- I

  }
  else{
    obs_s <- data$s[-num_experiment]
    obs_Iv <- data$Iv[-num_experiment]
    obs_I <- I [-num_experiment]
  }

  pred_initial_value <- sum(c(0,data$w[inds_experiment])[1:(obs_in_experiment+1)])

  pred_jumps <- cumsum(data$w[inds_prediction])+pred_initial_value

  return(
    list(observations=list(w= obs_w, Iv=obs_Iv, s=obs_s, I=obs_I),
         pred_initial_value = pred_initial_value, pred_jumps = pred_jumps)
  )
}

#' Perform a leave-one-out analysis on data
#'
#' A function that performs a leave-one-out analysis (LOO) on \code{data}.
#' The for \code{future_steps} determines whether predictions
#' should be made for waiting times/immediate successors of jump times
#' (if \code{future_steps=1}) or for a sum of the next \code{future_steps}
#' waiting times after each jump time.
#'
#'
#' @param data the entire data set as a list supported by
#' \code{\link{DataSanityCheck}}
#' @param future_steps an integer to determine the index difference between
#' observations and the predicted jump time, i.e. the predicted random variable
#' will be a sum of \code{future_steps} exponentially distributed random variables
#' @param alpha,alpha1,alpha2 error choices in the predictions; see
#' \code{\link{predfail}} for more details.
#' @param methods a vector of methods supported by \code{\link{predfail}}.
#'
#' @return a list which contains a matrix \code{result} and an integer
#' \code{n_observations}. The matrix result contains the average score, coverage
#' rate and interval length in each of the methods, averaged over all predictions
#' made during the LOO. the integer \code{n_observations} counts the total number
#' of predictions made during the LOO.
#'
#' @seealso \code{\link{predfail},\link{leaveOneOutData}}
#'
#' @export

predNext_LOO <- function(data, future_steps = 1,
                         alpha=0.1, alpha1=1-sqrt(1-alpha), alpha2=1-sqrt(1-alpha),
                         methods=c("3depth", "delta", "wald", "LR", "naive")){

  J <- length(data$Iv)

  if(length(data$I)==1 & J>1){
    data$I <- rep(data$I, J)
  }

  l <- length(methods)

  scoreSum <- rep(0,l)     # ordered as in methodlist
  coverageSum <- rep(0,l)  # ordered as in methodlist
  lengthSum <- rep(0,l)    # ordered as in methodlist

  loo_count <- 0

  for(i in 1:J){
    Iv_curr <- data$Iv[i]

    if(Iv_curr < future_steps){
      next
    }
    for(j in 0:(Iv_curr - future_steps)){
      loo_count <- loo_count +1

      loo <- leaveOneOutData(data, i, j)

      pred_jump <- loo$pred_jumps[future_steps]


      for(k in 1:length(methods)){
        pred <- predfail(data$I[i], data$s[i], pred_start=j+1,
                         pred_pos=j+future_steps,
                         pred_initial_value=loo$pred_initial_value,
                         alpha1=alpha1, alpha2=alpha2, alpha = alpha,
                         obs_data=loo$observations, method=methods[k]
                         )
        scoreSum[k] <- scoreSum[k]+getScore(pred_jump, pred, alpha = alpha)
        coverageSum[k] <- coverageSum[k]+getCovRate(pred_jump, pred)
        lengthSum[k] <- lengthSum[k]+getIntLength(pred)
      }
    }
  }

  dimnameX <- c("Score (millions)", "Coverage rate", "Length (millions)")

  result <- matrix(c(scoreSum /(1e6*loo_count), coverageSum/loo_count,
                      lengthSum/(1e6*loo_count)),
                   nrow=3, byrow=TRUE,
                   dimnames=list(dimnameX, methods))

  return(list(result=result, n_observations=loo_count))
}
