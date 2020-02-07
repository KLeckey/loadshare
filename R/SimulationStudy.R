
#' Simulation of exponentially distributed random variable with outliers
#'
#' Samples a vector of (by default) exponentially distributed random variables with rates
#' given in the vector \code{rates}. The sample can be contaminated with outliers
#' given as an alteration of single exponential samples by either the function
#' \code{LowerOutlierFunc} or by the function \code{UpperOutlierFunc}. A
#' contamination occurs with probability \code{LowerOutlierChance+UpperOutlierChance}
#' whereas the choice between \code{LowerOutlierFunc} and \code{UpperOutlierFunc}
#' is made with probabilities proportional to \code{LowerOutlierChance} and
#' \code{UpperOutlierChance}. The sample distribution can be altered by changing
#' \code{distr}. Note that the prediction methods in this paper are designed for
#' the exponential distribution only.
#'
#' @param rates a vector of rates. The i-th coordinate of the vector determines
#' the rate used for the rate of the i-th entry of the resulting sample vector
#' @param distr a sample distribution with two arguments: the number N of samples
#' and the N-dimensional vector \code{rates} of rates (cf. exponential distribution)
#' @param withOutliers a boolean variable. If set \code{TRUE}, the data is
#' contaminated with outliers created using the functions \code{LowerOutlierFunc}
#' and \code{UpperOutlierFunc}.
#' @param outlierExceptions an optional vector of positions in the resulting
#' vector which should never be altered by the outlier functions.
#' @param LowerOutlierChance,UpperOutlierChance the probabilities to apply the
#' functions \code{LowerOutlierFunc} and \code{UpperOutlierFunc} to create
#' outliers.
#' @param LowerOutlierFunc,UpperOutlierFunc functions used to create outliers.
#'
#' @return a vector of exponentially distributed random variables. If
#' \code{withOutliers} is \code{TRUE}, each random variable has a chance to be
#'  altered using the functions \code{LowerOutlierFunc} and \code{UpperOutlierFunc}
#'
#' @seealso \code{\link{singleSimulation}}
#'
#' @export

sample_w <- function(rates, distr=rexp, withOutliers = FALSE, outlierExceptions = NULL,
                     LowerOutlierChance = 0.1,
                     LowerOutlierFunc = function(x){return(1)},
                     UpperOutlierChance = 0.1,
                     UpperOutlierFunc = function(x){return(10*x)}){

  result <- distr(length(rates),rate=rates)

  if(withOutliers){
    N <- length(result)
    OutlierChance <- LowerOutlierChance + UpperOutlierChance

    # Sample for each position in result whether an outlier occurs (1) or not (0)
    OutlierPositions <- sample(c(0,1), N, replace=TRUE,
                               prob=c(1-OutlierChance, OutlierChance))
    OutlierPositions[outlierExceptions] <- 0

    # Sample which positions have upper outliers (1) and which have lower outliers (-1)
    OutlierType <- sample(c(-1,1), N, replace=TRUE,
                          prob=c(LowerOutlierChance, UpperOutlierChance))

    PositionLower <- which(OutlierPositions*OutlierType==-1)
    PositionUpper <- which(OutlierPositions*OutlierType==1)

    result[PositionLower] <- LowerOutlierFunc(result[PositionLower])
    result[PositionUpper] <- UpperOutlierFunc(result[PositionUpper])
  }
  return(result)
}

#' Generate sample distribution with white noise
#'
#' Generates a sample distribution with two arguments \code{N} and \code{rates}.
#' The distribution generates \code{N} random variables consisting of a
#' sum of exponentially distributed random variables with \code{rate=rates}
#' and normally distributed random variables with \code{mean} zero and a
#' standard deviation given by the parameter \code{sd}.
#'
#' @param sd the standard deviation used for the white noise
#' @param min.val the minimum value for the sampled random variabe, usually to
#' ensure nonnegativity. Default is set to 1.
#'
#' @return a function with arguments \code{N} and \code{rates} that generates
#' \code{N} random variables based on an exponential distribution with white noise
#'
#' @seealso \code{\link{sample_w}}
#'
#' @export

whiteNoiseSampler <- function(sd=1e4, min.val=1){

  return(function(N, rates){return(pmax(rexp(N, rate=rates)+rnorm(N, mean=0, sd=sd),min.val))})
}


#' Generate a lognormal sample distribution
#'
#' Generates a sample distribution with two arguments \code{N} and \code{rates}.
#' The distribution generates \code{N} random variables with a lognormal
#' distribution. The parameters \code{mean} and \code{sd} of the underlying
#' normal distribution are chosen in the following way: \code{mean} is
#' chosen such that the median of the lognormal distribution coincides with the
#' median of the exponential distribution with \code{rate=rates}.
#' and normally distributed random variables with \code{mean} zero and a
#' standard deviation given by the parameter \code{sd}. If \code{adjustSDtoExp=TRUE}
#' then \code{sd} is chosen in such a way that ratio betwenn the standard
#' deviation of the lognormal distribution and the one of the exponential
#' distribution is equal to the chosen value \code{SDExpRatio}. Otherwise,
#' the standard deviation is set to be the fixed value \code{sdFixed}
#'
#' @param adjustSDtoExp boolean variable to indicate whether the standard
#' deviation of the resulting distribution should be given by \code{SDExpRatio}-times
#' the standard deviation of the exponential distribution or should take the
#' fixed value \code{sdFixed}
#' @param SDExpRatio the ratio between the standard deviations if \code{adjustSDtoExp=TRUE}.
#' Default choice is 1.
#' @param sdFixed the value for the standard deviation of the lognormal distribution
#' if \code{adjustSDtoExp=FALSE}
#'
#' @return a function with arguments \code{N} and \code{rates} that generates
#' \code{N} random variables with lognormal distribution.
#'
#' @seealso \code{\link{sample_w}}
#'
#' @export

logNormSampler <- function(adjustSDtoExp=TRUE,
                           SDExpRatio=1, sdFixed=10){

  if(adjustSDtoExp){
    res <- function(N, rates){
      mean <- log(log(2))-log(rates)
      sd <- sqrt( 1/2+sqrt(1/4+SDExpRatio*rates/(log(2)^2)))
      return( exp(rnorm(N, mean=mean, sd=sd)))
    }
    return(res)
  }
  else{
    res <- function(N, rates){
      mean <- log(log(2))-log(rates)
      return( exp(rnorm(N, mean=mean, sd=sdFixed)))
    }
    return(res)
  }
}

#' Single simulation
#'
#' Perform a single simulation with specified parameters. The simulations
#' samples observations according to the parameters \code{obs_s,obs_Iv,obs_I}.
#' Based on the observation, predictions are produced using the
#' specified \code{methods}. Observed data can be contaminated with outliers
#' \code{withOutliers=TRUE}.
#'
#' @param theta a fixed model parameter used to determine the rates
#' of waiting times via \code{\link{getRates}} (in the Basquin Model)
#' @param obs_s,obs_Iv,obs_I parameters for the observed data. The observed
#' waiting times are exponentially distributed random variables with
#' rates obtained via \code{\link{getRates}}
#' @param pred_s,pred_start,pred_pos,pred_I parameter settings for the
#' experiment that is to be predicted. If \code{pred_start>1}, then
#' the first \code{pred_start-1} waiting times of the new experiment are
#' part of the observed data. The vector \code{pred_pos} contains the indices
#' of jump times that are to be predicted.
#' @param alpha,alpha1,alpha2 error choices for the predictions; see
#' \code{\link{predfail}} for more details
#' @param distr,withOutliers,LowerOutlierChance,LowerOutlierFunc,UpperOutlierChance,UpperOutlierFunc options
#' to add outliers to the data; see \code{\link{sample_w}} for more details
#' @param methods a vector of methods used for predictions. These methods
#' need to be valid choices in \code{\link{predfail}}
#'
#' @return a matrix that contains the average performance of each prediction method.
#' The i-th row of the matrix contains the average coverage rate, length, and
#' interval score of the method \code{methods[i]}
#'
#' @seealso \code{\link{sample_w}, \link{predfail}}
#'
#' @export

singleSimulation <- function(theta, obs_s, obs_Iv, obs_I,
                           pred_s, pred_start, pred_pos, pred_I,
                           alpha, alpha1, alpha2, distr, withOutliers,
                           LowerOutlierChance, LowerOutlierFunc,
                           UpperOutlierChance, UpperOutlierFunc,
                           methods){

  if(length(obs_I)==1){
    obs_I <- rep(obs_I, length(obs_s))
  }

  rates_obs <- getRates(theta, obs_s, obs_Iv, obs_I, h_basq)
  rates_pred <- getRates(theta, pred_s, max(pred_pos), pred_I, h_basq)
  rates <- c(rates_obs, rates_pred)

  if(withOutliers){
    outlierExceptions <- (length(rates_obs)+pred_start):length(rates)
  }
  else{
    outlierExceptions <- NULL
  }

  w <- sample_w(rates, distr, withOutliers, outlierExceptions,
                LowerOutlierChance, LowerOutlierFunc,
                UpperOutlierChance, UpperOutlierFunc)


  if(pred_start>1){
    obs_s <- c(obs_s, pred_s)
    obs_Iv <- c(obs_Iv, pred_start-1)
    obs_I <- c(obs_I, pred_I)

    pred_initial_value <- sum(w[(length(rates_obs)+1):
                                  (length(rates_obs)+pred_start-1)])
  }
  else{
    pred_initial_value <- 0
  }

  N_obs <- sum(obs_Iv)

  obs_data <- list(w=w[1:N_obs], s=obs_s, Iv=obs_Iv, I=obs_I)

  pred_times <- pred_initial_value +
                cumsum(w[(N_obs+1):(length(w))])

  pred_times <- pred_times[ pred_pos - pred_start+1]

  # Prediction intervals

  result <- matrix(0, nrow=length(methods), ncol=3,
                   dimnames = list(methods, c("coverage", "length", "score")))

  for(i in 1:length(methods)){
    prediction <- predfail(pred_I = pred_I, pred_s=pred_s,
                           pred_start = pred_start,
                           pred_pos=pred_pos,
                           pred_initial_value = pred_initial_value,
                           alpha=alpha, alpha1=alpha1, alpha2=alpha2,
                           obs_data = obs_data, method=methods[i])

    result[i,1] <- getCovRate(pred_times, prediction)
    result[i,2] <- getIntLength(prediction)
    result[i,3] <- getScore(pred_times, prediction, alpha = alpha)
  }
  return(result)
}

#' Repetition of \code{\link{singleSimulation}} with an increasing number of observations
#'
#' A function to apply \code{\link{singleSimulation}} several times
#' to data sets where \code{obs_s, obs_Iv, obs_I} are K-times
#' repetitions of the choice for \code{block_s, block_Iv, block_I}
#' for every K in \code{K_range}.
#'
#' @param K_range a vector of integers. Each element K of the vector
#' yields an application \code{\link{singleSimulation}} in which the length
#' of the observed data is proportional to K.
#' @param block_s,block_Iv,block_I choices which the parameters of
#' the observed data are based on.
#' @param theta,pred_s,pred_start,pred_pos,pred_I,distr,withOutliers,LowerOutlierChance,LowerOutlierFunc,UpperOutlierChance,UpperOutlierFunc,methods parameter
#' choices in \code{\link{singleSimulation}}
#' @param alpha the overall alpha-error in the prediction. Used in \code{\link{singleSimulation}}
#' @param alpha1Func a function on alpha and K. The return value is the
#' \code{alpha1}-value used in \code{\link{singleSimulation}}. The
#' \code{alpha2}-value in \code{\link{singleSimulation}} is set to be
#' \code{1-(1-alpha)/(1-alpha1)}
#'
#' @return a three-dimensional array \code{result} that contains the average performance of
#' each prediction method for each value of K in \code{K_range}.
#' The i-th row of the matrix contains the average coverage rate, length, and
#' interval score of the method \code{methods[i]}. The entry \code{result[ , ,K]}
#' is the return value of \code{\link{singleSimulation}} in the case
#' where the observed data is a \code{K_range[K]}-repetition of the parameters
#' \code{block_s, block_Iv, block_I}
#'
#'
#' @seealso \code{\link{singleSimulation}}
#'
#' @export


simulateKRepetition <- function(K_range, theta, block_s, block_Iv, block_I,
                                pred_s, pred_start, pred_pos, pred_I,
                                alpha, alpha1Func, distr, withOutliers,
                                LowerOutlierChance, LowerOutlierFunc,
                                UpperOutlierChance, UpperOutlierFunc,
                                methods){

  result <- array(0, dim=c(length(methods), 3, length(K_range)),
                   dimnames = list(methods, c("coverage", "length", "score"), K_range))

  if(length(block_I)==1){
    block_I <- rep(block_I, length(block_s))
  }

  for(i in 1:length(K_range)){

    obs_s <- rep(block_s, K_range[i])
    obs_Iv <- rep(block_Iv, K_range[i])
    obs_I <- rep(block_I, K_range[i])

    alpha1 <- alpha1Func(alpha, K_range[i])
    #alpha1 <- alpha - alpha2
    alpha2 <- 1-(1-alpha)/(1-alpha1) #choosen such that (1-alpha1)(1-alpha2)=(1-alpha)

    result[,,i] <- singleSimulation(theta=theta, obs_s=obs_s, obs_Iv=obs_Iv,
                                    obs_I=obs_I, pred_s=pred_s,
                                    pred_start=pred_start, pred_pos=pred_pos,
                                    pred_I=pred_I, alpha=alpha, alpha1=alpha1,
                                    alpha2=alpha2, distr=distr,
                                    withOutliers=withOutliers,
                                    LowerOutlierChance=LowerOutlierChance,
                                    LowerOutlierFunc=LowerOutlierFunc,
                                    UpperOutlierChance=UpperOutlierChance,
                                    UpperOutlierFunc=UpperOutlierFunc,
                                    methods=methods)

  }
  return(result)
}


#' Default choice for \code{alpha1Func}
#'
#' The default choice for the function \code{alpha1Func} used in simulations
#'
#' @param alpha,K input arguments (overall prediction error and number of repetitions)
#'
#' @return the choice for alpha1 in predictions
#'
#'
#' @seealso \code{\link{createSimulationFunction}, \link{createBatch}}
#'
#' @export

alpha1_default <- function(alpha,K){
  return(alpha/(2*(K)^(1/2)))
}


#' Create a wrapper function for simulations
#'
#' This function returns a function to start simulations. The simulation function
#' will have a single integer argument which either sets the seed
#' (if \code{setSeeds=TRUE}) or is a dummy argument with no effect.
#'
#' @param n_sims the number of simulations (iterations of \code{simulateKRepetition})
#' the returned function should perform.
#' @param setSeeds a boolean variable to choose whether the returned function
#' should use its input value to set a seed.
#' @param K_range,theta,block_s,block_Iv,block_I,pred_s,pred_start,pred_pos,pred_I,alpha,alpha1Func,distr,withOutliers,LowerOutlierChance,LowerOutlierFunc,UpperOutlierChance,UpperOutlierFunc,methods parameters
#' used in each of the \code{n_sims} iterations of \code{simulateKRepetition}
#'
#' @return a function with a single integer as an argument. The function
#' starts a simulation based on the chosen parameters.
#'
#' @seealso \code{\link{simulateKRepetition}, \link{createBatch}}
#'
#' @export

createSimulationFunction <- function(n_sims, setSeeds=TRUE, K_range=1:30,
                                    theta=c(28,3), block_s=c(200,100,80),
                                    block_Iv=c(6,6,6), block_I=35, pred_s=60,
                                    pred_start=4, pred_pos=6, pred_I=35,
                                    alpha=0.1, alpha1Func = alpha1_default,
                                    distr=rexp, withOutliers=FALSE,
                                    LowerOutlierChance=0.1,
                                    LowerOutlierFunc=function(x){return(1)},
                                    UpperOutlierChance=0.1,
                                    UpperOutlierFunc=function(x){return(10*x)},
                                    methods=c("3depth", "delta", "wald", "naive")){


  doSim <- function(placeholder){
    return(simulateKRepetition(K_range, theta, block_s, block_Iv, block_I,
                                pred_s, pred_start, pred_pos, pred_I,
                                alpha, alpha1Func, distr, withOutliers,
                                LowerOutlierChance, LowerOutlierFunc,
                                UpperOutlierChance, UpperOutlierFunc,
                                methods))
  }
  if(n_sims==1){
    simFunc <- function(seed){
      if(!missing(seed) & setSeeds){
        set.seed(seed)
      }
      return(doSim(1))
    }
  }
  else{
    simFunc <- function(seed){
      if(!missing(seed) & setSeeds){
        set.seed(seed)
      }
      sims <- lapply(1:n_sims, doSim)
      return(Reduce(x=sims, f="+")/n_sims)
    }
  }
  return(simFunc)
}

#' Create a registry for the \code{BatchJobs} package
#'
#' This function creates a registry to perform parallel simulations on a computer
#' cluster (LiDoNg at TU Dortmund). The function requires the \code{BatchJobs}
#' package
#'
#' @param id id of the registry; see \code{BatchJobs::makeRegistry}
#' @param seedlist a vector of seeds. The length of the vector also determines
#' the number of parallel jobs. The i-th jobs uses \code{seedlist[i]} as a seed
#' @param setSeeds boolean variable to choose whether seeds should be set or not
#' @param n_sims,K_range,theta,block_s,block_Iv,block_I,pred_s,pred_start,pred_pos,pred_I,alpha,alpha1Func,distr,withOutliers,LowerOutlierChance,LowerOutlierFunc,UpperOutlierChance,UpperOutlierFunc,methods parameters
#' used in \code{\link{createSimulationFunction}}. The i-th parallel job will
#' evaluate the simulation function at position \code{seedlist[i]}
#'
#' @return a registry based on \code{BatchJobs::makeRegistry}
#' containing the jobs for the parallel computation
#'
#' @seealso \code{\link{createSimulationFunction}, BatchJobs::makeRegistry}
#'
#' @export

createBatch <- function(id, seedlist=1:200, setSeeds=TRUE, n_sims=50,
                        K_range=1:20, theta=c(28,3), block_s=c(200,100,80),
                        block_Iv=c(6,6,6), block_I=10, pred_s=60, pred_start=4,
                        pred_pos=6, pred_I=10, alpha=0.1,
                        alpha1Func = alpha1_default,
                        distr=rexp, withOutliers=FALSE,
                        LowerOutlierChance=0.1,
                        LowerOutlierFunc=function(x){return(1)},
                        UpperOutlierChance=0.1,
                        UpperOutlierFunc=function(x){return(10*x)},
                        methods=c("3depth", "delta", "wald", "naive")){

  simFunc <- createSimulationFunction(n_sims, setSeeds, K_range, theta, block_s,
                                      block_Iv, block_I, pred_s, pred_start,
                                      pred_pos, pred_I, alpha, alpha1Func,
                                      distr, withOutliers,
                                      LowerOutlierChance, LowerOutlierFunc,
                                      UpperOutlierChance, UpperOutlierFunc,
                                      methods)

  reg <- makeRegistry(id = id)
  addRegistryPackages(reg, "loadshare")
  batchMap(reg, simFunc, seed=seedlist)

  return(reg)
}

#' Plot simulation results
#'
#' A function to plot the simulation  results
#'
#' @param simResult a three-dimensional array produced by
#' the simulation function returned by \code{\link{createSimulationFunction}}
#' or via the parallel computation induced by \code{\link{createBatch}}
#' @param type a vector with entries in {1,2,3}. Type 1 refers
#' to a plot of the coverage rates, type 1 yields a plot of the interval
#' lengths and type 3 provides a plot of the interval scores.
#' @param methods a vector of names of methods used to produce the simulation
#' results. These names will appear in the legend of the graphic.
#' @param ylabs,xlabs a vector of labels for y- and x-axis for
#' each of the plots produced by the function (the total number of plots
#' is equal to the length of \code{type})
#' @param ylims an optional parameter to specify the \code{ylim} in each plot.
#' If \code{length(type)>1}, then \code{ylims} needs to be a matrix with two
#' columns in which the i-th row represents the \code{ylim} choice for
#' \code{type[i]}
#' @param mSymbols a vector of \code{pch} choices for the points in each of the
#' \code{methods}. Has to have the same length as \code{methods}
#' @param lwd,cex additional plot options for the line width and the point size.
#' @param inMillions a vector of booleans with the same length as \code{type}.
#' The i-th entry specifies whether the y-axis in plot for \code{type[i]} should
#' be in millions. If \code{TRUE}, each y-value is divided by \code{1e6}
#'
#' used in \code{\link{createSimulationFunction}}. The i-th parallel job will
#' evaluate the simulation function at position \code{seedlist[i]}
#'
#' @return NULL
#'
#' @seealso \code{\link{createSimulationFunction}, \link{createBatch}}
#'
#' @export


plotSimResult <- function(simResult, type=c(1,3), methods=dimnames(simResult)[[1]],
                          ylabs=dimnames(simResult)[[2]][type],
                          xlabs=rep("Repetitions of stress ranges 200, 100, 80", length(type)),
                          ylims=NULL, alpha=0.1, mSymbols=1:length(methods),
                          lwd=2, cex=0.5, inMillions=c(FALSE,TRUE)){

  K_range=as.numeric(dimnames(simResult)[[3]])

  for(i in 1:length(type)){
    if(!is.null(ylims)){
      if(length(type)==1){
        ylim=ylims
      }
      else{
        ylim=ylims[i,]
      }
    }
    else{
      ylim=c(0.9*min(simResult[,type[i],]), 1.1*max(simResult[,type[i],]))/
                                                                (1e6)^(inMillions[i])
    }

    plot(NULL,NULL, xlim=c(0,max(K_range)), ylim=ylim, xlab=xlabs[i], ylab=ylabs[i])

    for(j in 1:length(methods)){
      points(K_range, simResult[j,type[i],]/(1e6)^(inMillions[i]), pch=mSymbols[j], cex=cex)
    }
    if(type[i]==1){
      abline(h=1-alpha, col="red", lwd=lwd)
      legend("bottomright", methods, pch=mSymbols)
    }
    else{
      legend("topright", methods, pch=mSymbols)
    }
  }
  return(invisible(NULL))
}



