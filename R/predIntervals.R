
#' Estimation for the information matrix
#'
#' A function to compute the estimated information matrix for a dataset. Note
#' that the resulting matrix is in a format suited for further calculations.
#' In order to transform the return value to a properly scaled estimator for
#' the information matrix, a rescaling by the factor \code{1/N} is required,
#' where \code{N} denotes the number of observations, i.e.\code{N=sum(Iv)}.
#'
#' @param thetahat The MLE for theta. Can be computed via \code{\link{estML}}
#' @param s,Iv,I parameter settings in the experiment
#' @param linkfun,gradlinkfun Linkfunction and its gradient (with respect to
#' theta) for the model. Only needs to be specified if \code{BasquinModel=FALSE}
#' @param BasquinModel a boolean variable. If \code{TRUE} the Basquin model is
#' used. Otherwise, the model is specified via the choice of \code{linkfun}
#'
#' @return The estimated information matrix (multiplied by the factor
#' \code{N=sum(Iv)})
#'
#' @seealso \code{\link{predImat_basq}}
#'
#' @export

predImat<-function(thetahat, s, Iv, I, linkfun, gradlinkfun, BasquinModel){
  if(BasquinModel){
    return(predImat_basq(s, Iv, I))
  }
  if(length(I)==1){
    I<-rep(I,length(s))
  }
  #Vector of arguments for linkfun and gradlinkfun
  x_vec<-rep(s*I,Iv)/
    (rep(I,Iv)-getbreakvec(Iv))
  #Vector of elements that appear in the sum computing I(theta)
  #The first half of the vector contains the vector of first entries of gradlinkfun(...)
  #divided by linkfun(...), the second have contains the vector of second entries of
  #gradlinkfun(...) divided by linkfun
  lambda_vec<-gradlinkfun(thetahat,x_vec)/linkfun(thetahat,x_vec)
  lambda1_vec<-lambda_vec[1:(length(lambda_vec)/2)]
  lambda2_vec<-lambda_vec[-(1:(length(lambda_vec)/2))]

  #Compute the entries a11, a12, a21, a22 of the information matrix I(theta)
  a11<-sum(lambda1_vec^2)
  a12<-a21<-sum(lambda1_vec*lambda2_vec)
  a22<-sum(lambda2_vec^2)
  return(matrix(c(a11,a21,a12,a22),nrow=2))
}


#' 3-depth test statistics
#'
#' A function to compute the normalized 3-depth of the residual vector
#' \code{residuals}. The normalized 3-depth is given by
#' \code{N(d3(residuals)-0.25)}, where \code{N} denotes the length of
#' \code{residuals}. The equivalent representation used in this function is
#' based on Equation (3.2) in Kustosz, Leucht, Mueller (2014).
#'
#' @param residuals a vector of residuals. Note that the vector can be replaced
#' by \code{sign(residuals)} since the 3-depth is only based on sign changes.
#'
#' @return the normalized 3-depth of \code{residuals}
#'
#' @export

normed3depth <- function(residuals){
  N <- length(residuals)
  S <- cumsum( sign(residuals))
  diffs <- (c( S, S[N]- S))^2

  result <- ( 3/4*diffs[N] - 3/2* sum(diffs)/N)/(N-2)*(1+1/(N-1)) +
    3/4*(1+1/(N-1))*(1+2/(N-2))
  return(result)
}

#' Quantiles for 3-depth
#'
#' The (asymptotic) quantiles of the normalized 3-depth if the signs of the
#' residuals are i.i.d. unifomly in {-1,1} distributed random variables. By
#' default, asymptotic quantiles are used for residuals of length larger than 100.
#' Otherwise, simulated quantiles from an unpublished work by M. Horn
#' and C.H.Müller (2019) are used.
#'
#' @param p the p-value for which the quantile is computed
#' @param N the length of residual vector considered in the normalized 3-depth.
#' This value is only relevant if \code{useAsymptotics=FALSE}
#' @param useAsymptotics a boolean variable. If \code{TRUE}, asymptotic
#' quantiles from the package \code{rexpar} are used (not published on CRAN,
#' only available via GitHub). Otherwise, simulated quantiles for the specified
#' N are used.
#'
#' @return the p-quantile of the normalized 3-depth
#'
#' @seealso \code{\link{normed3depth}}
#'
#' @export

q3depth <- function( p, N=101, useAsymptotic=(N>100) ){
  if(useAsymptotic){
    temp <- floor(p*1000)
    result <- (p*1000 - temp)*loadshare::SimQuants[temp+2,2] +
              (temp+1-p*1000)*loadshare::SimQuants[temp+1,2]
    return(result)
  }
  else{
    temp <- floor(p*10000)
    unscaledResult <- (p*10000 - temp)*loadshare::quants3[temp+2,N-9] +
      (temp+1-p*10000)*loadshare::quants3[temp+1,N-9]
    return(N*(unscaledResult-0.25))
  }
}

#' Prediction intervals based on confidence sets
#'
#' A function that transforms a confidence set \code{ConSet} with error
#' \code{alpha2} to (a matrix of) (1-\code{alpha})(1-\code{alpha2})-prediction
#' intervals for waiting times with indices \code{pred_pos} in an
#' experiment with parameters \code{pred_I} and \code{pred_s}.
#' If \code{pred_star>1}, the first \code{pred_start}-1 waiting times of
#' the new experiments are assumed to be part of the observed input. In this
#' case, the last jumptimes (i.e. the sum of waiting times 1,...,
#' \code{pred_start}-1) needs to be provided as another variable
#' \code{pred_initial_value}.
#'
#' @param pred_I,pred_s model parameters s and I for the experiment that is to
#' be predicted
#' @param pred_start the first waiting time that is not observed in the
#' predicted experiment.
#' @param pred_pos a vector of indices indicating which waitingtimes are to be
#' predicted
#' @param pred_initial_value a number representing the last jumptime observed
#' in the new experiments before the prediction starts
#' @param alpha the error for each prediction interval
#' @param ConSet the confidence set in the data format provided by the function
#' \code{\link{ConSet}}
#' @param linkfun the linkfunction used for the model
#' @param approx an additional parameter to reduce computation time by
#' adding another approximation error. instead of considering all columns in
#' \code{ConSet}, only a number of \code{approx} equidistant ones is used.
#' @param auto_matrix_rep a boolean variable used in \code{\link{predHypoexp}}
#'
#' @return the union(s) over all (1-\code{alpha})-prediction intervals taken over
#' all parameters in \code{ConSet}. The result is given as a matrix with
#' 2 columns in which the i-th row stores the prediction for the jump time
#' at \code{pred_pos[i]}
#'
#' @seealso \code{\link{ConSet}, \link{predfail}}
#'
#' @export

predUnionConSet <- function(pred_I, pred_s, pred_start=1, pred_pos,
                          pred_initial_value=0, alpha=0.05, ConSet,
                          linkfun=h_basq, approx=dim(ConSet)[2],
                          auto_matrix_rep=FALSE){

  N<-dim(ConSet)[2]
  #Index increment in the for loop - only indices of the form 1+k*increment will be considered
  #when computing the union of the prediction intervals
  increment<-floor(N/approx)
  #Initial upper and lower interval bounds for each prediction interval
  #(updated in the first iteration of the for loop below)
  upper<-rep(-Inf, length(pred_pos))
  lower<-rep(Inf, length(pred_pos))

  #For empty ConSet this function returns NaNs for upper/lower bounds
  if(is.null(ConSet)){
    return(matrix(rep(NaN,2*length(upper)),ncol=2))
  }

  for(i in 1:approx){
    #The two extremal values for the rate based on the lower and upper theta1 bound in ConSet
    rate1 <- getRates(ConSet[-2,i*increment], pred_s, max(pred_pos), pred_I, linkfun)
    rate2 <- getRates(ConSet[-1,i*increment], pred_s, max(pred_pos), pred_I, linkfun)

    #Compute new prediction intervals with these rate and update upper and lower if necessary
    for(j in 1:(length(pred_pos))){
      new_Int1<-predHypoexp(alpha,rate1[pred_start:(pred_pos[j])],
                            auto_matrix_rep = auto_matrix_rep)
      new_Int2<-predHypoexp(alpha,rate2[pred_start:(pred_pos[j])],
                            auto_matrix_rep = auto_matrix_rep)
      upper[j]<-max(upper[j],new_Int1[2],new_Int2[2])
      lower[j]<-min(lower[j],new_Int1[1],new_Int2[1])
    }
  }
  return(matrix(c(lower,upper),ncol=2,
                dimnames=list(pred_pos, c("lower", "upper") ))+pred_initial_value)
}



#' Prediction intervals based on the delta-method
#'
#' A function to compute \code{(1-alpha1)(1-alpha2)}-prediction intervals based
#' on the delta-method.
#'
#' @param pred_I,pred_s model parameters s and I for the experiment that is to
#' be predicted
#' @param pred_start the first waiting time that is not observed in the
#' predicted experiment.
#' @param pred_pos a vector of indices indicating which waitingtimes are to be
#' predicted
#' @param pred_initial_value a number representing the last jumptime observed
#' in the new experiments before the prediction starts
#' @param alpha1 error in the correction to the MLE of the asymptotic quantiles
#' used as prediction intervals
#' @param alpha2 choice for the quantiles estimated via MLE. Up to corrections
#' made according to the choice of \code{alpha1}, the prediction interval
#' coincides with a naive \code{(1-alpha2)}-prediction interval obtained via
#' \code{\link{predHypoexp}} with \code{alpha=alpha2}.
#' @param obsData observed data which the prediction is based on. The data
#' has to be formated as in \code{\link{DataSanityCheck}}
#' @param BasquinModel a boolean variable to indicate whether the Basquin model
#' is used. Note that the current version of this package only support the
#' Basquin model.
#' @param linkfun,gradlinkfun linkfunction and its gradient if a model
#' other than the Basquin model is used. Note that this has no effect in
#' the current version since only the Basquin model is supported.
#' @param auto_matrix_rep a boolean variable used in \code{\link{predHypoexp}}
#'
#' @return a matrix of \code{(1-alpha1)(1-alpha2)}-prediction intervals
#' in which the i-th row stores the prediction for the jump time
#' at \code{pred_pos[i]}
#'
#' @seealso \code{\link{predfail}}
#'
#' @export

pred_delta <- function(pred_I, pred_s, pred_start, pred_pos,
                       pred_initial_value, alpha1, alpha2, obs_data,
                       BasquinModel=TRUE, linkfun=NULL, gradlinkfun=NULL,
                       auto_matrix_rep = FALSE){

  if(BasquinModel){
    linkfun <- h_basq
    gradlinkfun <- grad_basq
  }
  else{
    stop("The delta method is only implemented for the Basquin model")
  }

  thetahat <- estML(data=obs_data, linkfun=linkfun, BasquinModel=BasquinModel)
  Imat <- predImat(thetahat=thetahat,s=obs_data$s,
                   Iv=obs_data$Iv, I=obs_data$I,
                   linkfun=linkfun, gradlinkfun=gradlinkfun,
                   BasquinModel=BasquinModel)
  invI <- inv_cramer(Imat)

  N <- length(pred_pos)
  predInts <- matrix(0,ncol=2,nrow=N, dimnames=list(pred_pos, c("lower","upper")))

  rates <- getRates(thetahat, pred_s, max(pred_pos), pred_I, linkfun)

  for(j in 1:N){

    naive <- predHypoexp(alpha2,rates[pred_start:(pred_pos[j])],
                         auto_matrix_rep = auto_matrix_rep)

    grad_b1 <- matrix(grad_qhypoexp_basq(thetahat, alpha2/2 , pred_start-1,
                                         pred_pos[j], pred_I, pred_s), ncol=1)

    grad_b2 <- matrix(grad_qhypoexp_basq(thetahat, 1-alpha2/2 , pred_start-1,
                                         pred_pos[j], pred_I, pred_s),ncol=1)

    v1 <- qnorm(1-alpha1/2)*sqrt(drop(t(grad_b1)%*%invI%*%grad_b1))
    v2 <- qnorm(1-alpha1/2)*sqrt(drop(t(grad_b2)%*%invI%*%grad_b2))

    predInts[j,] <-  naive +c(-v1,v2)+pred_initial_value
  }
  return(predInts)
}

#' Function to produce prediction intervals based on several methods
#'
#' A function to compute either \code{(1-alpha)}-prediction intervals
#' based on the naive approach or \code{(1-alpha1)(1-alpha2)}-prediction
#' intervals based on either one of the confidence set based approaches
#' (Wald, LR, 3-depth) or the delta-method.
#'
#' @param pred_I,pred_s model parameters s and I for the experiment that is to
#' be predicted
#' @param pred_start the first waiting time that is not observed in the
#' predicted experiment.
#' @param pred_pos a vector of indices indicating which waitingtimes are to be
#' predicted
#' @param pred_initial_value a number representing the last jumptime observed
#' in the new experiments before the prediction starts
#' @param alpha1 error in the confidence set (if "wald", "LR", or "3depth" are
#' chose as a method) or the correction to the MLE of the asymptotic quantiles
#' used as prediction intervals (if "delta" is chosen as a method);
#' see \code{\link{ConSet},\link{pred_delta}} for more details.
#' @param alpha2 choice for \code{alpha} in \code{\link{predHypoexp}}
#' if any method other than "naive" is chosen; see
#' \code{\link{predUnionConSet}, \link{pred_delta}} for more details.
#' @param auto_matrix_rep_hypoexp choice for \code{aut_matrix_rep} in
#' \code{\link{predHypoexp}}
#' @param obsData observed data which the prediction is based on.
#' The data should be a list that stores a vector of waiting times
#' in \code{obsData$w} and the experiment settings of the observed
#' data in \code{obsData$s}, \code{obsData$Iv}, \code{obsData$I};
#' see also \code{\link{DataSanityCheck}} for the data format check
#' @param BasquinModel a boolean variable to indicate whether the Basquin model
#' is used.
#' @param linkfun,gradlinkfun linkfunction and its gradient if a model
#' other than the Basquin model is used.
#' @param method a method choice. Supported choices are "naive", "wald", "LR",
#' "3depth" and "delta".
#' @param grid_xlim,grid_ylim,mesh_width,checkFullGrid,autoExpandGrid,expandX,expandY,maxExpansion_xlim,maxExpansion_ylim parameter
#' choices in \code{\link{ConSet}}
#'
#' @return a matrix of prediction intervals in which the i-th row stores the
#' prediction for the jump time at \code{pred_pos[i]}
#'
#' @seealso \code{\link{predHypoexp}, \link{pred_delta}, \link{predUnionConSet}, \link{ConSet}}
#'
#' @export

predfail <- function(pred_I, pred_s, pred_start=1, pred_pos,
                     pred_initial_value=0, alpha1=0.05, alpha2=0.05,
                     alpha = alpha1+alpha2, auto_matrix_rep_hypoexp=FALSE,
                     obs_data=list(w=NULL, s=NULL, Iv=NULL, I=NULL),
                     BasquinModel=TRUE,  linkfun=NULL, gradlinkfun=NULL,
                     method="3depth", grid_xlim=c(20,35), grid_ylim=c(0,5),
                     mesh_width=0.05, checkFullGrid=FALSE,
                     autoExpandGrid = TRUE, expandX=5, expandY=5,
                     maxExpansion_xlim = c(0,100),
                     maxExpansion_ylim = c(0,50)
                     ){

  DataSanityCheck(data=obs_data)

  if(BasquinModel){
    linkfun <- h_basq
    gradlinkfun <- grad_basq
  }

  ML <- estML(obs_data, linkfun, BasquinModel=BasquinModel)

  if(method=="naive"){
    pred_naive <- matrix(0,nrow=length(pred_pos),ncol=2,
                         dimnames=list(pred_pos, c("lower", "upper") ))

    rates <- getRates(ML, pred_s, max(pred_pos), pred_I, linkfun)

    for(i in 1:(length(pred_pos))){

      pred_naive[i,] <- predHypoexp(alpha,  rates[pred_start:(pred_pos[i])],
                                    auto_matrix_rep = auto_matrix_rep_hypoexp)+
                        pred_initial_value
    }
    return(pred_naive)
  }
  if(method=="3depth" | method=="wald" | method=="LR"){
    CS<- ConSet(data=obs_data, method=method,
                alpha=alpha1, BasquinModel=BasquinModel,
                linkfun=linkfun, gradlinkfun=gradlinkfun,
                grid_xlim=grid_xlim, grid_ylim=grid_ylim,
                mesh_width=mesh_width, checkFullGrid=checkFullGrid,
                autoExpandGrid = autoExpandGrid,
                expandX=expandX, expandY=expandY,
                maxExpansion_xlim = maxExpansion_xlim,
                maxExpansion_ylim = maxExpansion_ylim)

    pred <- predUnionConSet(pred_I, pred_s, pred_start,
                            pred_pos, pred_initial_value, alpha2,
                            CS, linkfun, auto_matrix_rep = auto_matrix_rep_hypoexp)
    return(pred)
  }

  if(method=="delta"){
    if(!BasquinModel){
      stop("delta method only implemented for the Basquin model")
    }

    pred <- pred_delta(pred_I, pred_s, pred_start, pred_pos,
                       pred_initial_value, alpha1, alpha2,
                       obs_data, BasquinModel, linkfun, gradlinkfun,
                       auto_matrix_rep=auto_matrix_rep_hypoexp)
    return(pred)
  }
  else{
    stop("Unknown method. Choose between: naive, 3depth, wald, LR, delta")
  }
}




