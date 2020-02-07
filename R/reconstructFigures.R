# Each figure in the article and the supplementary material can be reproduced
# with the corresponding functions below. In cases where the computation
# takes a long time (e.g. creating data in the simulation study), the data is
# provided in .RData files.

#' Compute predictions in Figure 2
#'
#' @export

get_log_wij_preds <- function(data=loadData(), alpha=0.1, alpha1=1-sqrt(1-alpha),
                              alpha2=alpha1, precision = 0.01){
  I <- data$I
  if(length(I)==1){
    I <- rep(I, length(data$s))
  }

  # Confidence set with 3-depth:
  CS_depth <- ConSet(data, method="3depth",alpha=alpha2, mesh_width = precision)

  # predicitions for each stress level:
  stress <- rep(data$s*I, data$Iv)/(rep(I, data$Iv) - getbreakvec(data$Iv))

  # remove unnecessary duplicants:
  stress <- sort(unique(stress))

  # based on the model we can predict each wire break independently from another
  # using the entries in stress for s and an arbitray number for I
  # since the first wire break is unaffected by the choice for I

  preds_delta <- preds_depth <- matrix(0, nrow=3, ncol=length(stress))

  for(i in 1:length(stress)){
    preds_delta[2:3, i] <- log(
      predfail(pred_I=35, pred_s=stress[i], pred_start=1, pred_pos=1,
               alpha1=alpha1, alpha2=alpha2, alpha = alpha,
               obs_data=data, method="delta"), base=10)
    preds_depth[2:3, i] <- log(
      predUnionConSet(pred_I = 35, pred_s = stress[i], pred_start=1, pred_pos=1,
                      alpha=alpha1, ConSet=CS_depth), base=10)

  }
  preds_delta[1,] <- preds_depth[1, ] <- stress

  return(list(delta=preds_delta, depth=preds_depth))
}


# Labels for Figure 1
Fig1_xlab <- "sI/(I-i)"
Fig1_ylab <- expression(paste("log"[10],"(w"["ij"],")",sep=""))
Fig1_legend1_txt <- c("200MPa", "455MPa", "200MPa", "150MPa", "98MPa", "200MPa",
                 "100MPa", "60MPa", "80MPa", "80MPa", "50MPa")
Fig1_legend2_txt=c("Fitted expectation of the exponential distribution",
                  expression(paste("Pointwise 90%-prediction intervals using ",delta,"-method", sep="")),
                  "Pointwise 90%-prediction intervals using 3-depth")

#' Reproduce Figure 1
#'
#' @export

getFigure1 <- function(data=loadData(), xlim=c(0,750),ylim=c(0,12), xlab=Fig1_xlab,
                       ylab=Fig1_ylab, legend1_txt=Fig1_legend1_txt,
                       legend2_txt=Fig1_legend2_txt, cex=1, lwd=1,
                       cex.lab=cex, pt.cex=cex, y.intersp=1){

  thetahat <- estML(data)
  preds <- get_log_wij_preds()

  ## Initialize plot ##
  plot(NULL,NULL, xlim=xlim, ylim=ylim, xlab=xlab,
       ylab=ylab, cex.lab=cex.lab)

  collist <- rep(rainbow(length(data$s)), data$Iv)

  stress <- rep(data$s*data$I, data$Iv)/(data$I-getbreakvec(data$Iv))

  ## Plot data points ##
  for(i in 1:length(stress)){
    points(stress[i],log(max(data$w[i],1),base=10),col=collist[i], pch=16, cex=pt.cex)
  }

  legend("topright", legend1_txt, col=rainbow(11), pch=rep(16,11),
         pt.cex=rep(pt.cex,11), cex=cex, y.intersp=y.intersp)


  ## Plot predictions ##
  lines(x=preds$delta[1,], y=preds$delta[2,],lty=2, lwd=lwd)
  lines(x=preds$delta[1,], y=preds$delta[3,], lty=2, lwd=lwd)

  lines(x=preds$depth[1,], y=preds$depth[2,],lty=3, lwd=lwd)
  lines(x=preds$depth[1,], y=preds$depth[3,], lty=3, lwd=lwd)

  ## Plot fitted expectations based on MLE ##
  x <- preds$delta[1,]
  y <- -log(h_basq(thetahat, x),base=10)
  lines(x=x, y=y, lty=1, lwd=lwd)


  legend("topleft",legend2_txt, lwd=lwd,
           lty=c(1,2,3), cex=cex, bty="n", y.intersp=y.intersp)

  return(invisible(NULL))
}

#' Reproduce Figure 2
#'
#' @export

getFigure2 <- function(n_sim=1e6, xlab="Quantiles of Exp(1)",
                       ylab="Sample Quantiles", ylim=c(0,15),
                       main="Q-Q-Plot with 90%-confidence band"){


  real_data <- loadshare::loadData()
  theta <- loadshare::estML(real_data)
  real_hat <- real_data$w*loadshare::getRates(theta, real_data$s, real_data$Iv,
                                              real_data$I)

  bands <- getConBands(alpha=0.1, n_sim=n_sim, theta=c(28,3), s=real_data$s,
                       Iv=real_data$Iv, I=real_data$I)

  plotConBands(bands=bands, alpha=bands$alpha, theo_q=qexp((1:137)/138),
               sample_q=sort(real_hat), ylim=ylim, xlab=xlab, ylab=ylab, main=main)


  return(invisible(NULL))
}



#' Reproduce Table 1
#'
#' Note that the leave-one-out cross validation is time consuming. Hence the
#' computation of the table can take around 30 minutes to finish.
#'
#' @export

getTable1 <- function(){
  data <- loadData()

  return(list(Next=predNext_LOO(data),
              Next5=predNext_LOO(data, future_steps = 5)))
}

# Labels for Figure 3a
Fig3a_legend_txt <- c(expression(paste(delta,"-method", sep="")), "3-depth",
                      "SB06")

#' Reproduce Figure 3a
#'
#' @export

getFigure3a <- function(xlim = c(0,3), ylim=c(0,150), legend_txt=Fig3a_legend_txt,
                        main="", ylab="Load cycles in millions",
                        xlab="Number of broken tension wires", lwd=1){


  obs_data <- loadData(excludeSB06a=FALSE)

  end_obs <- sum(obs_data$Iv[1:10])
  jumptimes_SB06 <- obs_data$w[end_obs+1]
  obs_data$Iv <- obs_data$Iv[1:10]
  obs_data$s <- obs_data$s[1:10]
  obs_data$w <- obs_data$w[1:end_obs]
  if(xlim[2]-xlim[1]<5){
    plot(NULL, NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
         main=main, xaxt = "n")
    axis(1, at = (xlim[1]):(xlim[2]))
  }
  else{
    plot(NULL, NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
         main=main)
  }

  #Plot SB06 data as a step function
  addstepfun(c(0,1,2,2), c(0, jumptimes_SB06/1e6, 108), col="red", lwd=lwd)

  #Add an arrowhead to the end of the step function
  arrows(x0=2, x1=2, y0=jumptimes_SB06/1e6, y1=108, col="red", lwd=lwd, code=2,
         length=0.1)

  # Compute and plot prediction intervals based on delta method #
  pred_del <- predfail(35, 50, pred_start = 1, pred_pos = 1:(xlim[2]-1),
                       obs_data=obs_data, method="delta")

  addstepfun(c(1,1:(xlim[2])), c(0, pred_del[,1]/1e6), lwd=lwd, lty=2)
  addstepfun(c(1,1:(xlim[2])), c(0, pred_del[,2]/1e6), lwd=lwd, lty=2)

  # Compute and plot prediction intervals based on 3-depth #
  pred_depth <- predfail(35, 50, pred_start = 1, pred_pos = 1:(xlim[2]-1),
                         obs_data=obs_data, method="3depth")

  addstepfun(c(1,1:(xlim[2])), c(0, pred_depth[,1]/1e6), lwd=lwd, lty=3)
  addstepfun(c(1,1:(xlim[2])), c(0, pred_depth[,2]/1e6), lwd=lwd, lty=3)


  legend("topleft", legend=legend_txt,
         lty=c(2,3,1), col=c("black", "black", "red"), lwd=rep(lwd,3), bty="n")

  return(invisible(NULL))
}

# Labels for Figure 3b
Fig3b_legend_txt <- c(expression(paste(delta,"-method", sep="")), "3-depth",
                      "SB06a")

#' Reproduce Figure 3b
#'
#' @export

getFigure3b <- function(xlim = c(0,17), ylim=c(0,30), legend_txt=Fig3b_legend_txt,
                        main="", ylab="Load cycles in millions",
                        xlab="Number of broken tension wires", lwd=1){


  obs_data <- loadData(excludeSB06a=FALSE)

  end_obs <- sum(obs_data$Iv[1:11])
  jumptimes_SB06a <- cumsum(obs_data$w[(end_obs+1):(sum(obs_data$Iv[1:12]))])
  obs_data$Iv <- obs_data$Iv[1:11]
  obs_data$s <- obs_data$s[1:11]
  obs_data$w <- obs_data$w[1:end_obs]

  plot(NULL, NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)

  #Plot SB06a data as a step function
  addstepfun(0:18, c(0, jumptimes_SB06a/1e6), col="red", lwd=lwd)

  # Compute and plot prediction intervals based on delta method (note that only 34 wires remain) #
  pred_del <- predfail(34, 120, pred_start = 1, pred_pos = 1:17,
                       obs_data=obs_data, method="delta")

  addstepfun(c(1,1:18), c(0, pred_del[,1]/1e6), lwd=lwd, lty=2)
  addstepfun(c(1,1:18), c(0, pred_del[,2]/1e6), lwd=lwd, lty=2)

  # Compute and plot prediction intervals based on 3-depth #
  pred_depth <- predfail(34, 120, pred_start = 1, pred_pos = 1:17,
                         obs_data=obs_data, method="3depth")

  addstepfun(c(1,1:18), c(0, pred_depth[,1]/1e6), lwd=lwd, lty=3)
  addstepfun(c(1,1:18), c(0, pred_depth[,2]/1e6), lwd=lwd, lty=3)


  legend("topleft", legend=legend_txt,
         lty=c(2,3,1), col=c("black", "black", "red"), lwd=rep(lwd,3), bty="n")

  return(invisible(NULL))
}


#' Reproduce Figure 4a
#'
#' @export

getFigure4a <- function(cex=0.5){
  plotSimResult(loadshare::resNo, type=1, methods=c("3-depth", "Wald",
                                         expression(paste(delta,"-method", sep="")),
                                         "Naive"), cex=cex,
                ylabs="Coverage rate", ylims=c(0.7,1), inMillions=FALSE)
}

#' Reproduce Figure 4b
#'
#' @export

getFigure4b <- function(cex=0.5){
  plotSimResult(loadshare::resNo, type=3, methods=c("3-depth", "Wald",
                                         expression(paste(delta,"-method", sep="")),
                                         "Naive"), cex=cex,
                ylabs="Interval score (millions)", ylims=c(10,40), inMillions=TRUE)
}

#' Reproduce Figure 5a
#'
#' @export

getFigure5a <- function(cex=0.5){
  plotSimResult(loadshare::resSym, type=1, methods=c("3-depth", "Wald",
                                         expression(paste(delta,"-method", sep="")),
                                         "Naive"), cex=cex,
                ylabs="Coverage rate", ylims=c(0.5,1), inMillions=FALSE)
}

#' Reproduce Figure 5b
#'
#' @export

getFigure5b <- function(cex=0.5){
  plotSimResult(loadshare::resSym, type=3, methods=c("3-depth", "Wald",
                                          expression(paste(delta,"-method", sep="")),
                                          "Naive"), cex=cex,
                ylabs="Interval score (millions)", ylims=c(10,50), inMillions=TRUE)
}



################################################################################
############## Figures in the supplementary material ###########################
################################################################################

#' Reproduce Figure 1a
#'
#' @export

getSupFigure1a <- function(cex=0.5){
  plotSimResult(loadshare::resLow, type=1, methods=c("3-depth", "Wald",
                                          expression(paste(delta,"-method", sep="")),
                                          "Naive"), cex=cex,
                ylabs="Coverage rate", ylims=c(0.6,1), inMillions=FALSE)
}

#' Reproduce Figure 1b
#'
#' @export

getSupFigure1b <- function(cex=0.5){
  plotSimResult(loadshare::resLow, type=3, methods=c("3-depth", "Wald",
                                          expression(paste(delta,"-method", sep="")),
                                          "Naive"), cex=cex,
                ylabs="Interval score (millions)", ylims=c(10,40), inMillions=TRUE)
}




#' Reproduce Figure 2a
#'
#' @export

getSupFigure2a <- function(cex=0.5){
  plotSimResult(loadshare::resUp, type=1, methods=c("3-depth", "Wald",
                                         expression(paste(delta,"-method", sep="")),
                                          "Naive"), cex=cex,
                ylabs="Coverage rate", ylims=c(0,1), inMillions=FALSE)
}

#' Reproduce Figure 2b
#'
#' @export

getSupFigure2b <- function(cex=0.5){
  plotSimResult(loadshare::resUp, type=3, methods=c("3-depth", "Wald",
                                         expression(paste(delta,"-method", sep="")),
                                          "Naive"), cex=cex,
                ylabs="Interval score (millions)", ylims=c(10,90), inMillions=TRUE)
}


################################################################################
######### Additional figures not included in the paper #########################
################################################################################


#' Additional Predictions for the last failure of SB04
#'
#' @export

getPred19SB04 <- function(data=loadData(), alpha=0.1, alpha1 = 1-sqrt(1-alpha), alpha2 = alpha1){

  preds_delta <- preds_3depth <-  matrix(0, nrow=data$Iv[9], ncol=2)

  for(i in 0:(data$Iv[9]-1)){
    data_loo <- leaveOneOutData(data, 9,i)


    preds_3depth[i+1,] <- predfail(pred_I=35, pred_s=data$s[9], pred_start=i+1,
                                   pred_pos=data$Iv[9],
                                   pred_initial_value=data_loo$pred_initial_value,
                                   alpha1=alpha1, alpha2=alpha2, alpha = alpha,
                                   obs_data=data_loo$observations, method="3depth")

    preds_delta[i+1, ] <- predfail(pred_I=35, pred_s=data$s[9], pred_start=i+1,
                                   pred_pos=data$Iv[9],
                                   pred_initial_value=data_loo$pred_initial_value,
                                   alpha1=alpha1, alpha2=alpha2, alpha = alpha,
                                   obs_data=data_loo$observations, method="delta")
  }
  return( list(depth = preds_3depth, delta = preds_delta))
}

# Labels for Figure A1
Fig3_xlab <- "Number of observed wire breaks"
Fig3_ylab <- expression(paste("Log. prediction of the 19"^"th"," wire break", sep=""))
Fig3_legend_txt <- c(expression(paste("time of the 19"^"th"," break", sep="")),
                     expression(paste(delta,"-method", sep="")),
                     "3-depth")


#' Reproduce Figure A1
#'
#' @export

getFigureA1 <- function(xlab=Fig3_xlab, ylab=Fig3_ylab, legend_txt=Fig3_legend_txt,
                       ylim = c(7,8.1), lwd = 2, dist_ints = 0.2,
                       data=loadData(), preds=getPred19SB04(data)){

  pred_jumptime <- sum(data$w[(sum(data$Iv[1:8])+1):sum(data$Iv[1:9])])

  N <- dim(preds$delta)[1]

  plot(NULL, NULL, xlim=c(-1,N), ylim= ylim, xlab=xlab, ylab=ylab)

  # Add prediction intervals to the plot
  for(i in 1:N){
    pos <- i-1-dist_ints/2

    temp <- log(preds$delta[i,], base=10)
    lines(x=c(pos+dist_ints, pos+dist_ints), y=temp, lwd=lwd, lty=2)

    temp <- log(preds$depth[i,], base=10)
    lines(x=c(pos,pos), y=temp, lwd=lwd, lty=3)
  }

  # Add time of 19th wire break in SB04 to the plot
  log_jump <- log(pred_jumptime, base=10)
  abline(h=log_jump, lwd=lwd)

  # Add legend to the plot
  legend("topright", legend_txt, col="black", lty=c(1,2,3), lwd=rep(lwd,3))

  return(invisible(NULL))
}
