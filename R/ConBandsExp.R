#' Sample of standardized quantiles
#'
#' Computes sample quantiles of standardized exponentially distributed random
#' variables. The standardization works as follows: First compute the MLE for
#' \code{theta} based on the samples; see \code{\link{estML}} for the computation
#' of the MLE. Then compute the estimated rates by using the MLE for \code{theta}
#' in the model; see \code{\link{getRates}}. The standardized exponentially distributed
#' random variables then are defined as the product of the samples with their
#' estimated rates. Note that this yields asymptotically standard exponentially distributed
#' random variables if the estimation for \code{theta} is consistent.
#'
#' @param theta the true model parameter in the underlying model. Note that
#' the distribution of the resulting quantiles does not depend on \code{theta}
#' if the default \code{linkfun} is used. Hence the input is irrelevant in
#' the default setting.
#' @param s,Iv,I parameter settings in the experiment.
#' @param linkfun link function used for the rate. Default choice is
#' \code{\link{h_basq}}.
#'
#' @return a vector of sample quantiles of the standardized exponential distribution.
#'
#' @seealso \code{\link{getRates}, \link{h_basq}, \link{estML}, \link{getConBands}}
#'
#' @export

sampleStandQuants <- function(theta, s, Iv, I, linkfun=h_basq){
  rates <- getRates(theta, s, Iv, I, linkfun=linkfun)
  w <- sample_w(rates=rates)
  theta_hat <- estML(data=list(w=w, s=s, Iv=Iv, I=I), linkfun = linkfun)

  w_hat <- w*getRates(theta_hat, s, Iv, I, linkfun=linkfun)

  q_hat <- sort(w_hat)

  return(q_hat)
}

#' Confidence band for exponential distribution
#'
#' Compute a simultaneous \code{(1-alpha)}-confidence band for the sample quantiles
#' in the Basquin model. The confidence band compares the rescaled waiting times
#' with the standard exponential distribution. This rescaling is done as in
#' \code{\link{sampleStandQuants}}. The confidence band is computed based on the
#' Kolmogorov-Smirnov test statistic, i.e. by computing a band which contains
#' all sample quantiles for which the sup distance between
#' their empirical distribution function and the standard exponential distribution
#' is below the corresponding \code{alpha}-quantile.
#'
#' @param theta,s,Iv,I,linkfun parameter choices in \code{\link{sampleStandQuants}}
#' @param alpha the error in the confidence band.
#' @param n_sim the number of repetitions (simulations) on which the computation
#' of the confidence band is based. Increase for more precision, decrease to speed
#' up the computation. Default choice is one million.
#'
#' @return a list with three entries: \code{upper_band} and \code{lower_band}
#' are two vectores in which the ith entry contains the upper/lower bound for
#' the ith quantile of the exponential distribution. The third entry \code{alpha} is solely for
#' convenience and returns the \code{alpha} chosen as input.
#'
#' @seealso \code{\link{sampleStandQuants}}
#'
#' @export

getConBands <- function(alpha, n_sim=1e6, theta=c(28,3), s, Iv, I, linkfun=h_basq){

  N <- sum(Iv)

  res_matrix <- matrix(0, nrow=n_sim, ncol=N) #Matrix to store simulated waiting times
  rates <- getRates(theta=theta, s=s, Iv=Iv, I=I, linkfun=linkfun)

  for(i in 1:n_sim){
    res_matrix[i,] <- sampleStandQuants(theta=theta, s=s, Iv=Iv, I=I, linkfun=linkfun)
  }

  distKS <- rep(0,n_sim) #Initialize vector storing the Kolmogorov-Smirnov distances

  for(i in 1:n_sim){
    distKS[i] <- max(c(abs(1-exp(-res_matrix[i, ])-(1:N/N)),
                       abs(1-exp(-res_matrix[i, ])- (0:(N-1))/N)))
  }
  #Only keep samples below the (empirical) (1-alpha)-quantile:
  conbandFuns <- res_matrix[which(distKS <= quantile(distKS, probs=1-alpha)), ]

  up <- down <- rep(0, N) #Initialize upper/lower bound of the confidence band

  for(i in 1:N){
    up[i] <- max(conbandFuns[, i])
    down[i] <- min(conbandFuns[, i])
  }

  return(list(upper_band=up, lower_band=down, alpha=alpha))
}


#' Q-Q-plot with confidence bands
#'
#' Creates a Q-Q-plot of the theoretical quantiles \code{theo_q} and the sample
#' quantiles \code{sample_q}. Also adds a confidence band given as a list \code{bands}.
#' This confidence band can be computed via \code{\link{getConBands}}.
#'
#' @param bands a list providing the confidence band, see \code{\link{getConBands}}.
#' @param alpha the error chosen for \code{bands}. Only used for the default
#' \code{main} of the plot
#' @param theo_q,sample_q two vectors storing the theoretical- and sample quantiles
#' for the QQ-plot.
#' @param q_type the type for the QQ-plot. Default is a linear interpolation of
#' the quantiles, see \code{type} in \code{\link{plot}} for more options.
#' @param col_band,fill_bands,fill_col options for plotting the confidence band.
#' \code{col_band} and \code{fill_col} provide the colors for the border of the
#' confidence band and the fill color of the bands. The fill color is only used if
#' \code{fill_bands} is \code{TRUE}.
#' @param add_diagonal,lty_diag If \code{add_diagonal=TRUE}, the diagonal \code{y=x}
#' is added to the plot. The line type of the diagonal can be adjusted via \code{lty_diag}.
#' @param ylim,xlab,ylab,main parameter choices in \code{\link{plot}}.
#' The default setting \code{ylim=NULL} yields an automated choice for \code{ylim}.
#'
#' @return \code{invisible(NULL)}.
#'
#' @seealso \code{\link{sampleStandQuants},\link{plot}}
#'
#' @export

plotConBands <- function(bands, alpha=bands$alpha, theo_q, sample_q, q_type='l',
                         col_band="gray", fill_bands=TRUE, fill_col="gray",
                         add_diagonal=TRUE, lty_diag=2, ylim=NULL,
                         xlab="Quantiles of Exp(1)", ylab="Sample Quantiles",
                         main=paste("Q-Q-Plot with ", 100*(1-alpha),"%-confidence band", sep="")){
  if(is.null(ylim)){
    plot(theo_q, sample_q, xlab=xlab, ylab=ylab, main=main, type=q_type)
  }
  else{
    plot(theo_q, sample_q, xlab=xlab, ylab=ylab, main=main, type=q_type, ylim=ylim)
  }
  if(fill_bands){
    polygon(x=c(0,theo_q, rev(theo_q)), y=c(0,bands$upper_band, rev(bands$lower_band)),
            col=fill_col, border=col_band)
    points(theo_q, sample_q, type=q_type)
  }
  else{
    lines(theo_q, bands$upper_band, col=col_band)
    lines(theo_q, bands$lower_band, col=col_band)
  }
  if(add_diagonal){
    abline(a=0, b=1, lty=lty_diag)
  }
  return(invisible(NULL))
}


#' Test on exponential assumption
#'
#' A \code{alpha}-level significance test to determine whether the waiting times
#' of the given \code{data} follow exponential distributions with rates
#' determined via the link function \code{linkfun}, see \code{\link{getRates}}.
#' The test also produces a Q-Q-plot with a simultaneous \code{(1-alpha)}-confidence
#' band if \code{makePlot} is \code{TRUE}. Since prediction intervals based on
#' \code{\link{predfail}} assume that the waiting times follow exponential distribution,
#' we recommend using this test before applying \code{\link{predfail}}.
#' Note that the confidence bands are computed via simulation hence causing
#' this function to have a slow runtime. Reducing \code{n_sim} speeds up the
#' function but increases the approximation error due to the simulation.
#' Note that due to discrete time measurements and the continuity of the
#' exponential distribution, the test would always reject if any of the waiting
#' times equal zero. We thus recommend to increase zero values slightly, e.g.
#' by setting them to one instead (this is only recommended if the typical
#' value of a waiting time is way above one). This will be done if \code{increaseZeros}
#' is \code{TRUE}. In that case, all zeros are replaced by \code{zeroReplacement}.
#'
#' @param data a list storing four vectors: the waiting times \code{w},
#' the stresses \code{s} in each experiment, the number of failures \code{Iv}
#' in each experiment, the total number \code{I} of components in each experiment.
#' @param alpha the error chosen for the confidence bands computed via \code{\link{getConBands}}
#' @param linkfun the link function for the rates in the model. Default is the
#' Basquin link function \code{\link{h_basq}}.
#' @param increaseZeros a boolean variable. If \code{TRUE}, then all zeros among
#' the waiting times in \code{data} are replaced by the value \code{zeroReplacement}
#' @param zeroReplacement a positive number as a replacement for all zeros in
#' among the waiting times in \code{data}
#' @param n_sim the number of repetitions in the simulation of the confidence
#' band, see \code{\link{getConBands}}.
#' @param makePlot a boolean variable. If true, the function produces a QQ-plot
#' with a confidence band via \code{\link{plotConBands}}.
#' \code{col_band} and \code{fill_col} provide the colors for the border of the
#' confidence band and the fill color of the bands. The fill color is only used if
#' \code{fill_bands} is \code{TRUE}.
#' \@param q_type,cold_band,fill_bands,fill_col,xlab,ylab,main parameter choices
#' in \code{\link{plotConBands}}.
#'
#' @return a list with three entries: the first entry is a boolean variable
#' \code{rejectNull} which is \code{TRUE} if the test rejects the Null
#' Hypothesis, i.e. if the exponential assumption is rejected. The
#' second entry is a vector \code{pos_outside} which provides the positions
#' in the sample quantile vector which lie outside the confidence band. The
#' last entry \code{distance} yields the (directed) distances of the sample quantiles (with
#' index in \code{pos_outside}) to the confidence band. Negative distances
#' correspond to quantiles below the lower bound of the confidence band.
#'
#' @seealso \code{\link{getConBands},\link{plotConBands},\link{predfail}}
#'
#' @export


testExpDist <- function(data, alpha, linkfun=h_basq,
                        increaseZeros=TRUE, zeroReplacement=1,
                        n_sim=1e6, makePlot=TRUE, q_type='l', col_band="gray",
                        fill_bands=TRUE, fill_col="gray",
                        xlab="Quantiles of Exp(1)", ylab="Sample Quantiles",
                        main=paste("QQ-Plot with ", 100*(1-alpha),"% confidence band", sep="")){

  if(increaseZeros){
    data$w[which(data$w==0)] <- zeroReplacement
  }

  MLE <- estML(data, linkfun = linkfun)
  bands <- getConBands(alpha=alpha, n_sim=n_sim, theta=MLE, s=data$s, Iv=data$Iv,
                       I=data$I, linkfun=linkfun)
  # Compute standardized sample quantiles
  rates_hat <- getRates(theta=MLE, s=data$s, Iv=data$Iv, I=data$I, linkfun=linkfun)
  sample_q <- sort(data$w * rates_hat)
  N <- length(sample_q)

  if(makePlot){
    plotConBands(bands, alpha=alpha, theo_q= qexp( (1:N)/(N+1)), sample_q=sample_q,
                 q_type=q_type, col_band=col_band, fill_bands=fill_bands,
                 fill_col=fill_col, xlab=xlab, ylab=ylab, main=main)
  }
  where_above <- which(sample_q>bands$upper)
  where_below <- which(sample_q < bands$lower)
  where <- c(where_above, where_below)

  return(list(rejectNull= (length(where)>0), pos_outside=where,
              distances = c(sample_q[where_above]- bands$upper[where_above],
                            sample_q[where_below]- bands$lower[where_below])))
}

