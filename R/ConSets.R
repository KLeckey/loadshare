#' Computation of confidence sets
#'
#' Lists the boundary points of confidence sets computed for two-dimensional
#' parameters. Possible methods for the computation are the robust 3-depth method
#' and the classical methods based on Wald's test and the Likelihood ratio test.
#' The computation is based on a grid search in a specified region. The region
#' is expanded automatically if necessary. This auto expansion can be turned off
#' by setting \code{autoExpandGrid} to \code{FALSE}.
#'
#' @param data a list containing the observed data. The list has to be given in
#' format \code{list(w=waitingtimes, s=stresses, Iv=maxFailures, I=totalWires)}
#' @param method method to obtain a confidence set. Valid choices are "3depth", "wald" and "LR"
#' @param alpha the setting for the confidence level. The output has an (asymptotic) \code{(1-alpha)} confidence level.
#' @param BasquinModel a boolean variable indicating whether the default Basquin model
#' should be used. If set \code{TRUE} then \code{linkfun} and \code{gradlinkfun} will be changed
#' to the corresponding functions in the Basquin model.
#' @param linkfun the linkfunction used for the model
#' @param gradlinkfun the derivative of \code{linkfun} with respect to the parameter \eqn{\theta}
#' @param grid_xlim the horizontal grid limits for the first coordinate of the parameter \eqn{\theta}
#' @param grid_ylim the vertical grid limits for the second coordinate of the paramter \eqn{\theta}
#' @param mesh_width the mesh widh of the grid
#' @param checkFullGrid a boolean variable. If \code{TRUE} the algorithm will scan the entire grid
#' to compute the confidence set. If \code{FALSE} the algorithm will stop once
#' a full connected component of the confidence set is found.
#' @param autoExpandGrid a boolean variable. If \code{TRUE} the grid will be expanded
#' if the algorithm could not find a full connected component of the confidence set.
#' The grid will be expanded in all direction if no point in the set was found before.
#' Otherwise, expansion occur in the directions where the boundary of the current part
#' of the set overlaps with the margin of the grid.
#' @param expandX,expandY the horizontal and vertical increments if the grid is
#' expanded according to \code{autoExpandGrid}
#' @param maxExpansion_xlim,maxExpansion_ylim the maximal horizontal  and vertical
#' limit that the grid may expand to.
#'
#' @return a matrix with three rows and a variable number of columns. Each column
#' \eqn{(x,y,z)^\top} represents the leftmost point (x,z) and the rightmost point (y,z)
#' of the confidence set on the horizontal level \eqn{\theta_2=z}.
#'
#' @examples
#' stresses <- c(200,100)
#' max_failures <- c(7,3)
#' totalWires <- 35
#' rates <- h_basq(theta=c(28,3), x=c(35*200/(35-0:6), 35*100/(35-0:2)))
#' waitingtimes <- rexp(10, rate=rates)
#' data <- list(w=waitingtimes, s=stresses, Iv=max_failures, I=totalWires)
#' ConSet(data, method="3depth")
#' ConSet(data, method="wald")
#'
#' @seealso \code{\link{GridSearch}, \link{ConSet2d_wald}}
#'
#' @export

ConSet<-function(data, method="3depth",
                 alpha=0.05, BasquinModel=TRUE, linkfun=NULL, gradlinkfun=NULL,
                 grid_xlim=c(20,35), grid_ylim=c(0,5),
                 mesh_width=0.05, checkFullGrid=FALSE,
                 autoExpandGrid = TRUE, expandX=5, expandY=5,
                 maxExpansion_xlim = c(0,100),
                 maxExpansion_ylim = c(0,50)
                ){

  if(missing(data)){
    stop("Data is missing with no default.")
  }
  DataSanityCheck(data)


  if(BasquinModel){
    linkfun <- h_basq
    gradlinkfun <- grad_basq
  }
  if(method=="wald"){
    return(ConSet2d_wald(data=data, alpha=alpha, BasquinModel=BasquinModel,
                         linkfun=linkfun, gradlinkfun=gradlinkfun,
                         precision=mesh_width,
                         thetahat=estML(data,linkfun, BasquinModel=BasquinModel))
           )
  }

  N<-length(data$w)
  I <- data$I
  if(length(I)==1){
    I <- rep(I, length(data$Iv))
  }

  #Couple data in matrix: first col=waitingtimes, second col= argument x for linkfun
  Z<-matrix(0,nrow=N,ncol=2)
  Z[,1] <- data$w
  Z[,2] <- rep(data$s*I,data$Iv)/
    (rep(I,data$Iv)-getbreakvec(data$Iv))

  #Sort rows of Z in such a way that the second column is increasing:
  Z<-Z[order(Z[,2]),]
  if(method!="3depth"){
    #ML estimator required
    thetahat <- estML(data=data,linkfun=linkfun, BasquinModel=BasquinModel)
  }
  if(method=="3depth"){
    #Get the asymptotic alpha quantile of the 3depth from Kustosz et al.
    qalpha <- q3depth(alpha)
  }

  #Boolean Function to decide whether theta is in the confidence set
  if(method=="3depth"){
    InConSet<-function(theta){
      depth <- normed3depth(Z[,1]-log(2)/linkfun(theta,Z[,2]))
      return(depth>=qalpha)
    }
  }
  if(method=="LR"){
    InConSet<-function(theta){
      lam1<-linkfun(theta,Z[,2])
      lam2<-linkfun(thetahat, Z[,2])
      temp<- -2*sum(log(lam1/lam2))+2*sum( Z[,1]*(lam1-lam2))
      return(temp<=qchisq(1-alpha,df=length(theta)))
    }
  }
  if(method!="3depth" & method!="LR"){
      stop("Unknown method. Choose type between 3depth, wald and LR.")
  }


  expanding <- TRUE
  while(expanding){
    expanding <- FALSE
    temp <- GridSearch(InSetFun=InConSet, grid_xlim=grid_xlim, grid_ylim=grid_ylim,
                       mesh_width=mesh_width, checkFullGrid=checkFullGrid)

    expReq <- temp$expansionRequirements
    if(autoExpandGrid){
      if(expReq$left & (grid_xlim[1]>maxExpansion_xlim[1])){
        grid_xlim[1] <- max(grid_xlim[1] - expandX, maxExpansion_xlim[1])
        expanding <- TRUE
      }
      if(expReq$right & (grid_xlim[2]<maxExpansion_xlim[2])){
        grid_xlim[2] <- min(grid_xlim[2] + expandX, maxExpansion_xlim[2])
        expanding <- TRUE
      }
      if(expReq$bottom & (grid_ylim[1]>maxExpansion_ylim[1])){
        grid_ylim[1] <- max(grid_ylim[1] - expandY, maxExpansion_ylim[1])
        expanding <- TRUE
      }
      if(expReq$top & (grid_ylim[2]<maxExpansion_ylim[2])){
        grid_ylim[2] <- min(grid_ylim[2] + expandY, maxExpansion_ylim[2])
        expanding <- TRUE
      }
    }
  }
  return(temp$result)
}


#' Grid search for confidence sets
#'
#' A grid search that finds points of a set on a given grid. The set is represented
#' by a boolean function \code{InSetFun} which returns true if and only if the input
#' is part of the set.
#'
#' @param InSetFun a boolean function to represent the set
#' @param grid_xlim,grid_ylim,mesh_width,checkFullGrid options for the grid search.
#' Details can be found in \code{\link{ConSet}}
#'
#' @return a matrix with three rows and a variable number of columns.
#' Details can be found in \code{\link{ConSet}}
#'
#' @seealso \code{\link{ConSet}, \link{GridSearchX}}
#'
#' @export

GridSearch <- function(InSetFun, grid_xlim, grid_ylim,
                       mesh_width, checkFullGrid){
  x_candidates <- seq(grid_xlim[1], grid_xlim[2], by=mesh_width)
  y_candidates <- seq(grid_ylim[1], grid_ylim[2], by=mesh_width)

  foundStart <- FALSE
  result <- NULL

  for(y in y_candidates){
    temp <- GridSearchX(InSetFun=InSetFun, y=y, x_candidates=x_candidates)
    if(!is.null(temp)){
      result <- c(result, temp,y)
      if(!foundStart){
        foundStart <- TRUE
      }
    }
    else{
      if(foundStart & (!checkFullGrid)){
        break
      }
    }
  }

  if(!is.null(result)){
    result <- matrix(result, nrow=3)
    expansionRequirements <- list(left=any(result[1,]==x_candidates[1]),
                                  right=any(result[2,]==x_candidates[2]),
                                  top=(result[3,dim(result)[2]]==y_candidates[length(y_candidates)]),
                                  bottom=(result[3,1]==y_candidates[1])
                                  )
  }
  else{
    expansionRequirements <- list(left=TRUE, right=TRUE, top=TRUE, bottom=TRUE)
  }

  return(list(result=result, expansionRequirements=expansionRequirements))
}

#' Auxiliary function for grid search
#'
#' A function to find the minimal- and maximal x-coordinate of points in a set
#' for a fixed y coordinate. This function is used repeatedly to perform
#' \code{\link{GridSearch}}.
#'
#' @param InSetFun a boolean function. Details can be found in \code{\link{GridSearch}}.
#' @param y a real number for the fixed y-coordinate considered in the set.
#' @param x_candidates a vector of x-coordinates which are tested during
#' grid search on level \code{y}.
#'
#' @return a vector of length two containing the minimal- and maximal x values
#' of points inside the set that have y-coordinate \code{y}. The function returns
#' \code{NULL} if none of the pointes listed in \code{x_candidates} belong to the set.
#'
#' @seealso \code{\link{GridSearch}}

GridSearchX <- function(InSetFun, y, x_candidates){
  InFun <- function(x){
    return(InSetFun(c(x,y)))
  }
  temp <- sapply(x_candidates, InFun)
  if(!any(temp)){
    return(NULL)
  }
  else{
    pos <- which(temp)
    return(c(x_candidates[min(pos)], x_candidates[max(pos)]))
  }
}


#' Computation of confidence sets based on Wald's test
#'
#' An alternative implementation for confidence sets avoiding grid search.
#' This faster implementation makes use of the fact that the confidence set
#' based on Wald's test is an ellipse.
#'
#' @param data,alpha,BasquinModel,linkfun,gradlinkfun parameters chosen in \code{\link{ConSet}}
#' @param precision the precision in the representation of the ellipse. This value
#' corresponds to the \code{mesh_width} in \code{\link{ConSet}}
#' @param thetahat the maximum likelihood estimation (MLE) for the parameter
#' theta. The MLE can be found via \code{\link{estML}}.
#'
#' @return a matrix with three rows and a variable number of columns.
#' Details can be found in \code{\link{ConSet}}
#'
#' @seealso \code{\link{ConSet}}
#'
#' @export

ConSet2d_wald <- function(data, alpha, BasquinModel, linkfun, gradlinkfun,
                          precision, thetahat){

  if(BasquinModel){
    Imat <- predImat_basq(data$s, data$Iv, data$I)
    linkfun <- h_basq
    gradlinkfun <- grad_basq
  }
  else{
    Imat <- predImat(thetahat, data$s, data$Iv, data$I, linkfun, gradlinkfun,
                     BasquinModel=FALSE)
  }
  q <- qchisq(1-alpha,df= 2)
  theta2_diff <- sqrt(q*Imat[1,1]/(Imat[2,2]*Imat[1,1]-(Imat[1,2])^2))
  theta2_vec <- seq(-theta2_diff, theta2_diff, precision)


  # Check whether first/last element of theta2_vec needs to be removed due to
  # rounding errors
  if( q/Imat[1,1]+(theta2_vec[1])^2*
      (( (Imat[1,2])^2 - Imat[2,2]*(Imat[1,1]))/((Imat[1,1])^2)) < 0){
    theta2_vec <- theta2_vec[-c(1, length(theta2_vec))]
  }

  temp <- sqrt(q/Imat[1,1]+theta2_vec^2*(( (Imat[1,2])^2 - Imat[2,2]*(Imat[1,1]))/((Imat[1,1])^2)))


  theta1_lower_vec <- thetahat[1] - theta2_vec * Imat[1,2]/Imat[1,1] - temp
  theta1_upper_vec <- thetahat[1] - theta2_vec * Imat[1,2]/Imat[1,1] + temp

  theta2_vec <- theta2_vec + thetahat[2]

  return( matrix( c(theta1_lower_vec, theta1_upper_vec, theta2_vec), nrow=3, byrow = TRUE))
}


#' A plot function for confidence sets
#'
#' A simple plot function for confidence sets given in the return format of \code{\link{ConSet}}.
#'
#' @param set a matrix with three rows and a variable number of columns representing
#' the confidence set. See the return value of \code{\link{ConSet}} for more details.
#'
#' @param alpha the alpha value chosen in \code{\link{ConSet}}. The chosen value
#' will be shown in the default plot title.
#' @param add,col,xlim,ylim,main,xlab,ylab plot options. See \code{\link{plot}} and \code{\link{par}}.
#'
#' @return \code{NULL}
#'
#' @seealso \code{\link{ConSet}}
#'
#' @export

plotConSet<-function(set, add=FALSE, col="blue", alpha=0.05,
                     xlab=expression(theta[1]),ylab=expression(theta[2]),
                     xlim=c(0.99*min(set[1:2,]),1.01*max(set[1:2,])),
                     ylim=c(0.99*min(set[3,]),1.01*max(set[3,])),
                     main=paste(100*(1-alpha),'% Confidence set for \\Theta',sep="")){
  if(dim(set)[1]!=3){
    stop("set has to be given as a matrix with 3 rows;
         see the output of the function ConSet for more details.")
  }
  if(!add){
    plot(NULL,NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  }
  if(dim(set)[2]>1){
    prec<-set[3,2]-set[3,1]
    for(i in 1:length(set[1,])){
      rect(set[1,i],set[3,i]-prec/2,set[2,i],set[3,i]+prec/2,col=col,border=col)
    }
  }
  else{
    lines(c(set[1,1],set[2,1]), c(set[3,1],set[3,1]), col=col, lwd=5)
  }
  return(invisible(NULL))
}




