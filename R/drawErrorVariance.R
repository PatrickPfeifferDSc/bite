#' Draw Error Variance
#'
#' Internal.
#'
#' @param eps2 vector of squared residuals
#' @param indt indicators where panel times start.
#' @param indy0 index y vector
#' @param sn previous matrix to determine correct dimension of variance matrix
#' @param prior object containing prior information, prior standard deviance is used
#'

drawErrorVariance <- function(eps2,indt,indy0,sn, prior, fix.sigma = FALSE){
  if (fix.sigma){
    sgma2 <-  matrix(0, Tmax, 2)
    sgma2[,1] <- 0.25
    sgma2[,2] <- 1
    var_yt <- matrix(0,length(indy0),1)
    indy1 <- !indy0
    var_yt[indy0] <- indt[indy0,]%*%sgma2[,1]
    var_yt[indy1] <- indt[indy1,]%*%sgma2[,2]
    return(list(sgma2 = sgma2, var_yt = var_yt))
  }

  Tmax <- dim(sn)[1]

  sum_eps2 <- matrix(0, Tmax, 2)
  var_yt <- matrix(0,length(indy0),1)

  indy1 <- !indy0

  if (length(indy0) == (2*dim(indt)[1])){ # for the model with imputation
  indt <- rbind(indt, indt)
  }

  for (t in 1:Tmax){
  sum_eps2[t,1] <- sum(eps2[indt[,t] & indy0])
  sum_eps2[t,2] <- sum(eps2[indt[,t] & indy1])
  }

  Sn <- prior$S0 + sum_eps2/2
  sgma2 <- matrix(1/rgamma(dim(sn)[1]*dim(sn)[2], shape = as.vector(sn), scale = as.vector(1/Sn)), dim(sn)[1], dim(sn)[2])

  var_yt[indy0] <- indt[indy0,]%*%sgma2[,1]
  var_yt[indy1] <- indt[indy1,]%*%sgma2[,2]

  return(list(sgma2 = sgma2, var_yt = var_yt))
}
