#' Compute Negative Log-Likelihood
#'
#' Computes normal, conditional negative Loglikelihood for all x observations with treatment j
#'
#' @param jtrt treatment j
#' @param epsy residual vector of y \eqn{\epsilon = y - \mu_y}
#' @param mu_xst mean vector of x*
#' @param start_x Tmax x trt matrix, containing indicators for
#'                panel time - treatment combinations in data vector
#' @param nx number of observations
#' @param start_y Tmax x trt matrix, containing indicators for
#'                panel time - treatment combinations in outcome vector
#' @param sgmaj variance under treatment j
#' @param rhoj correlation coefficient under treatment j
#'
#' @return the log-Likelihood value for given point
#'
#' @family loglikelihoods

# computes the negativ log-likelihood for all  x-observations with treatment jtrt
#  mxst..n x 1 mean and variance of the latent utilities x_star
#  nu      ...   degrees of freedom for the t-distribution
#
# Author: Helga Wagner
# Last Change: Sep, 10,2012
##########################################################################


lcondx <- function(jtrt, epsy, mu_xst, start_x, nx, start_y, sgmaj, rhoj){

  Tmax <- length(sgmaj)
  Tmin <- which(start_x !=0)[1]

  if ((0.999- sum(rhoj^2)) < 0){
    loglikx <- -Inf
  }  else{
    loglikx <- 0
    for (t in Tmin:Tmax){
      sgmat <- sgmaj[1:t]
      rhot <- rhoj[1:t]
      n_xtj <- nx[t]
      ind <- start_x[t] + (0:(n_xtj-1))
      indy <- start_y[t] + (0:(t*n_xtj-1))
      epsy_tj <- matrix(as.vector(epsy[indy]), t,n_xtj)  #n_xtjx t vector
      # Likelihood contribution of x|y
      os2 <- (rhot/sgmat)
      sxc <- sqrt(1-sum(rhot^2))
      mxc <- t(mu_xst[ind]) +t(os2) %*% epsy_tj
      if (jtrt==0){
        likx <- pnorm(-(mxc/sxc))
      } else{
        if (jtrt==1){
          likx <- pnorm(mxc/sxc)
        }
      }
      loglikx <- loglikx + sum(log(likx))
    }
  }
  return(loglikx)
}
