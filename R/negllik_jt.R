#' Compute Negative log-Likelihood
#'
#' [Internal Function]
#'
#' @param theta parameter of variance/covariance matrix of dimension 2Tmax x 1
#' @param j treatment
#' @param t panel time
#' @param epsy residuals of y
#' @param mu_xst mean vector of x*
#' @param start_x Tmax x trt matrix, containing indicators for
#'                panel time - treatment combinations in data vector
#' @param nx dimension of x
#' @param start_y Tmax x trt matrix, containing indicators for
#'                panel time - treatment combinations in outcome vector
#' @param sgmaj variance for treatment j
#' @param rhoj correlation for treatment j
#'
#' @return A scalar which is the negative log-Likelihood at the given point.
#'
#' @family loglikelihoods

negllik_jt <- function(theta, j, t, epsy, mu_xst, start_x, nx, start_y, sgmaj, rhoj){

  Tmax <- length(sgmaj)
  Tmin <- which(start_x!=0)[1]

  sgmaj[t]<- exp(theta)  # scalar, but may be updated simultaneouly?

  if ((0.999 - sum(rhoj^2)) < 0) {
    llikn <- Inf
  } else{
    logliky <- matrix(0, Tmax,1)
    for (t in Tmin:Tmax){
      sgmat <- sgmaj[1:t]
      n_xtj <- nx[t]
      indy <- start_y[t] + (0:(t*n_xtj - 1))
      epsy_tj <- matrix(as.vector(epsy[indy]), t, n_xtj)
      Sigmainv <- diag(1/sgmat^2)
      sin_llhnorm <- loglik_mvnorm_balpan(epsy_tj, matrix(0, t, n_xtj), Sigmainv) # single loglikelihood contributions of obs to theta (function has default mu = zero matrix)
      logliky[t] <- sum(sin_llhnorm)
    }
  loglikx <- lcondx(j-1, epsy, mu_xst, start_x, nx, start_y, sgmaj, rhoj)
  llikn <- -(sum(logliky)) - loglikx
  }
  return(llikn)
}
