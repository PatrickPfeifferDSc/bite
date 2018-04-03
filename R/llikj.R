#' Compute log-Likelihood
#'
#' Computes the negativ log-likelihood for all observations with treatment j
#'
#' @param sgmaj,rhoj vector of \eqn{2*Tmax} parameter of the variance/covariance matrix
#' @param j treatment
#' @param epsy residuals on y
#' @param mu_xst mean vector of x*
#' @param start_x Tmax x trt matrix, containing indicators for
#'                panel time - treatment combinations in data vector
#' @param nx number of observations
#' @param start_y Tmax x trt matrix, containing indicators for
#'                panel time - treatment combinations in outcome vector
#'
#' @return The log-likelihood value corresponding to input values.
#'
#' @family loglikelihoods

llikj <- function(sgmaj, rhoj,j,epsy,mu_xst,start_x,nx,start_y){
  Tmax <- length(sgmaj)
  Tmin <- which(start_x!=0)[1]
  if (0.999 - sum(rhoj^2) < 0){
    llik = Inf
  }else{
    logliky <- matrix(0, Tmax,1)
    for(t in Tmin:Tmax){
      sgmat <- sgmaj[1:t]
      n_xtj <- nx[t]
      indy <- start_y[t] + (0:(t*n_xtj - 1))
      epsy_tj <- matrix(as.vector(epsy[indy]), t, n_xtj)
      Sigmainv <- diag(1/sgmat^2)
      logliky[t] <- sum(loglik_mvnorm_balpan(epsy_tj, matrix(0, t, n_xtj), Sigmainv))
    }
    loglikx <- lcondx(j-1, epsy,mu_xst,start_x, nx, start_y, sgmaj, rhoj)
  }
  llik <- loglikx + sum(logliky)
  return(llik)
}
