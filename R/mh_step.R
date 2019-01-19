#' Computes a Metropolis Hastings Step
#'
#' [Internal Function] Performs Metropolis Hastings steps for
#' internal MCMC procedure. lpost_old/lpost_new are log-posteriors at
#' theta_old/theta_new
#'
#' @param theta_old old proposal parameter value
#' @param theta_new new proposal parameter value
#' @param lpost_old old log posterior
#' @param lpost_new new log posterior
#' @param lq_new log q(theta_new|theta_old) the proposal log-densities
#' @param lq_old log q(theta_old|theta_new) the proposal log-densities
#'
#' @return A list of the new theta accepted value, or the old
#' value if rejected and an acceptance indicator.

mh_step <- function(theta_old, theta_new, lpost_old, lpost_new, lq_new, lq_old){
  lold <- lpost_old - lq_old
  lnew <- lpost_new - lq_new
  u <- log(runif(1))
  acc <- (u <= min(0, lnew - lold))
  if  (acc){
    theta <- theta_new
  }else{
      theta <- theta_old
  }
  return(list(theta=theta, acc=acc))
}
