#' Compute Normal log-Likelihood
#'
#' Internal.
#' \code{lnorm} computes log-Likelihood of a Normal distributed sample
#'
#' @param x sample
#' @param mu,varinv mean, variance parameters
#'
#' @return The negative log-Likelihood of given parameters.
#'
#' @family loglikelihoods

lnorm <- function(x, mu, varinv){
  n <- length(as.vector(x))
  # if ((k > 1) || (dim(mu)[1] != 1) || (dim(mu)[1] > 1) || ((dim(varinv)[2] > 1))){
  #   stop("error: wrong dimensions of input variables")
  # }
  if (is.vector(varinv)) return(-0.5 * ( log(2*pi) - log(varinv) + (t(x-mu) %*% varinv %*% (x-mu))))
  if (dim(varinv)[1] == n){
    return(-0.5 * (n * log(2*pi) - log(varinv) + t(x-mu)%*%varinv %*%(x-mu)))
  }else{
    if (dim(varinv)[1] == 1){
    return(-0.5 * (n * log(2*pi) - n*log(varinv) + t(x-mu) %*% varinv %*% (x-mu)))
    }
    stop("error: wrong dimensions of input variables")
  }
}
