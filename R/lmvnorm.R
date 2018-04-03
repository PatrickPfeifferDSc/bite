#' Compute log-Likelihood of multivariate Normals
#'
#' @param y outcome vector
#' @param mu vector of means
#' @param Sigmainv inverse of variance matrix
#'
#' @return Log-Likelihood value under multivariate Normal assumption.
#'
#' @family loglikelihoods


lmvnorm <- function(y, mu, Sigmainv){
  # for y Nx1
  n <- length(y)
  ycen <- y -mu
  return( -0.5 * (log(2*pi)*n - log(det(Sigmainv)) + ycen%*%Sigmainv%*%ycen))
}
