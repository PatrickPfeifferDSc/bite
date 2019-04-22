#' Compute log-Likelihood of Balanced Panels
#'
#' Internal. \code{loglik_mvnorm_balpan} computes for balanced panel data
#' \eqn{dim(y) = t x N} with equal \eqn{\Sigma} across observations
#'
#' @param y data
#' @param mu mean matrix
#' @param Sigmainv inverse variance matrix
#'
#' @return A scalar value, the log-Likelihood for given data, mu and Sigma.
#'
#' @family log-Likelihoods

loglik_mvnorm_balpan <- function(y, mu =  matrix(0, dim(y)[1] , dim(y)[2]), Sigmainv){  # idea of default matrix 0
  t <- dim(y)[1]
  ycen <- y - mu
  llik <- (-0.5)* (log(2*pi)*t - log(det(Sigmainv)) + rowSums((t(ycen) %*% Sigmainv) * t(ycen)))
  return(llik)
}
