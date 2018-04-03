#' Compute Maximum log-Likelihood
#'
#' [Internal Function] \code{comp_logml} will return a numeric value, the maximum log-likelihood
#' y = X*alpha+epsilon,  epsilon N(0,S), alpha ~ N(0,inv(iA0))
#'
#' ySy ...  y'S^{-1}y
#' XSy ...  X'S^{-1}y
#' invBN=X'S^{-1}X+ iB0
#' iB0     k x  k prior precision matrix
#' indic    k indicators for regression effects
#'
#' @return Returns an object containing the numeric value of the marginal Likelihood and beta parameters
#'
#' @param ySy matrix product of y vector normalized by standard deviation
#' @param XSy matrix product of data, normalizing matrix and y outcome
#' @param invBN inverse Beta matrix
#' @param indic indices of likelihood to be accessed
#' @param iB0 previous precision matrix
#'
#' @seealso \code{\link{lcondx}}

comp_logml <- function(ySy, XSy, invBN, indic, iB0){

  k <- length(indic)
  index <- which(indic!=0)
  invB0 <- iB0[index,index]  # get entries from the precision matrix

  if (invB0[1,1]==0){
    lDet_B0 <- - log(det(invB0[2:k,2:k]))
  } else{lDet_B0 <- - log(det(invB0))}

  BN <- solve(invBN[index,index])
  lDet_BN <- log(det(BN))
  bNh <- XSy[index]
  Q <- ySy - (t(bNh)%*% BN %*% bNh)
  h <- lDet_BN - lDet_B0
  lmarlik <- 0.5 * (h-Q)

  return(list(lmarlik=lmarlik, BN=BN, bNh=bNh))
}



