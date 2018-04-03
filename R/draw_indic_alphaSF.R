#' Draw indicators and regression coefficients
#'
#' [Internal Function] \code{draw_indic_alphaSF} draws indicator values and
#' regression coefficients for the selection model (probit model of
#' choosing treatment)
#'
#' The procedure also includes variable selection process, utilizing a
#' Dirac spike prior distribution for variable selection
#' single move update
#' y... n x 1 response from fixed effects model
#' X... n x k regressor matrix
#'
#' @param y outcome vector
#' @param X data matrix
#' @param delta delta index of selection
#' @param deltafix fixed delta indices, not subject to variable selection
#' @param omega correlation parameter
#' @param invA0 inverse alpha entries matrix
#' @param isel indicator whether selection is performed
#'
#' @return delta ...  k x 1 indicators
#' @return alpha ... k x 1  regression effects
#'
#' @import stats

draw_indic_alphaSF <- function(y, X, delta, deltafix, omega, invA0, isel){

  df <- length(delta)
  nd <- df - sum(deltafix)

  invAN <- (t(X)%*%X) + invA0
  Xy <- t(X)%*%y
  yy <- t(y)%*%y

  logml.obj <- comp_logml(yy,Xy,invAN,delta,invA0)
  lmlik_old <- logml.obj$lmarlik; AN <- logml.obj$BN; aNh<- logml.obj$bNh;

  if (isel==1){
    idvar <- which(deltafix==0)
    delta_ranord <- sample(1:nd, nd, replace=FALSE)

    for (i in 1:nd){

      j <- idvar[delta_ranord[i]] # Selection of j-th not fixed effect
      dj <- delta[j]
      lpri <- dj*log(omega)+(1-dj)*(log(1-omega))
      lold <- lmlik_old + lpri
      # compute marginal likelihood for  model, where delta(j) is changed
      delta_new <- delta
      delta_new[j] <- 1-dj
      lpri_new <- (1-dj)*log(omega)+dj*(log(1-omega))
      logml.obj <- comp_logml(yy,Xy,invAN,delta_new,invA0)
      lmlik_new <- logml.obj$lmarlik; AN_new <- logml.obj$BN; aNh_new <- logml.obj$bNh;
      lnew <- lmlik_new + lpri_new
      # Draw new indicator delta(j)
      # compute posterior probabilities
      lh <- c(lnew,lold)
      maxl <- max(lh)
      prh <- exp(lh - maxl)
      postprob <- prh/sum(prh)
      # Update indicator
      if (runif(1) <= postprob[1]){
        delta[j] <- delta_new[j]
        lmlik_old <- lmlik_new
        aNh <- aNh_new
        AN <- AN_new
      }
    }
  }

  # Draw the regression effects
  aN <- AN %*% aNh
  ha <- cbind(rnorm(length(aN)))
  alpha_new <- t(chol(AN))%*%ha + aN

  alpha <- matrix(0, df, 1)
  idx <- as.logical(delta)
  alpha[idx] <- alpha_new
  return(list(delta=delta, alpha = alpha))
}
