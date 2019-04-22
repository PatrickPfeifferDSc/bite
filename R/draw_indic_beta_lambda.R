#' Draw Indices for beta and lambda coefficients
#'
#' Internal. \code{draw_indic_beta_lambda} draws beta coefficients for the outcome
#' models (y) and lambda parameters for shared factor modelling.
#'
#' y....Tn x1 ...response from fixed effects model
#' x .... nx 1 ... treatments
#' X ....Tn x df regressor matrix
#'
#' @param y outcome vector
#' @param X data matrix
#' @param sgma2 variance
#' @param delta vector of indices which model parameters are exposed to
#'              the spike and slab variable selection process
#' @param deltafix vector indicating which variables are fixed
#' @param omega1 omega parameter of beta
#' @param omega2 omega parameter of lambda
#' @param dy dimension of outcome
#' @param invB0 startvalue matrix of beta coefficients
#' @param isel indicator whether selection is performed
#'
#' @import Matrix

drawBetaLambdaIndicesSF <- function(y, X, sgma2, delta, deltafix, omega1, omega2, dy, invB0, isel, fix.beta=FALSE){

  if (fix.beta){
    delta <- c(1)

    return(list(delta=delta, beta=beta))
  }

  df <- length(delta)
  nd <- df - sum(deltafix)

  ysy <- t(y/sgma2) %*% y
  XS <- t(X) %*% Matrix::Diagonal(x = as.vector(1/sgma2))   # pck use 'Matrix' transforms XS into dgCMatrix class
  invBN <- XS %*% X + diag(invB0)
  Xy_full <- XS %*% y
  iB0 <- diag(invB0)
  logml.obj <- comp_logml(ysy, Xy_full, invBN, delta, iB0)
  lmlik_old <- as.numeric(logml.obj$lmarlik); BN <- as.matrix(logml.obj$BN); bNh <- logml.obj$bNh;
  # casting types on return values for further use, usually would be class dgeMatrix
  # print(c(ysy, invBN[1:3,1:3], head(Xy_full, 3), lmlik_old, BN[1:2,1:2], bNh[1:3]))

  if (isel==1){
    idvar <- which(!deltafix)                   # draw indices which are not fixed
    delta_ranord <- sample(1:nd,length(1:nd), replace=FALSE)   # draw the next indices in random order
    for (i in 1:nd){
      j <- idvar[delta_ranord[i]] # Selection of j-th not fixed effect
      # delta_new: change deltaj
      dj <- delta[j]
      delta_new <- delta
      delta_new[j] <- 1 - dj
      if (j <= dy){
        lpri <- dj*log(omega1) + (1-dj)*(log(1-omega1))
        lpri_new <- (1-dj)*log(omega1) + dj*(log(1-omega1))
      } else{
        lpri <- dj*log(omega2) + (1-dj)*(log(1-omega2))
        lpri_new <- (1-dj)*log(omega2)+dj*(log(1-omega2))
      }
      # compute log posterior
      lold <- lmlik_old + lpri
      # print(lold)
      logml.obj <- comp_logml(ysy, Xy_full, invBN, delta_new, iB0)
      lmlik_new <- as.numeric(logml.obj$lmarlik); BN_new <- as.matrix(logml.obj$BN); bNh_new <- logml.obj$bNh;
      # print(c("new maxlik:", lmlik_new, "BN_new", BN_new[1:2,1:2],"bnhnew", bNh_new[1:3]))
      lnew <- lmlik_new + lpri_new
      # compute posterior probabilities and update delta(j)
      lh <- c(lnew,lold)
      maxl <- max(lh)
      prh <- exp(lh - maxl)
      post.prob <- prh/sum(prh)
      # print(c("post.prob" , post.prob))
      # Update indicator
      if (is.finite(post.prob[1])){
        if (runif(1) <= post.prob[1]){
          delta[j] <- delta_new[j]
          lmlik_old <- lmlik_new
          bNh <- bNh_new
          BN <- BN_new
        }
      }
    }
  }
  # Draw the regression effects
  bN <- BN %*% bNh
  beta <- matrix(0, df, 1)
  hb <- rnorm(length(bN))

  beta_new <- t(chol(BN)) %*% hb + bN
  beta[as.logical(delta)] <- beta_new
  return(list(delta=delta, beta=beta))
}
