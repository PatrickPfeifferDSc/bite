#'
#' Draw indices for selection upon sigma and rho parameters for SRF/SRI models
#'
#' Internal.

draw_indic_coeff_switchreg <- function(x, y, Wx, Wy, start_x, nx, start_y, sgma,rho, lambda, D, delta, deltafix,
                           omega1, omega2, pri_invB0, isel){

  df <- length(delta)
  nd <- df - sum(deltafix)
  dx <- dim(Wx)[2]

  # help quantities for posterior moments
  W.obj <- postmom_help_reg(x, y, Wx, Wy, start_x, nx, start_y, sgma, rho, lambda, D)
  WSW <- W.obj$WSW; WSy <- W.obj$WSy; ySy <- W.obj$ySy

  invB0 <- diag(pri_invB0)
  invBN <- invB0 + WSW

  bn.obj <- comp_logml(ySy,WSy,invBN,delta,invB0)
  lmlik_old <- bn.obj$lmarlik; BN <- bn.obj$BN; bNh<- bn.obj$bNh;

  if (isel == 1){
    idvar <- which((!deltafix)!=0)
    delta_ranord <- sample(1:nd, length(1:nd))
    for (i in 1:nd){
      j <- idvar[delta_ranord[i]]    # Selection of j-th not fixed effect
      # delta_new: change deltaj
      dj <- delta[j]
      delta_new <- delta
      delta_new[j] <- 1-dj
      if (j <= dx){     # whether reg soefs for x component or for y component
        lpri <- dj*log(omega1)+(1-dj)*(log(1-omega1))
        lpri_new <- (1-dj)*log(omega1)+dj*(log(1-omega1))
      }else{
        lpri <-dj*log(omega2)+(1-dj)*(log(1-omega2))
        lpri_new <- (1-dj)*log(omega2)+dj*(log(1-omega2))
      }
      lold <- lmlik_old+lpri

      loglik.obj <- comp_logml(ySy,WSy,invBN,delta_new,invB0)
      lmlik_new <- loglik.obj$lmarlik; BN_new <- loglik.obj$BN; bNh_new <- loglik.obj$bNh
      lnew <- lmlik_new+lpri_new

      # Draw new indicator delta(j)
      # compute posterior probabilities

      lh <- c(lnew,lold)
      maxl <- max(lh)
      prh <- exp(lh-maxl)
      postprob <- prh/sum(prh)
      # Update indicator
      if (runif(1)<= postprob[1]){
        delta[j] <- delta_new[j]
        lmlik_old <- lmlik_new
        bNh <- bNh_new
        BN <- BN_new
      }
    }
  }

  # Draw the regression effects
  bN <- BN %*% bNh
  beta <- matrix(0,df,1)
  hb <- matrix(rnorm(length(bN)), length(bN),1)

  beta_new <- (t(chol(BN)) %*% hb) + bN
  beta[as.logical(delta)] <- beta_new  # only delta indicated are updated? check whether this is a good solution or all are to be updated.
  return(list(delta=delta, beta=beta))
}
