#' Compute treatment effects
#'
#' [Internal Function] computes treatment effects from the coefficients
#' sampled in the MCMC process.
#'

comp_trt_effects <- function(mcmc_select, thin = 1){
  data <- mcmc_select$data
  model <- mcmc_select$model
  file.name <- paste(mcmc_select$data$name, "_", mcmc_select$model$type)
  Tmax <- max(data$Ti)
  a <- mcmc_select$mcmc$burnin + 1
  e <- mcmc_select$mcmc$burnin + mcmc_select$mcmc$M
  indmc <- seq(a,e,by=thin)
  nmc <- length(indmc)
  alpha <- mcmc_select$alpha[indmc,]
  lambdax <- mcmc_select$lambdax[indmc]
  # sgnlx <- sign(lambdax)

  beta <- mcmc_select$beta[indmc, ]
  lambda <- mcmc_select$lambda[indmc, , ]

  #  computing average treatment effects

  k <- model$dys-Tmax + 1   # number of effects of covariates
  kall <- model$dys

  # ATE  are computed for comparison
  ate_log <- matrix(0,nmc,Tmax)
  ate_perc <- matrix(0, nmc,Tmax)
  y0 <-  matrix(0,nmc,Tmax)

  W <- model$W[data$Tvec==1, 1:k]

  for (imc in 1:nmc){
    theta <- beta[imc,(model$dys+(1:model$dys))]
    diff_t1 <- W %*% theta[1:k]
    fac_t1 <- exp(diff_t1)
    panel_eff<- t(c(0, theta[(k+1):kall]))
    ate_log[imc,] <- mean(diff_t1) + panel_eff
    ate_perc[imc,] <- mean(fac_t1)*exp(panel_eff)-1
    y0[imc,] <- mean(t(W*beta[imc,1:k])) + panel_eff
  }
  mat <- cbind(colMeans(ate_log), apply(ate_log,2,sd),  colMeans(ate_perc), apply(ate_perc,2,sd))
  mat
  varnames <- as.character((1:Tmax))
  rm(W)

  # compute treatment effect on untreated

  Wx <- model$Wx[1:data$nx0, ]
  indy <- data$indy0 & (data$Tvec == 1)
  W <- model$W[indy, 1:k]   # only take rows where t=1

  tu_log <- matrix(0, nmc,Tmax)
  tu_perc <- matrix(0,nmc,Tmax)

  # difference in lambda from trt to untrt
  diff_lambda <- lambda[, , 2] - lambda[ , , 1]

  for (imc in 1:nmc){
    theta <- beta[imc,model$dys + (1:model$dys)]
    diff_t1 <- W %*% theta[1:k]
    mux <- Wx %*% alpha[imc,]
    sx <- sqrt(lambdax[imc]^2+1)    # lambda^2 = variance of errors in SF model
    a <- mux/sx
    h0 <- dnorm(-a)/pnorm(-a)
    panel_eff <- c(0, theta[(k+1):kall])
    tu_logh <- do.call(cbind, replicate(Tmax, diff_t1, simplify=FALSE)) - ((h0 %*% lambdax[imc] %*% diff_lambda[imc,]) / sx)
    tu_log[imc,] <-  mean(tu_logh) + panel_eff
    tu_perc[imc,] <- mean(exp(tu_logh))*(exp(panel_eff)-1)
  }

  mattu <- cbind(colMeans(tu_log), apply(tu_log,2,sd),  colMeans(tu_perc), apply(tu_perc,2,sd))
  varnames <- as.character(1:Tmax)

  # compute treatment effect on treated

  Wx <- model$Wx[(data$nx0+1):data$n, ]
  indy <- data$indy1 & (data$Tvec==1)
  W <- model$W[indy,1:k]

  tt_log <- matrix(0,nmc,Tmax)
  tt_perc <- matrix(0,nmc,Tmax)
  for (imc in 1:nmc){
    theta <- beta[imc, model$dys + (1:model$dys)]
    diff_t1 <- W %*% theta[1:k]
    mux <- Wx %*% alpha[imc,]
    sx <- sqrt(lambdax[imc]^2 + 1)
    a <- mux/sx
    h1 <- dnorm(-a)/pnorm(a)
    panel_eff <-  c(0, theta[(k+1):kall])
    tt_logh <- do.call(cbind, replicate(Tmax, diff_t1, simplify=FALSE)) + ((h1 %*% lambdax[imc] %*% diff_lambda[imc,]) / sx)
    tt_log[imc,] <- mean(tt_logh) + panel_eff
    tt_perc[imc,] <- mean(exp(tt_logh))*exp(panel_eff) - 1
  }
  mattt<- cbind(colMeans(tt_log), apply(tt_log,2,sd),  colMeans(tt_perc), apply(tt_perc,2,sd))
  results <- list(ate = mat, tu = mattu, tt = mattt)
  return(results)
}
