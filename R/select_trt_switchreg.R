#' 
#' Run treatment computations with switching regression model
#'
#' Internal.


select_trt_switchreg <- function(data, model, prior, mcmc){
  Tmax <- max(data$Ti)
  y <- data$y
  dn <- c(0,0)
  if (model$type == 'SRF'){
    indt <- matrix(0, nrow=data$Tn, ncol=Tmax)
    for (t in 1:Tmax){
      indt[,t] <- data$Tvec == t
      }
  }else{
    dn[1] <- prior$d + (data$nx0 / 2)
    dn[2] <- prior$d + (data$nx1 / 2)
  }
  dfx <- sum(model$deltax_fix)             # number of fixed delta coefs
  dvx <- model$dx - dfx                    # variable(= free - fixed) number of effects (for treatment) in the model?

  dfy <- sum(model$deltay_fix)
  dvy <- model$dy - dfy                    # number of variable effects (free - fixed) in the model for y (= earnings)?
  dall <- model$dx + model$dy              # number of all free variables in the model

  indx1 <- as.logical(data$x)              # index for mothers which chose treatment
  pri_invB0 <- c(prior$par$invA0, prior$par$invB0)

  nmc <- mcmc$M + mcmc$burnin
  mcmc_select <- list(prior = prior, model = model, data = data, mcmc = mcmc)
  mcmc_select$deltax <- mcmc_select$alpha <- matrix(0, nrow=nmc, ncol=model$dx)
  mcmc_select$acc_alpha <- mcmc_select$omega_alpha <- matrix(0, nmc, 1)
  mcmc_select$deltay <- mcmc_select$beta <- matrix(0, nmc, model$dy)
  mcmc_select$omega_beta <- matrix(0, nmc, 1)
  mcmc_select$sgma2 <- mcmc_select$rho <- mcmc_select$accrs <-  array(0, dim=c(nmc, Tmax, 2))  # in these arrays the mcmc draws will be saved (6x2 matrix for each MCMC run)
  if (model$type == "SRF"){
    mcmc_select$lambda <- array(0, dim = c(nmc, Tmax, 2))             # 3 dimensional array for lambdas for every coef draw, panel times, trt x
  }else{
    mcmc_select$D <- matrix(0, nmc, 2)
  }# for the Random Intercept model we assume no cov. between the factors, therefore we need only draws and trtm x

  # Starting values for mcmc
  deltax <- mcmc$start$deltax
  deltay <- mcmc$start$deltay
  delta <- c(deltax, deltay)
  delta_fix <- c(model$deltax_fix, model$deltay_fix)
  # inclusion probabilities
  omega_alpha <-  runif(1)
  omega_beta <-  runif(1)
  start.obj <- startingValues(model, data$x, data$y, data$Tn)
  alphav <- as.vector(start.obj$alphav); beta0 <- start.obj$beta0; res_var <- start.obj$res_var;
  # variable selection model
  mu_xst <- model$Wx %*% alphav
  xst <- drawUtility(data$x, mu_xst, 1)
  muy_fix <- model$W %*% beta0
  a <-  runif(1)
  sgma <- sqrt(a*res_var)* matrix(1, Tmax, 2)
  resy <- y - muy_fix
  epsy <- resy
  rho <- mcmc$start$rho
  if (model$type == "SRF"){
    D <- matrix(1, 2 ,1)
    lambda <- matrix(0, Tmax, 2)
  }else{
    D <- (1-a) * res_var * matrix(1,2,1)
    lambda <- matrix(1, Tmax, 2)
  }

  isel <- 0

  ############ Run MCMC #######################################################

  for (imc in 1:(mcmc$burnin + mcmc$M)){                       # burnin period + wanted draws M
    if (mod(imc,round((mcmc$burnin + mcmc$M)/4)) == 0){        # about 4 times throughout the process save intermediate result, display iteration
      print(imc)
    }
    if (imc > mcmc$start_select && sum(!(model$deltay_fix)) > 0 ) isel <- 1

    # STEP I: draw mixture weights
    if (isel==1){
      # Update mixture weights

      dx1 <- sum(deltax == 1) - dfx
      omega_alpha <- rbeta(n = 1, prior$model$ax + dx1, prior$model$bx + dvx - dx1)
      dy1 <- sum(deltay == 1) - dfy
      omega_beta <- rbeta(n = 1, prior$model$ay + dy1, prior$model$by + dvy - dy1)
    }
    # Draw utility conditional on latent factor/random intercept
    # in case of random intercept the conditional variance is always 1 and rho is 0
    # the latent factor model can have covariance matrix != I
    mxc.obj <- condmoms_util(mu_xst, epsy, data$start_x, data$nx, data$start_y, sgma, rho)
    mxc <- mxc.obj$mxc; var_xc <- mxc.obj$var_xc;
    xst <- drawUtility(data$x, mxc, sqrt(var_xc))

    # STEP IIa Sample indicators and fixed effects  marginalising over
    # latent factors or random intercept

    delta.obj <- draw_indic_coeff_switchreg(xst, data$y, data$Wx, model$W,
                                            data$start_x, data$nx, data$start_y,
                                            sgma, rho, lambda, D, delta, delta_fix,
                                            omega_alpha, omega_beta, pri_invB0, isel)
    delta <- delta.obj$delta; coeff <- delta.obj$beta

    deltax <- delta[1:model$dx]
    alphav <- coeff[1:model$dx]
    mu_xst <- data$Wx %*% alphav
    resx <- xst - mu_xst

    deltay <- delta[(model$dx+1):dall]
    beta <- coeff[(model$dx+1):dall]
    muy_fix <- model$W %*% beta
    resy <- y - muy_fix

    # STEP IIb : draw the random factors/random intercepts

    f <- draw_factor_switchreg(resy, resx, data$start_x, data$nx, data$start_y, sgma, rho, lambda, D)
    fyi <- f[data$indy]

    if (model$type == 'SRF'){
      F0 <- do.call(cbind, replicate(Tmax, fyi*data$indy0, simplify=FALSE)) * indt
      F1 <- do.call(cbind, replicate(Tmax, fyi * data$indy1, simplify=FALSE)) * indt

      # STEP III  : draw the factor loadings of the random factor
      Xf <- cbind(F0, F1)
      lambda_new <- draw_loadings(resy, Xf, resx, data.start_x, data.nx, data.start_y, sgma, rho, prior)
      fcontr <- Xf %*% lambda_new
      lambda <- matrix(as.vector(lambda_new), Tmax, 2)
      rs <- do.call(rbind, replicate(Tmax, sign(runif(2)-0.5), simplify=FALSE))
      lambda <- lambda * rs

      epsy <- resy - fcontr
    }else{
        Dn <- c(0,0)
        # STEP III  : draw the random intercept variances
        Dn[1] <- prior$D + (t(f[!indx1]) %*% f[!indx1] / 2)
        Dn[2] <- prior$D + (t(f[indx1]) %*% f[indx1] / 2)
        D <- 1 / (rgamma(2, shape = dn, scale = 1/Dn))
        epsy <- resy - fyi
    }

    # STEP IV : draw rho and sigma
    covparms.obj <- draw_cov_parms(epsy, mu_xst, data$start_x, data$nx, data$start_y, sgma,
                                   rho, prior)
    sgma <- covparms.obj$sgma; rho <- covparms.obj$rho ; acc <- covparms.obj$acc;

    # Save the draws

    mcmc_select$deltax[imc,] <- deltax
    mcmc_select$alpha[imc,] <- alphav
    mcmc_select$omega_alpha[imc] <- omega_alpha
    # mcmc_select$acc_alpha(imc) <- acc_alpha

    mcmc_select$beta[imc,] <- beta
    mcmc_select$deltay[imc,] <- deltay
    mcmc_select$omega_beta[imc] <- omega_beta

    mcmc_select$sgma2[imc, , ] <- sgma^2
    mcmc_select$rho[imc, , ] <- rho
    mcmc_select$accrs[imc, , ] <- acc

    if (model$type == "SRF"){
      mcmc_select$lambda[imc, , ] <- lambda
    }else{
      mcmc_select$D[imc, ] <- D
    }
  }
  endtime <- Sys.time()
  mcmc_select$etime <- endtime - starttime
  return(mcmc_select)
}
