#' Draw covariates for SRF/SRI models
#'
#' Internal.

draw_cov_parms <- function(epsy, mu_xst, start_x, nx,start_y, sgma_old, rho_old, pri){

  df <- 10
  Tmax <- dim(sgma_old)[1]
  ntrt <- dim(sgma_old)[2]
  sgma <- rho <- acc <- matrix(0, Tmax, ntrt)
  # Draw sigma

  for (j in 1:ntrt){
    # Proposal: maximization with  starting value theta_old
    sgmaj <- sgma_old[,j]
    rhoj <- rho_old[,j]
    start_xj <- start_x[,j]
    n_xj <- nx[,j]
    start_yj <- start_y[,j]
    tperm <- sample(1:Tmax, length(1:Tmax), FALSE)

    for (it in 1:Tmax){
      t <- tperm[it]
      thetat_old <- log(sgmaj[t])
      # Use sequential Nonlinear Programming (SQP-Algorithm) to find theta that minimizes LLH
      sqp <- solnp(pars = thetat_old, fun = negllik_jt, LB=-10, UB=10, j=j, t=t, epsy=epsy, mu_xst=mu_xst, start_x=start_xj,
                   nx=n_xj, start_y=start_yj, sgmaj=sgmaj, rhoj=rhoj) # starting value for theta = old value  (variance parameter)?
      mtheta <- sqp$pars
      hessian <- sqp$hessian
      print(hessian)

      ## Use sequential quadratic programming slsqp from package nloptr
      # qpres <- slsqp(x0 = thetat_old, fn = negllik_jt(j=j, t=t, epsy=epsy, mu_xst=mu_xst, start_x=start_x,
      #                                        nx=nx, start_y=start_y, sgmaj=sgmaj, rhoj=rhoj),
      #       lower = -10, upper = 10)
      # mtheta <- qpres$value

      # Generate proposal theta_new
      Vtheta <- 1/hessian
      h <- rt(1,df)
      thetat_new <- (t(sqrt(Vtheta)) %*% h) + mtheta
      sgmaj_new <- sgmaj
      sgmaj_new[t] <-  exp(thetat_new)

      # theta is Normal(mu_old, sigma_old) -> updated via Likelihood
      # Logposterior at theta_old
      lpri_old <- lnorm(thetat_old, pri$lnsig$c0[t,j], pri$lnsig$C0inv[t,j])
      lpost_old <- llikj(sgmaj,rhoj,j,epsy,mu_xst,start_xj,n_xj,start_yj) + lpri_old

      # Logposterior at theta_new
      lpri_new <- lnorm(thetat_new, pri$lnsig$c0[t,j], pri$lnsig$C0inv[t,j])
      lpost_new <- llikj(sgmaj_new,rhoj ,j, epsy, mu_xst, start_xj, n_xj, start_yj) + lpri_new

      # Log-transition density  and MH
      lq_new <- log(pt((thetat_new - mtheta)/sqrt(Vtheta),df=df))
      lq_old <- log(pt((thetat_old - mtheta)/sqrt(Vtheta),df=df))

      mh.obj <- mh_step(thetat_old, thetat_new, lpost_old, lpost_new, lq_new, lq_old)
      theta_tj <- mh.obj$theta;  acc_tj <- mh.obj$acc;
      sgmaj[t] <- exp(theta_tj[1])
    }
    sgma[,j] <- sgmaj
  }

# same steps for drawing correlation rho:

for (j in 1:ntrt){
  # Proposal: maximization with  starting value theta_old;
  sgmaj <- sgma[,j]
  rhoj <- rho_old[,j]
  start_xj <- start_x[,j]
  n_xj <- nx[,j]
  start_yj <- start_y[,j]
  tperm <- sample(1:Tmax, length(1:Tmax), FALSE)

  for (it in 1:Tmax){
    t <- tperm[it]
    thetat_old <- rhoj[t]
    bd_rhot <- sqrt(0.999 - sum(rhoj^2) + rhoj[t]^2)  #boundaries for rho_t

    # Use sequential Nonlinear Programming (SQP-Algorithm) to find theta that minimizes LLH
    sqp <- solnp(pars = thetat_old, fun = negllik_jt, LB=-bd_rhot, UB=bd_rhot, j=j, t=t, epsy=epsy, mu_xst=mu_xst, start_x=start_xj,
                 nx=n_xj, start_y=start_yj, sgmaj=sgmaj, rhoj=rhoj) # starting value for theta = old value  (variance parameter)?
    mtheta <- sqp$pars
    hessian <- sqp$hessian

    # Generate proposal theta_new
    Vtheta <- 1/hessian
    iconstr <- FALSE
    while (!iconstr){
      h <- rt(1, df=df)
      thetat_new <- sqrt(Vtheta)*h + mtheta
      iconstr <- (abs(thetat_new) < bd_rhot)
    }
    rhoj_new <- rhoj
    rhoj_new[t] <- thetat_new

    # Logposterior at theta_old

    lpri_old <- lnorm(thetat_old, pri$rho$c0[t,j], pri$rho$C0inv[t,j])
    lpost_old <- llikj(sgmaj, rhoj, j, epsy, mu_xst, start_xj, n_xj, start_yj) + lpri_old

    # Logposterior at theta_new
    lpri_new <- lnorm(thetat_new, pri$rho$c0[t,j], pri$rho$C0inv[t,j])
    lpost_new <- llikj(sgmaj, rhoj_new, j, epsy, mu_xst, start_xj, n_xj, start_yj) + lpri_new

    # Log-transition density  and MH
    lq_new <- log(pt((thetat_new-mtheta)/sqrt(Vtheta),df=df))
    lq_old <- log(pt((thetat_old-mtheta)/sqrt(Vtheta),df=df))

    mh.obj <- mh_step(thetat_old, thetat_new, lpost_old, lpost_new, lq_new, lq_old)
    theta_tj <- mh.obj$theta; acc_tj <- mh.obj$acc;

    rhoj[t] <- theta_tj
    acc[t,j] <- acc_tj
  }
  rho[,j] <- rhoj
}
  return(list(sgma=sgma, rho=rho, acc=acc))
}
