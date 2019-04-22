#' Modelling Treatment Effects via Shared Factor Model
#'
#' Internal. \code{select_trt_sharedfac} implements the shared
#' factor model. It is the core function of wrapper \code{\link{bayesTrtEffects}}
#' which passes its parameters on, concerning the MCMC process, model and data.
#'
#' @param data takes data object from the function bayesTrtEffects
#' @param model takes model parameters from function bayesTrtEffects
#' @param prior specified prior parameters
#' @param mcmc mcmc parameters specified by bayesTrtEffects, potentially subject to change
#'
#' @return Returns a MCMC selection list object containing all estimated coefficients and
#' other meassures.

select_trt_sharedfac <- function(data, model, prior, mcmc, control){

  # y...  observed data
  # y0...Potential outcomes 0  n x 1
  # y1...Potential outcomes 1  n x 1
  # x  ... treatments (0/1) n x 1
  # W ...Matrix  for fixed effects n x nw
  # T_i  ... nx1 vector with number of repeated obs. per subjectnu=ones(data.n,1)
  # indx0,indx1... indicates wheter x-obs. is 0or 1
  # indy  ...data.Tn x 1.... number of x observation to which y belongs
  # indy0, indy1.. data.Tn x 1 vectors .. indicates response belonngs to treatmt 0/1
  # beta Normal(b0,B0)  1x model$df contains(beta0,beta1,alpha,gamma)
  # dx,dfx ... number of free/fixed to 1 effects in model for x
  # dx,dfx ... number of free/fixed to 1 effects in model for y

  filetext <- paste0(model$name,'_',date())
  savefile <- paste0('save ', filetext, '_mcmc_select')

  # Create necessary variables
  Tmax <- max(data$Ti)            # maximum panel time
  y <- data$y                     # outcome vector (y for each subject and each panel time)

  dx <- model$dx + 1              # delta x (dx) represents number of covariates (potential) for selection model
  dfx <- sum(model$deltax_fix)    # dfx is the number of fixed effects (not subject to variable selection)
  dvx <- dx - dfx                 # dvx are variable covars.

  dfy <- sum(model$deltay_fix[1:model$dy])    # dy is the number of covars in y model (outcome model)
  dvy <- model$dy - dfy
  dyall <- model$dy + 2*Tmax                  # delta y all is the number of all covars for outcome model + the covariance parameters ??

  dfl <- sum(model$deltay_fix[(model$dy+1):dyall])  # delta fix l is the number of fixed effects of lambdas (latent factor weights) on each panel time ?
  dvl <- 2 * Tmax - dfl

  nmc <- mcmc$M + mcmc$burnin   # number of mcmc runs

  mcmc_select <- list(prior = prior, model = model, data = data, mcmc = mcmc)   # object to be returned

  mcmc_select$deltax <- matrix(0, nmc, dx)                # holds for each mcmc iteration the delta (variable selection indicator) of trt selection coefficients
  mcmc_select$alpha <- matrix(0, nmc, model$dx)           # holds the trt select model coefficient values
  mcmc_select$alpha_probit <- matrix(0, nmc, model$dx)    # ?
  mcmc_select$omega_alpha <- matrix(0, nmc, 1)            #

  mcmc_select$lambdax <- matrix(0, nmc, 1)                # holds the lambda parameter of trt selection model for nmc iterations

  mcmc_select$deltay <- matrix(0, nmc, model$dy)          # holds indicator values for each y model covar.
  mcmc_select$beta <- matrix(0, nmc, model$dy)            # holds the coefficients values of outcome models (both under treatment and without)
  mcmc_select$omega_beta <- matrix(0, nmc, 1)             # holds the estimated covariance values of outcome models

  mcmc_select$lambda <- array(0, dim = c(nmc,Tmax,2))       # holds lambda parameters for outcome models (all panel times and every treatment)
  mcmc_select$delta_lambda <- matrix(0, nmc, 2*Tmax)      # heolds indicator values for lambda parameters if included, for outcome models (lambda0 , lambda1) and lambdax is automaticly included?
  mcmc_select$omega_lambda <- matrix(0, nmc, 1)           # some correlation coefficient ?

  mcmc_select$sgma2 <- array(0, dim= c(nmc,Tmax,2))       # holds variance estimates


  #  compute posterior parameters which are independent of draws

  sn <- matrix(0, Tmax,2)                    # ?
  indt <- matrix(0, data$Tn, Tmax)           # index panel time (is 1 for the paneltime in which its active, 0 otherwise, therefore 6 columns here)

  for (t in 1:Tmax){
    indt[,t] <- (data$Tvec == t)
    sn[t,1] <- prior$s0[t,1] + sum(indt[,t] * data$indy0)/2
    sn[t,2] <- prior$s0[t,2] + sum(indt[,t] * data$indy1)/2
  }

  # indt <- matrix(as.vector(indt, mode='logical'), data$Tn, Tmax)
  invA0 <- diag(prior$par$invA0)                      # alpha start parameter values (prior)
  invB0 <- c(prior$par$invB0, prior$invL0)            # beta and lambda start parameter values (prior)

  # Starting values for  mcmc
  deltax <- mcmc$start$deltax       # starting deltas for x model (inclusion in probit model)
  deltay <- mcmc$start$deltay       # starting deltas for outcome model (16*2 + 6*2)

  # initial inclusion probabilities
  omega_alpha <- runif(1)             # omega_alpha shows inclusion probabilities for coefficients into the probit model (trt selection)
  omega_beta <- runif(1)
  omega_lambda <- runif(1)

  start.obj <- startingValues(model, data$x, y, data$Tn)    # get starting values for alpha_nu, beta param
  alphav <- start.obj$alphav; beta0 <- start.obj$beta0; res_var <- start.obj$res_var;
  rm(start.obj)
  # selection model
  mu_xstar <- data$Wx %*% alphav    # mean value of x parameters for given data (Wx should be feature vector of subject before treatment)
  xstar <- drawUtility(data$x, mu_xstar , 1)     # given parameters draw a latent utility x'star from
  lambdax <- 1          # initial lambda_x value

  # outcome model
  muy_fix <- model$W %*% beta0        # first mean outcome (estimates based on data and initial beta coeffs)

  # latent factors
  lambda <- matrix(0, 2*Tmax, 1)        # factor loadings

  f <- matrix(rnorm(data$n), data$n, 1)    # draw shared factors for each subject from N(0,1)
  fy <- f[data$indy]                     # multiply drawn factors to each panel outcome of subjects

  Fact <- fy * indt

  ny0 <- sum(data$indy0)                 # number of panel outcomes observed without trt
  ny1 <- data$Tn - ny0                   # number of panel outcomes under trt

  F01 <- Matrix::Matrix(0, ny0, Tmax)            # Sparse Zero matrix (a 1 for which panel time will be true)??
  F10 <- Matrix::Matrix(0, ny1, Tmax)

  Wfupper <- cbind(Fact[1:ny0,], F01)
  Wflower <- cbind(F10 , Fact[(ny0+1):data$Tn,])
  Wf <- rbind(Wfupper, Wflower)

  fcontr <- tcrossprod(Wf,t(lambda) )              # the factors drawn for each subject times their loadings lambda

  # create design matrix with latent factors
  Wy <- cbind(model$W, Wf)

  ###########################################################################
  # Run MCMC

  isel <- 0                   # index, is 1 as soon as selection process starts
  start_time <- Sys.time()
  print("starting MCMC process...")
  for (imc in 1:(mcmc$burnin + mcmc$M)){
    if (imc == mcmc$start_select)  print(paste0("starting selection of covariates at iteration ", imc))
    if (imc == mcmc$burnin) print(paste0("end of burnin phase at iteration ", imc))
    if (imc == mcmc$burnin + ceiling(mcmc$M/2)) print(paste0("iteration ", (mcmc$burnin + ceiling(mcmc$M/2)), " of ", (mcmc$burnin + mcmc$M), " reached."))

    # STEP I : Sample the idio-syncratic variances
    resy <- y - muy_fix                   # outcome residual vector
    epsy2 <- (resy - fcontr)^2              # residuals of outcome minus the factor loading squared equals residuals of model (11) and (12), the pure error ?

    error.obj <- drawErrorVariance(epsy2,indt, data$indy0, sn, prior, control$fix.sigma)  # draw variance values for idiosyncratic errors e Sigma_1 and Sigma_2
    sgma2 <- error.obj$sgma2; var_yt <- error.obj$var_yt;  # var_yt keeps variance estimates of every observed outcome value (169539)

    # STEP II: Sample the latent factors
    resx <- xstar - mu_xstar
    f <- drawSharedFactor(resx,lambdax,resy,sgma2,lambda,data$start_x,data$nx,data$start_y,data$Tn, control$fix.f)
    Fact <- rep(f[data$indy],Tmax) * indt

    Wfupper <- cbind(Fact[1:ny0,], F01)
    Wflower <- cbind(F10 , Fact[(ny0+1):data$Tn,])
    Wf <- rbind(Wfupper, Wflower)

    # STEP III: Sample the latent utilities
    xstar <- drawUtility(data$x, mu_xstar + lambdax*f, 1)

    # STEP IV Sample coefficients of treatment equation
    if ((imc > mcmc$start_select) & (sum(!model$deltax_fix) > 0)) isel <- 1

    Wx <- cbind(data$Wx,f)

    # y = xstar; X=Wx; delta=deltax; deltafix=model$deltax_fix; omega=omega_alpha; invA0; isel;

    indic.obj <- drawAlphaIndicesSF(xstar, Wx, deltax, model$deltax_fix, omega_alpha, invA0, isel, control$fix.alpha)
    deltax <- indic.obj$delta; alpha <- indic.obj$alpha;

    alphav <- alpha[1:model$dx]
    lambdax <- alpha[model$dx+1]
    mu_xstar <- tcrossprod(data$Wx,t(alphav))

    # STEP V  Sample indicators and fixed effects
    Wy[,(model$dy+1):dyall] <- Wf

    # y=y;X=Wy;sgma2=var_yt;delta=deltay;deltafix=model$deltay_fix;omega1=omega_beta;omega2=omega_lambda;dy=model$dy;invB0=invB0;isel;

    indbeta.obj <- drawBetaLambdaIndicesSF(y, Wy, var_yt, deltay, model$deltay_fix, omega_beta, omega_lambda,
                                          model$dy, invB0, isel)
    deltay <- indbeta.obj$delta; betav <- indbeta.obj$beta;

    beta <- betav[1:model$dy]
    muy_fix <- tcrossprod(model$W, t(beta))
    lambda <- betav[(model$dy + 1):dyall]
    # lambda <- c(0.6, 0.6, 0.5, 0.5, -0.6, -0.6, -0.5, -0.5) # delete afterwards

    fcontr <- Wf %*% lambda

    #signswitch
    rs <- sign(runif(1, -0.5, 0.5))
    lambda <- lambda*rs
    lambdax <- lambdax*rs

    if  (isel==1){
    # Update mixture weights
    dx1 <- sum(deltax==1)-dfx
    omega_alpha <- rbeta(1, prior$model$ax + dx1, prior$model$bx+dvx-dx1)

    dy1 <- sum(deltay[1:model$dy]==1) - dfy
    omega_beta <- rbeta(1, prior$model$ay + dy1,prior$model$by + dvy - dy1)

    dl1 <- sum(deltay[(model$dy+1):dyall]==1) - dfl
    omega_lambda <- rbeta(1, prior$model$al+dl1, prior$model$bl+dvl-dl1)
    }

  # store sampled values
  mcmc_select$deltax[imc,] <- deltax
  mcmc_select$alpha_probit[imc,] <- alphav / sqrt(lambdax^2 + 1)
  mcmc_select$alpha[imc,] <- alphav
  mcmc_select$lambdax[imc] <- lambdax
  mcmc_select$omega_alpha[imc] <- omega_alpha

  mcmc_select$beta[imc,] <- beta
  mcmc_select$deltay[imc,] <- deltay[1:model$dy]
  mcmc_select$omega_beta[imc] <- omega_beta

  mcmc_select$lambda[imc,,] <- matrix(as.vector(lambda), Tmax,2)
  mcmc_select$delta_lambda[imc,] <- deltay[(model$dy+1):dyall]
  mcmc_select$omega_lambda[imc] <- omega_lambda
  mcmc_select$sgma2[imc, , ] <- sgma2
  }

  mcmc_select$etime <- Sys.time() - start_time  # registers passed time for estimation of model coefficients with x mcmc iterations
  savefile <- paste('save ', filetext, ' mcmc_select' )
  return(mcmc_select)
}
