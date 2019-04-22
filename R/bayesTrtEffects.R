#' Bayesian Estimation of Treatment Effects in Panel Setting
#'
#' \code{bayesTrtEffects} performs an estimation of treatment effects on panel
#' structured data. The dataset may be unbalanced, but needs to fullfill
#' certain conditions to be properly useable. The function offers
#' some model choices, which provide the framework in which the treatment effect
#' is to be estimated.
#'
#' @param base.mat a data frame or matrix file capturing the baseline (at paneltime 0)
#' @param panel.mat a data frame file. Needs to be already correctly structured.
#' @param type contains one of currently one possible Strings to define the model under which
#' coefficients and treatment effects are computed ("SF" - Shared Factor Model)
#' @param covars contains a character vector, naming the covariates of the dataset.
#' @param model.name String which names the model, defaults to "model1"
#' @param sort.data boolean which indicates whether the baseline data and panel data
#' should be organised by treatment and panel time. If the dataset is not yet sorted
#' this option may sort data in correct format.
#' @param mcmc.control a list of items controlling the mcmc parameters, see examples
#' @param control a list with some general control parameters, the most important is
#' control$sort.data which sorts the data for ids, treatment and panel times.
#'
#' There are 2 main objects necessary for computation of treatment effects, base.mat and panel.mat.
#' base.mat contains a dataframe at the baseline (panel time 0) with specific features for latent utility computation.
#' panel.mat is a dataframe containing information at panel times, including features that influence
#' the panel outcomes. The datasets have to have certain properties.
#' See example dataset for more information on how to prepare the dataset
#' accordingly. For the datasets, the function first sets a lot of parameters
#' influencing estimation and MCMC sampling process. For the time being this
#' package includes the functionality of the Shared Factor Model.
#' It may also be of interest to look at \code{\link{select_trt_sharedfac}} to
#' see how the function works.
#'
#' @return For an adequate input file, \code{bayesTrtEffects} returns output a list
#' object containing all coefficients of the model.
#'
#' @export

bayesTrtEffects <- function(base.mat, panel.mat, type='SF', covars = NULL, model.name = "model1",
                            mcmc.control = list(burnin = 1000, select = 500, M = 1000),
                            control = list(fix.alpha=FALSE, fix.beta=FALSE, fix.sigma=FALSE,
                                           fix.f=FALSE, sort.data=FALSE)){

  data.list <- readPanelUb(base.mat, panel.mat, type = type, covars = covars,
                           name = model.name, sort.data = control$sort.data)
  data <- data.list$data; model <- data.list$model;
  rm("data.list")
  Tmax <- dim(data$start_x)[1]
  # Define variables subject to selection and starting values for indicators
  # always keep intercept
  model$deltax_fix <- c(1, data$covx_fix)
  model$deltay_fix <- c(1, data$covy_fix[1:(model$dys-1)], 0, data$covy_fix[,1:(model$dys-1)])
  if (model$type == "SF"||model$type == "SFimp"){
    model$deltax_fix <- c(model$deltax_fix, 1)
    model$deltay_fix <- c(model$deltay_fix, rep(1, 2*Tmax))  # for SF model: no selection on factors
  }
  prior <- list()
  prior_var_sel <- 5    # higher values lead to smaller inclusion probabilities !!! maybe move this into the parameters of the function, with default value = 5
  prior$m <- 0.001
  prior_var_fix <- 0.1
  prior$type <- 'conjugate'
  prior$par$invB0 <- c(prior$m, 1/prior_var_sel * (1 - model$deltay_fix[2:model$dy]) + (1/prior_var_fix) * model$deltay_fix[2:model$dy])
  prior$par$invA0 <- 1/prior_var_sel * rep(1, model$dx)
  prior$model$ax <- 1
  prior$model$bx <- 1
  prior$model$ay <- 1
  prior$model$by <- 1
  # Prior for variance covariance parameters
  if (model$type == "SRF"){
    prior$lnsig$c0 <- matrix(0, Tmax,2)
    prior$lnsig$C0inv <- matrix(1,Tmax,2)
    prior$rho$C0inv <- matrix(1,Tmax,2)
    prior$rho$c0 <- matrix(0,Tmax,2)
    prior.invL0 <- matrix(1,(2*Tmax),1)
  }
  if (model$type == "SRI"){
    prior$lnsig$c0 <- matrix(0,Tmax,2)
    prior$lnsig$C0inv <- matrix(1,Tmax,2)
    prior$rho$C0inv <- matrix(1,Tmax,2)
    prior$rho$c0 <- matrix(0,Tmax,2)
    prior$d <- 3
    prior$D <- 1
  }
  if (model$type== "SF" | model$type == "SFimp"){
    prior$par$invA0 <- c(prior$par$invA0,1)
    prior$invL0 <- rep(1,(Tmax*2))
    prior$s0 <- matrix(0,Tmax,2)
    prior$S0 <- matrix(0,Tmax,2)
    prior$model$al <- 1
    prior$model$bl <- 1
  }

  if (all(model$type != c("SRF", "SRI", "SF", "SFimp"))){
    stop('Model type not correctly specified')
  }

  mcmc <- list()
  mcmc$burnin <- mcmc.control$burnin
  mcmc$M <- mcmc.control$M
  mcmc$start_select <- mcmc.control$select     # iteration at which variables are susceptible to selection
  if (mcmc$start_select < mcmc$burnin){
    model$name <- paste0(model$name, '_varsel')
  }

  starttime <- Sys.time()
  mcmc$start <- list()

  if (any(model$type == c('SRF', 'SRI'))){
    mcmc$start$deltax <- rep(1,model$dx)
    mcmc$start$deltay <- rep(1,model$dy)
    mcmc$start$lnsig <- mcmc$start$rho <- matrix(0, Tmax, 2)
    mcmc_select <- select_trt_switchreg(data, model, prior, mcmc)
  }else{
    mcmc$start$deltax <- rep(1,model$dx + 1)
    mcmc$start$deltay <- rep(1,model$dy + 2*Tmax)

    if(model$type== 'SF'){
      mcmc_select <- select_trt_sharedfac(data, model, prior, mcmc, control)
    }else{
      mcmc_select <- select_trt_SFimp(data, model, prior, mcmc)
    }
  }
  trt_effect <- comp_trt_effects(mcmc_select = mcmc_select)
  return(list(markovchain = mcmc_select, effects = trt_effect))
}
