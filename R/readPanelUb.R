#' Read Data file and transform
#'
#' [Internal] \code{readPanelUb} takes as arguments the baseline matrix and panel matrix
#' seperately and transforms the (unbalanced) panel data in useable format.
#' The matrices should xalready be in the suggested format for this
#' function. (See example dataset)
#'
#' Base matrices have 1 subject in each row and a unique corresponding ID.
#' It is used to capture the base features of a subject, which are used for
#' estimating selection model effects. Also included are columns with
#' the maximum panel time observed, treatment indicator (0/1).
#' The panel matrix represents an expansion on panel times of the bases
#' matrix with some additional information. Covariates for the treatment
#' model can be included, because they can vary over panel times,
#' yielding the heterogenous effect of treatment. Panel matrix must
#' include the columns ID, panel time, target outcome, some covariates
#' for the models and panel dummy columns (if no panel dummy are included
#' you can specify p.dummy = TRUE)
#' Strictly keep to this structure, the function provides some quality
#' checks, but may prove insufficient for some wierd structured individual
#' chases and proceed to calculate, potentially giving bad quality estimate
#' results. The base and panel matrices should also be in ascending order
#' by subject panel times, and by maximum panel time, meaning subjects with
#' max panel time 3 are listed before subjects with maximum panel time
#' 4, 5, etc. This is necessary to use quicker evaluation functions for the
#' estimation procedures.
#' Also at the moment the base matrix needs to have the subject ID
#' in the first column, panel times of each subject in second column, the
#' treatment in the third and all the covariates for the selection model in
#' subsequent columns.
#'
#' @param infile a list of base matrix and panel matrix
#' @param type contains a string of possible model types, must be one of: "SF"
#' @param covars contains additional information about covariate names (optional)
#' @param name contains name of the model, which will be saved
#' @param sort_data control parameter passed down in funcional hierarchy to
#' determine whether observations should be sorted for by ID and panel time
#'
#' @return A data object for further computation, containing of the modified
#' dataset and a model object, which contains covariate information
#' amongst other specifics.
#'

readPanelUb <- function(base_mat, panel_mat, type, covars, name, sort_data, control_test, data_name){

  data <- list(name = data_name)
  n <- dim(base_mat)[1]
  k <- dim(base_mat)[2]
  # n ... number of subjects
  # k ... number of subject features at baseline panel time t = 0 including ID, trt, and maxPaneltime
  id <- base_mat[,1]
  Ti <- base_mat[,2]
  x <- base_mat[,3]
  # check for unique IDs
  if (length(id) != length(unique(id))) stop("IDs (column 1) seem not to be unique.")
  # ensure that treatment is binary 0/1
  if ((max(x) - min(x))!=1) stop("Treatment specifications might be wrong. Has to be coded
                                 binary, the higher index being the treated group,
                                 the lower index being the untreated group.")
  if (min(x) != 0) x <- x - min(x)
  # displace the distribution of treatment values to be 0 and 1
  Wx_dim <- k-2       # number of regressors used: columns without columns 1,2,3 + intercept (dimension Wx)
  # read and check panel matrix
  # panel matrix consists of subject features at panel times (variables that may vary within panel,
  # are recorded at each panel time, and panel time dummies
  k_pan <- dim(panel_mat)[2]
  W_dim <- k_pan - 3                 # dimension W = columns - id - panelT - outcome
  Tn <- sum(Ti)           # sum of panel observations across all subjects

  # create the finally ordered data set
  data$x <- rep(0, n)
  data$Ti <- rep(0, n)
  data$indy <- rep(0, Tn)
  data$Tvec <- rep(0, Tn)
  data$y <- rep(0, Tn)
  data$Wi <- matrix(0, Tn, W_dim)
  data$n <- n
  data$ntrt <- length(unique(x))
  data$Tn <- Tn
  data$dw <- W_dim
  data$dx <- Wx_dim
  Tmin <- min(Ti)
  Tmax <- max(Ti)
  # create starting indices for treatment group with t panel observations
  # number of observations for certain panel time and treatments
  data$start_x <- matrix(0,Tmax, data$ntrt)
  data$start_y <- matrix(0,Tmax, data$ntrt)
  data$nx <- matrix(0,Tmax, data$ntrt)
  data$Wx <- matrix(NA, n, k-2)

  T_end <- cumsum(Ti)                    # indices where individual panel observed ends
  T_start <- c(1, T_end[1:data$n-1] + 1) # indices where individual panel starts (in panel matrix)

  # creates new help columns trt and Ti in panel matrix for sorting, which will be dropped afterwards

  panel_mat$trt <- 0
  panel_mat$trt[panel_mat$ID %in% base_mat$ID[which(base_mat$trt==1)]] <- 1
  panel_mat$Ti <- Tmin
  for (i in Tmin:Tmax){
    panel_mat$Ti[panel_mat$ID %in% base_mat$ID[which(base_mat$Ti==i)]] <- i
  }
  if (sort_data){
    base_mat <- base_mat[order(base_mat$trt, base_mat$Ti, base_mat$ID), ]
    panel_mat <- panel_mat[order(panel_mat$trt, panel_mat$Ti, panel_mat$ID, panel_mat$panelT), ]
  }
  index_concat <- vector("numeric")
  for (i in 0:(data$ntrt-1)){
    for (t in Tmin:Tmax){
      base_ind <- which(base_mat$trt==i & base_mat$Ti == t)
      panel_ind <- which(panel_mat$trt == i & panel_mat$Ti == t)
      index_concat <- c(index_concat, rep(base_ind, each = t))
      data$start_x[t,i+1] <- min(base_ind)
      data$start_y[t,i+1] <- min(panel_ind)
      data$nx[t,i+1] <- length(base_ind)
    }
  }
  data$x <- base_mat$trt
  data$Tvec <- panel_mat$panelT
  data$y <- panel_mat$y
  data$Wx <- as.matrix(cbind(1, base_mat[,4:k]))
  data$Wi <- as.matrix(panel_mat[,4:k_pan])
  data$indy <- index_concat
  data$Ti <- base_mat$Ti
  id <- base_mat[,1]
  upper_ind <- 0
  if (control_test){
    for (i in 1:n){
      Tiy <- Ti[i]
      lower_ind <- upper_ind + 1
      upper_ind <- upper_ind + Tiy
      id_panel <- panel_mat[lower_ind:upper_ind, 1]
      if (!all(id_panel==id[i])){ stop("subjects are not in the same order as in base
                                            file or observations are missing") }
      tpanel <- panel_mat[lower_ind:upper_ind, 2]
      if (!all(tpanel == 1:Tiy)) stop("subjects are not in the same order as in the base file or obs. missing")
    }
  }
  data$nx0 <- sum(1-data$x)
  data$nx1 <- sum(data$x)
  data$indy0 <- as.logical(panel_mat$trt == 0)
  data$indy1 <- as.logical(panel_mat$trt == 1)
  data$ny0 <- sum(data$indy0 == 1)
  data$ny1 <- sum(data$indy1 == 1)
  panel_mat$trt <- panel_mat$Ti <- NULL
  cov_x <- names(base_mat)[4:dim(base_mat)[2]]
  cov_y <- names(panel_mat)[4:dim(panel_mat)[2]]
  # read covariate names if available
  covarsnames <- c("x_fix", "y_fix", "y_common", "y_nb", "y_np")
  if (!is.null(covars) & all(covarsnames %in% names(covars))){               # checks correctness of covariates parameters
    names(covars) <- paste0("cov_", names(covars))
    data <- c(data, list(cov_x=cov_x, cov_y=cov_y), covars)
  }else{
    message("argument covars is NULL, all covariates will be subject to variable selection.")
    data <- c(data, list(cov_x_fix = matrix(0, 1, data$dx-1),
                         cov_y_fix = matrix(0, 1, data$dw), cov_y_common = matrix(0, 1, data$dw)))
  }
  model <- list(type = type, name = paste0(data$name,"-model: ",type))
  model_object <- make_regmat_unbal(data, model)
  # build covariate names of models for output
  # add intercept to covariate names
  x <- c("intercept", cov_x)
  y <- c("intercept", cov_y)
  model_object$model$covnames_trt <- paste0("x: ", x)
  covars_y_select <- y[!c(0, covars$cov_y_common)]  # Add 0 for the intercept, which is always selected, never common
  covars_y_common <- y[as.logical(c(0, covars$cov_y_common))]
  model_object$model$covnames_resp <- c(paste0("trt:0 ", covars_y_select),
                                        paste0("trt:1 ", covars_y_select), covars_y_common)
  return(model_object)
}
