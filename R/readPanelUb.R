#' Read Data file and transform
#'
#' Internal. \code{readPanelUb} takes as arguments the baseline matrix and panel matrix
#' seperately and transforms the (unbalanced) panel data in useable format.
#' The matrices should xalready be in the suggested format for this
#' function. (See example dataset)l
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
#' max panel time 3 are listed before subjects with maximum panel time 4, 5
#' ,etc. This is necessary to use the quicker evaluation functions for the
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
#' @param sort.data control parameter passed down in funcional hierarchy to
#' determine whether observations should be sorted for by ID and panel time
#'
#' @return A data object for further computation, containing of the modified
#' dataset and a model object, which contains covariate information
#' amongst other specifics.
#'

readPanelUb <- function(base.mat, panel.mat, type, covars, name, sort.data, control.test=FALSE){

  data.list <- list()
  n <- dim(base.mat)[1];  k <- dim(base.mat)[2]
  # n ... number of subjects , k ... number of subject features at baseline panel time t = 0

  id <- base.mat[,1]
  Ti <- base.mat[,2]            # max. panel time of each subject
  x <- base.mat[,3]
  # ensure that treatment is binary 0/1
  if ((max(x) - min(x))!=1) stop("Treatment specifications might be wrong. Has to be coded
                                 binary, optimally 0/1")
  if (min(x) != 0) x <- x - min(x)
  # displace the distribution of treatment values to be (0/1) binary
  Wx.dim <- k-2            # number of regressors used: columns without columns 1 to 3
  # + intercept (dimension Wx)
  # read and check panel matrix
  # panel matrix consists of subject features at panel times (variables that may change over time,
  # are recorded at each panel time, and panel time dummies
  k.pan <- dim(panel.mat)[2]
  W.dim <- k.pan - 3                 # dimension W = columns - id - panelT - outcome
  Tn <- sum(Ti)           # sum of panel observations across all subjects

  # create the finally ordered data set

  data.list$x <- rep(0, n)
  data.list$Ti <- rep(0, n)
  data.list$indy <- data.list$Tvec <- data.list$y <- rep(0, Tn)
  data.list$Wi <- matrix(0, Tn, W.dim)
  data.list$Wx <- matrix(NA, ncol=k-2, nrow=n)
  data.list$n <- n
  data.list$ntrt <- 2
  data.list$Tn <- Tn
  data.list$dw <- W.dim
  data.list$dx <- Wx.dim
  Tmin <- min(Ti)
  Tmax <- max(Ti)

  ##### START sort data ##################

  # create starting indices for treatment group with t panel observations
  # number of observations for certain panel time and treatments
  data.list$start_x <- data.list$start_y <- data.list$nx <- matrix(0,Tmax, data.list$ntrt)

  T_end <- cumsum(Ti)                    # indices where individual panel observed ends
  T_start <- c(1, T_end[1:data.list$n-1] + 1) # indices where individual panel starts (in panel matrix)

  # Start indices for all vectors
  index.base <- 1
  index.panel <- 0

  # creates new help columns trt and Ti in panel matrix for sorting, which will be dropped afterwards

  panel.mat$trt <- 0
  panel.mat$trt[panel.mat$ID %in% base.mat$ID[which(base.mat$trt==1)]] <- 1
  panel.mat$Ti <- Tmin
  for (i in Tmin:Tmax){
    panel.mat$Ti[panel.mat$ID %in% base.mat$ID[which(base.mat$Ti==i)]] <- i
  }
  if (sort.data){
    base.mat <- base.mat[order(base.mat$trt, base.mat$Ti, base.mat$ID), ]
    panel.mat <- panel.mat[order(panel.mat$trt, panel.mat$Ti, panel.mat$ID, panel.mat$panelT), ]
  }
  for (i in 0:1){
    for (t in Tmin:Tmax){
      base.ind <- which(base.mat$trt==i & base.mat$Ti == t)
      panel.ind <- which(panel.mat$trt == i & panel.mat$Ti == t)
      data.list$start_x[t,i+1] <- min(base.ind)
      data.list$start_y[t,i+1] <- min(panel.ind)
      data.list$nx[t,i+1] <- length(base.ind)
    }
  }
  data.list$x <- base.mat$trt
  data.list$Tvec <- panel.mat$panelT
  data.list$y <- panel.mat$y
  data.list$Wx <- as.matrix(cbind(1, base.mat[,4:k]))
  data.list$Wi <- as.matrix(panel.mat[,4:k.pan])
  data.list$indy <- panel.mat[,1]
  data.list$Ti <- Ti

  id <- base.mat[,1]
  upper.ind <- 0
  if (control.test){
    for (i in 1:n){
      Tiy <- Ti[i]
      lower.ind <- upper.ind + 1
      upper.ind <- upper.ind + Tiy
      id.panel <- panel.mat[lower.ind:upper.ind, 1]
      if (!all(id.panel==id[i])){ stop("subjects are not in the same order as in base
                                            file or observations are missing") }
      tpanel <- panel.mat[lower.ind:upper.ind, 2]
      if (!all(tpanel == 1:Tiy)) stop("subjects are not in the same order as in the base file or obs. missing")
    }
  }
  panel.mat$trt <- panel.mat$Ti <- NULL
  data.list$nx0 <- sum(1-data.list$x)
  data.list$nx1 <- sum(data.list$x)
  data.list$indy0 <- as.logical(1-data.list$x[data.list$indy])
  data.list$indy1 <- as.logical(data.list$x[data.list$indy])
  data.list$ny0 <- sum(data.list$indy0 == 1)
  data.list$ny1 <- sum(data.list$indy1 == 1)

  # read covariate names if available
  covarsnames <- c("covx", "covx_fix", "covy", "covy_fix", "covy_common", "covy_nb", "covy_np")
  if (is.null(covars) || !(all(covarsnames %in% names(covars))) ){               # checks correctness of covariates parameters
    warning("Argument covars empty or NULL, number of fixed covariates will be set to 0.")
    data <- c(data.list, list(covx_fix = matrix(0, 1, data.list$dx-1), covy_fix = matrix(0, 1, data.list$dw)))
  } else{
    data <- c(data.list, covars)
  }
  model <- list(type = type, name = paste0(name,"-",type))

  data.object <- make_regmat_unbal(data, model)
  # data_object$model <- build_covnames(data_object)
  return(data.object)
}
