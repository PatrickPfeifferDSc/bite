#' Transscribes Basic Object to Internal Matrix
#'
#' Internal. Takes the basic objects after read-in of base.mat
#' and panel.mat and restructures them to an internal data object to proceed.
#' Constructs a proper model matrix, which encapsulates the outcome model
#' covariates for each panel time
#'
#' @param data The collected information and data in one internal object
#' @param model Information on the model type, etc.
#'
#' @return A modelling object with W, the regressor matrix of dimension n*t+1 x p
#'

make_regmat_unbal <- function (data, model){
  if(is.null(data$cov_y_common)){
    data$cov_y_common <- rep(0, data$dw)
    print("No common covariates for y defined, normal, treatment-separated estimation for all covariates")
  }
  model$Wx <- data$Wx
  model$dx <- data$dx
  model$dys <- dys <- sum(!data$cov_y_common) + 1     # dimension of covariates in heterogeneous treatment + intercept
  W0 <- cbind(rep(1,data$Tn), data$Wi[ ,!data$cov_y_common])
  Wcom <- data$Wi[ ,as.logical(data$cov_y_common)]

  if (model$type %in% c("SF", "SRF", "SRI")){
    W <- matrix(0, data$Tn, dys*2)
    W[ ,(1:dys)] <- W0
    W[data$indy1, (dys+1):(2*dys)] <- cbind(rep(1, data$ny1), data$Wi[data$indy1, !data$cov_y_common])
    model$W <- cbind(W, Wcom)
  }else{
    stop("invalid model type has been specified")
  }
  model$dy <- dim(model$W)[2]
  data_object <- list(data=data, model=model)
  return(data_object)
}
