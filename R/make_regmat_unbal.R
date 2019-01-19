#' Transscribes Basic Object to Internal Matrix
#'
#' [Internal Function] Takes the basic objects after read-in of base.mat
#' and panel.mat and restructures them to an internal data object to proceed.
#' Constructs a proper model matrix, which encapsulates the outcome model
#' covariates for each panel time
#'
#' @param data The collected information and data in one internal object
#' @param model Information on the model type, etc.
#'
#' @return A modelling object with W, the regressor matrix of dimension
#' n*t+1 x p
#'

make_regmat_unbal <- function (data, model){
  # This method constructs model matrices for given base and panel covariables

  if(!exists("covy_common")){
    data$covy_common <- rep(0, data$dw)
    print("common covariates for y not found -> fixed at 0")
  }
  model$Wx <- data$Wx
  model$dx <- data$dx
  dys <- model$dys <- sum(!data$covy_common) + 1

  W0 <- cbind(rep(1,data$Tn), data$Wi[,!data$covy_common])
  Wcom <- data$Wi[ ,as.logical(data$covy_common)]
  if (model$type == "SRF" || model$type == "SRI" || model$type == "SF"){
    W <- matrix(0, data$Tn, dys*2)
    W[ ,1:dys] <- W0
    W[data$indy1, (dys+1):(2*dys)] <- cbind(rep(1, data$ny1), data$Wi[data$indy1, !data$covy_common])
    model$W <- cbind(W, Wcom)
  } else{
    if (model$type=="SFimp"){
      model$W <- cbind(W0, matrix(0, nrow=data$Tn, ncol=model$dys), Wcom, W0,W0,Wcom)
    } else{ stop("invalid model type has been specified")}
  }
  model$dy <- dim(model$W)[2]
  data_object <- list(data=data, model=model)
  return(data_object)
}
