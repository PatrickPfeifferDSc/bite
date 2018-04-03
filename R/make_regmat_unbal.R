#' Transscribes Basic Object to Internal Matrix
#'
#' [Internal Function] Takes the basic objects after read-in of base.mat
#' and panel.mat and restructures them to an internal data object to proceed.
#' Constructs a proper model matrix.
#'
#' @param data The collected information and data in one internal object
#' @param model Information on the model type, etc.
#'
#' @return Internal data object, containing all information about parameters, model, etc.
make_regmat_unbal <- function (data, model){
  # make_regmat_unbal method (for method read_data_unbal)
  # This method constructs a proper Model Matrix of the subject features

  if(!exists("covy_common", where=data)){
    data$covy_common <- matrix(0, ncol=1, nrow=data$dw)
  }
  model$Wx <- data$Wx
  model$dx <- data$dx
  dys <- model$dys <- sum(!data$covy_common)+1

  W0 <- cbind(rep(1,data$Tn), data$Wi[,!data$covy_common])
  Wcom <- data$Wi[,as.logical(data$covy_common)]
  if (model$type == "SRF" || model$type == "SRI" || model$type == "SF"){
    W <- matrix(0, ncol=dys*2, nrow=data$Tn)
    W[, 1:dys] <- W0
    W[data$indy1, (dys+1):(2*dys)] <- cbind(matrix(1, nrow=data$ny1, ncol=1), data$Wi[data$indy1, !data$covy_common])
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
