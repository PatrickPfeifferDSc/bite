#' Read Data file and transform
#'
#' \code{read_data_unbal} takes object of base matrix and full panel matrix
#' to read, structure and organise the data into a useable data.frame R object
#'
#' @param infile a list of base matrix and panel matrix
#' @param type contains a string of possible model types ("SF")
#' @param covars contains additional information about covariate names (optional)
#' @param ... control parameter passed down in funcional hierarchy
#'
#' @return A data object for further computation.
#'
#' @export
#'

read_data_unbal <- function(infile, type, covars = list(NULL), ...){
  # read_data_unbal method

  # function takes infile(=list of base matrix and panel matrix), type(=model type), covars(=list of covariates which are to be initialised)

  # at the moment this function takes an R list as argument infile, before a .mat file is read to a list type object
  #type needs to be one of ("SF", "SFimp", "SRI", "SRF")

  if (!(any(names(infile)== "base.mat") && any(names(infile)=="panel.mat"))) stop("must include list objects called 'base.mat' and 'panel.mat', representing the base matrix and panel matrix of the panel")
  # checking whether the object is structured correctly, at the moment must include 3 objects: base.mat, panel.mat, covars
  # read and check base matrix
  ldata <- list()
  ldata$name <- deparse(substitute(infile))
  n <- dim(infile$base.mat)[1];  k <- dim(infile$base.mat)[2]
  id <- infile$base.mat[,1]
  Ti <- infile$base.mat[,2]
  x <- infile$base.mat[,3]
  # ensure that treatment is binary 0/1
  if ((max(x) - min(x))!=1) stop("Treatment specifications might be wrong, must be binary coded, 0/1 optimally).")
  if (min(x) != 0) x <- x - min(x)
  # displace the distribution of treatment values to be (0/1) binary
  dwx <- k-2  # potential number of regressors used further: intercept + column 4 - column k from base.mat

  # read and check panel matrix (panel matrix includes 1 panel time per observation)
  kp <- dim(infile$panel.mat)[2]
  dw <- kp - 3
  Tn <- sum(Ti)   #sum over panel times
  ind.up <- 0

  for (i in 1:n){
    Tiy <- Ti[i]
    ind.low <- ind.up + 1
    ind.up <- ind.up + Tiy
    id.panel <- infile$panel.mat[ind.low:ind.up, 1]
    if (sum(id.panel==id[i])<Tiy) stop("subjects are not in the same order as in base file or obs. missing")
  }
  # end

  tpanel <- infile$panel.mat[ind.low:ind.up, 2]
  if (sum(t(tpanel)==1:Tiy)<Tiy) stop("subjects are not in the same order as in the base file or obs. missing")
  # end

  # create the finally ordered data set

  ldata$x <- rep(0, 5)
  ldata$Ti <- rep(0, n)
  ldata$indy <- ldata$Tvec <- ldata$y <- rep(0, Tn)
  ldata$Wi <- matrix(0, ncol=dw, nrow=Tn)
  ldata$Wx <- matrix(NA, ncol=k-2, nrow=n)
  ldata$n <- n
  ldata$ntrt <- 2
  ldata$Tn <- Tn
  ldata$dw <- dw
  ldata$dx <- dwx
  Tmin <- min(Ti)
  Tmax <- max(Ti)

  # create starting indices for treatment group with t panel observations
  ldata$start_x <- ldata$start_y <- matrix(0, ncol=ldata$ntrt, nrow=Tmax)

  Tend <- cumsum(Ti)
  Tstart <- c(1, Tend[1:ldata$n-1] +1)
  # Start indices for all vectors
  ix <- 1
  iy <- 0
  ldata$nx <- matrix(0, ncol=2, nrow=Tmax)
  # finds indices of rows, which have the same paneltimes and treatment x, this is used for sorting the observations

  for (j in 0:1){
    for (t in Tmin:Tmax){
      indjt <- which(x == j & Ti == t)
      hn <- length(indjt)
      ldata$nx[t, j+1] <- hn
      ldata$start_x[t, j+1] <- ix
      ldata$start_y[t, j+1] <- iy+1
      it <- 1:t
      for (i in 1:hn){
        indx <- indjt[i]
        ldata$x[ix] <- x[indx]
        ldata$Ti[ix] <- Ti[indx]
        ldata$Wx[ix,] <- c(1, infile$base.mat[indx, 4:k])
        indy <- Tstart[indx]:Tend[indx]
        ldata$indy[iy+it] <- ix
        ldata$Tvec[iy+it] <- infile$panel.mat[indy, 2]
        ldata$y[iy+it] <- infile$panel.mat[indy, 3]
        ldata$Wi[iy+it,] <- infile$panel.mat[indy, 4:kp]
        ix <- ix+1
        iy <- iy+t
      }
    }
  }
  ldata$nx0 <- sum(1-ldata$x)
  ldata$nx1 <- sum(ldata$x)
  ldata$indy0 <- as.logical(1-ldata$x[ldata$indy])
  ldata$indy1 <- as.logical(ldata$x[ldata$indy])
  ldata$ny0 <- sum(ldata$indy0==1)
  ldata$ny1 <- sum(ldata$indy1==1)
  # read covariate names if available
  covarsnames <- c("covx", "covx_fix", "covy", "covy_fix", "covy_common", "covy_nb", "covy_np")
  if (is.null(covars) || !(all(covarsnames %in% names(covars))) ){               # checks correctness of covariates parameters
    warning("No covars object found, fixed values will be set to 0.")
    data <- c(ldata, list(covx_fix=matrix(0, ncol=(ldata$dx-1), nrow=1), covy_fix= matrix(0, ncol=(ldata$dx), nrow=1)))
  } else{
    data <- c(ldata, covars)
  }
  model <- list(type=type, name=paste(data$name,type))

  data_object <- make_regmat_unbal(data,model)
  # data_object$model <- build_covnames(data_object)
  return(data_object)
}
