#' Conditional Moments Computation
#'
#' Internal. \code{condmoms_util} returns an object consisting of conditional moments
#'
#' Computes conditional mean and variance in the utility x* model.
#' Construct array resyt (Residuals of y_t):
#' n x Tmax with zeros for not observed values
#'
#' @param mu_xst mean vector of latent utility values
#' @param resy residual vector on y
#' @param start_x matrix of Tmax x trtmax holding information at which indicator
#' value observations with panel time t and treatment j start
#' @param start_y analoguous to start_x for outcome values
#' @param nx T x trtmax matrix containing number of observations for each combination of t and j
#' @param sgma variance matrix
#' @param rho covariance matrix
#' @param lambda (optional) reserved parameter for further programming
#'
#' @return Returns a list of moments of centered x:
#'         \enumerate{
#'         \item conditional mean
#'         \item conditional variance
#' }

condmoms_util <- function(mu_xst,resy,start_x,nx,start_y,sgma,rho,lambda=NULL){
  n <- length(mu_xst)
  Tmax <- dim(start_x)[1]
  Tmin <- which(start_x[,1]==1)
  mxc <- var_xc <- matrix(0, n, 1)         # column vectors
  if (is.null(lambda)){
    for (j in 1:2){
      rhoj <- rho[,j]
      cj <- rhoj/sgma[,j]
      for (t in Tmin:Tmax){
        condvar <- 1-sum(rhoj[1:t]^2)
        nx_tj <- nx[t,j]
        ind <- start_x[t,j] + (0:(nx_tj-1))    # indices in the x-vector
        indy <- start_y[t,j] + (0:((t*nx_tj)-1))
        resy_tj <- matrix(as.vector(resy[indy,1]),nrow=t, ncol=nx_tj)
        mxc[ind] <- mu_xst[ind] + matrix(cj[1:t], nrow=1, ncol=t) %*% resy_tj
        var_xc[ind] <- condvar
      }
    }
  }  # there might be an additional condition, but till now it does not converge well: lambda given to condmoms_util()
  return(list(mxc=mxc,var_xc=var_xc))
}

