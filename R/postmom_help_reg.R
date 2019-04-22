#' postmom_help_reg method
#'
#' Internal. Computes help quantities for posterior moments of regression coefficients.
#'
#' @param x data
#' @param y outcome
#' @param Wx design matrix for selection model
#' @param Wy design matrix for outcome data model
#' @param start_x matrix containing indices where certain observations
#'          begin in the vector sorted by treatment and panel times
#' @param nx dimension of x
#' @param start_y matrix containing indices where certain outcomes
#'          begin in the vector sorted by treatment and panel times
#' @param sgma variance parameter matrix
#' @param rho correlation parameter matrix
#' @param lambda matrix of factor loadings
#' @param D initial inclusion for treatment j
#'
#' @return List of Matrices used for computation of posterior moments.

postmom_help_reg <- function(x, y, Wx, Wy, start_x, nx, start_y, sgma, rho, lambda, D){

  Tmax <- dim(start_x)[1]
  Tmin <- which(start_x[,1]==1)

  WSW <- 0
  WSy <- 0
  ySy <- 0

  dx <- dim(Wx)[2]
  dy <- dim(Wy)[2]

  for (j in 1:2){
    sgmaj <- sgma[,j]
    rhoj <- rho[,j]
    lambdaj <- lambda[,j]
    Dj <- D[j]
    for (t in Tmin:Tmax){
      sgma_tj <- sgmaj[1:t]
      omega_tj <- rhoj[1:t]*sgma_tj
      lambda_tj <- lambdaj[1:t]
      Lambdaj <- cbind(c(1, omega_tj), rbind(omega_tj, (diag(sgma_tj^2)+(t(t(lambda_tj))%*%lambda_tj)*Dj)))
      ILambdaj <- solve(Lambdaj)
      nx_tj <- nx[t,j]
      indx <- start_x[t,j] + (0:(nx_tj-1))      # indices in the x-vector
      indy <- start_y[t,j] + (0:(t*nx_tj-1))
      yh_tj <- rbind(x[indx], matrix(y[indy], t , nx_tj))
      Wi <- matrix(0, t+1, dx+dy)
      for (i in 1:nx_tj){
        yhi <- cbind(yh_tj[,i])
        Ly <- ILambdaj %*% yhi
        ySy <- ySy + (t(yhi) %*% Ly)
        Wi[1,1:dx] <- Wx[indx[i], ]
        Wi[2:(t+1), dx + (1:dy)] <- Wy[indy[(i-1)*t + (1:t)], ]
        WiV <- t(Wi) %*% ILambdaj
        WSW <- WSW + WiV %*% Wi
        WSy <- WSy + WiV %*% yhi
      }
    }
  }
  return(list(WSW=WSW, ySy=ySy, WSy=WSy))
}

