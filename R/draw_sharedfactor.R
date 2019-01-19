#' Draw Shared Factor
#'
#' [Internal Function] \code{draw_sharedfactor} draws the latent factors \eqn{f_i} for
#' SF model.
#'
#' @param epsx residual vector on x
#' @param lambdax lambda parameter of selection model x*
#' @param epsy residuals on y
#' @param sgma2 variance parameters
#' @param lambda lambda parameters on outcome for each panel time and treatment
#' @param start_x indicator matrix
#' @param nx number of observations
#' @param start_y indicator matrix on y
#' @param Tn number of observations times panel t
#'
#'@return f ... shared latent factor

drawSharedFactor <- function(epsx,lambdax,epsy,sgma2,lambda,start_x,nx,start_y,Tn, fix.f=FALSE){

  n <- length(epsx)

  Tmax <- dim(start_x)[1]
  Tmin <- which(start_x[,1]==1)

  lambdah <- matrix(as.vector(lambda), Tmax, 2)

  f <- matrix(0, n, 1)

  # is the first if statement to check for epsilon y being a vector of all observed outcomes in panel (twice for treatment=0/1)
  # in which format should the data come here?

  if (length(epsy)==(2*Tn)){
    for (t in Tmin:Tmax){
      # Compute Posterior variance matrix of latent factors
      X <- cbind(c(lambdax, lambdah[1:t,1], lambdah[1:t,2]))
      Sinv <- diag(c(1,1/sgma2[1:t,1],1/sgma2[1:t,2]))
      Xs <- t(X)%*%Sinv
      Ft <- solve(1 + (Xs%*%X))

      # Build regressor
      indx <- c( (start_x[t,1] + (0:(nx[t,1] - 1))), (start_x[t,2] + (0:(nx[t,2]-1))) )
      indy <- cbind(matrix(as.vector(start_y[t,1] + (0:(t*nx[t,1]-1))), t, nx[t,1]) , matrix(as.vector(start_y[t,2]+(0:(t*nx[t,2]-1))), t, nx[t,2]) )
      y <- rbind(t(epsx[indx]), epsy[indy], epsy[indy + Tn])
      ft <- Ft%*%(Xs%*%y)
      nt <- sum(nx[t,])
      f[indx,1] <- (sqrt(Ft) %*% rnorm(nt)) + ft
    }
  }else{
    # if statement for data being just panel times, then going through for all the treatment
    if (length(epsy)==Tn){
      for (j in 1:2){
        for (t in Tmin:Tmax){
          # Compute Posterior variance matrix of latent factors
          X <- c(lambdax,   lambdah[1:t, j]); dim(X) <- c(length(lambdah[1:t, j])+1, 1);
          Sinv <- diag(c(1,1/sgma2[1:t,j]))
          Xs <- t(X) %*% Sinv
          Ft <- solve(1 + Xs%*%X)
          # Build regressor
          indx <- start_x[t,j] + (0:(nx[t,j]-1))    # evaluate correct index for different panel times
          indy <- matrix(as.vector(start_y[t,j] + (0:(t*nx[t,j]-1))),t , nx[t,j])
          y <- epsx[indx]
          for (i in 1:t){
            y <- rbind(y, epsy[indy[i,]])
          }
          ft <- Ft%*%(Xs%*%y)
          nt <- sum(nx[t,j])
          f[indx,1] <- (sqrt(Ft)%*%rnorm(nt)) + ft
        }

      }
    }else{stop('dimension error in draw_shared factor')}
  }
return(f)
}
