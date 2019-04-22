#'
#' Draw factors for SRF/SRI
#'
#' Internal.

draw_factor_switchreg <- function(resy, resx, start_x, nx , start_y, sgma, rho, lambda, D){

  n <- length(resx)
  f <- matrix(0, n, 1)
  Tmax <- dim(start_x)[1]
  Tmin <- which(start_x[,1]==1)

  # Compute conditional mean and conditional variance
  # construct array resyt: n x Tmax with zeros for not observed values
  for (j in 1:2){
    sgmaj <- sgma[,j]
    rhoj <- rho[,j]
    lambdaj <- lambda[,j]
    Dinvj <- 1/(D[j])
    for (t in Tmin:Tmax){
      sgma_tj <- sgmaj[1:t]
      omega_tj <- cbind(rhoj[1:t] * sgma_tj)
      lambda_tj <- lambdaj[1:t]
      inv_sig_condtj <- solve(diag(sgma_tj^2) - omega_tj %*% t(omega_tj))

      Fh <- t(lambda_tj) %*% inv_sig_condtj
      F_tj <- solve(Dinvj + (Fh %*% lambda_tj))

      nx_tj <- nx[t,j]
      indx <- start_x[t,j]+(0:(nx_tj-1)) # indices in the x-vector
      indy <- start_y[t,j] + (0:(t*nx_tj-1))
      res_tj <- matrix(as.vector(resy[indy]), t, nx_tj) - omega_tj %*% t(resx[indx])
      f_tj = F_tj %*% Fh %*% res_tj
      f[indx] <- t(f_tj) + (sqrt(as.numeric(F_tj)) * matrix(rnorm(n=nx_tj), nx_tj, 1))
    }
  }
  return(f)
}
