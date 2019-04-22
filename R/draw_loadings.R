#' Helpfunction for Model SRF latent Factor Switching Regression
#'
#' Internal.

draw_loadings <- function(resy, Xf, resx, start_x, nx, start_y, sgma, rho, prior){

  Tmax <- dim(start_x)[1]
  Tmin <- which(start_x[,1]==1)

  invLn <- diag(prior.invL0)
  lnh <- 0

  # Compute conditional mean and conditional variance
  # construct array resyt: n x Tmax with zeros for not observed values
  for (j in 1:2){
    sgmaj <- sgma[,j]
    rhoj <- rho[,j]

    for (t in Tmin:Tmax){
      sgma_tj <- sgmaj[1:t]
      omega_tj <- rhoj[1:t] * sgma_tj

      inv_sig_condtj <- solve(diag(sgma_tj^2)-(omega_tj%*%t(omega_tj)))

      nx_tj <- nx[t,j]
      indx <- start_x[t,j] + (0:(nx_tj -1)) #indices in the x-vector
      indy <- start_y[t,j] + (0:((t*nx_tj)-1))
      res_tj <- matrix(as.vector(resy[indy]), t, nx_tj) - (omega_tj %*% t(resx[indx]))

      for (i in 1:nx_tj){
        Xfi <- Xf[indy[(i-1) * t + (1:t)], ]
        invLn <- invLN + ( t(Xfi) %*% inv_sig_condtj %*% Xfi )
        lnh <- lnh + t(Xfi) %*% inv_sig_condtj %*% res_tj[,i]
      }
    }
  }
  Ln <- solve(invLn)
  ln <- Ln %*%lnh
  lambda <- (chol(Ln) %*% matrix(rnorm(2*Tmax), 2*Tmax, 1)) + ln
  return(lambda)
}
