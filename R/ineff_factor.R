#' Inefficiency Factor
#'
#' \code{ineff_factor} calculates inefficiency facors for certain time series
#' with autocorrelation after Geyer (1992).
#'
#' @param x a time series of length k
#'
#' @return empirical inefficiency factor

ineff_factor <- function(x){

  # x time series of length (k)
  # inefficiency factor : initial monotonoe sequence Geyer (1992)

  nlag <- 1000  # maximum_lag
  if(length(x)<=1001) nlag <- length(x)-2
  n <- length(x)
  mx <- mean(x)
  vx <- var(x)
  c <- rep(0, nlag)

  # empirical autocovariances
  for (i in 1:nlag){
    x1 <- x[(i+1):n] - mx
    x2 <- x[1:(n-i)] - mx
    c[i] <- (x1%*%x2)/(n*vx)
  }
  index <- seq(2,nlag-1, by=2)
  gam <- c[index] + c[index + 1]
  h1 <- which(gam < 0)

  if (length(h1) < 1){
    m1 <- (nlag-1)/2
  }else{
    m1 <- min(h1)-1 }

  if(m1 > 0){
    dGam <- diff(gam)
    h2 <- which(dGam > 0)
    if (length(h2) <1){m <- m1
    }else{ m <- min(c(m1, min(h2)-1)) }

  }else{
    m <- m1
  }
  tau <- 1 + 2*sum(c[1:(2*m+1)])
  return(list(tau=tau, m=m))
}

