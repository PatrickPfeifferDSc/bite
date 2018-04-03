#' Draw latent utilities
#'
#' [Internal Function] \code{draw_utility} draws latent utilities for the
#' probit model, which represents the selection of treatment. Latent
#' utilities (x*) range from -Inf to +Inf. For x* being positive, the
#' subject selects treatment, subject does not decide for treatment otherwise.
#'
#' @param x vector of subject features
#' @param mxst vector of means of utilities x*
#' @param sxst vector of standard deviations of utilities x*
#'
#' @importFrom truncnorm rtruncnorm
#'
#' @return xst x*, the latent utility for selection model


draw_utility <- function(x, mxst, sxst){
  # draws latent utilities  for the probit model
  # Author: Helga Wagner
  # Last Version : Jan 4, 2012
  #########################################################################
  # Input
  # x                n x 1         treatment effects, numbered 1..ntreat
  # mu_xst=mu(x)     1 x n         fixed effects mean for the utilities xst
  ##########################################################################
  # Output
  # xst              1xn           latent utilities
  n <- length(x)
  xst <- matrix(0,n,1)

  for (j in 0:1){
    idx <- (x==j)
    if (j == 0){
      tau_low <- -Inf
      tau_high <- 0
    }else{
      tau_low <- 0
      tau_high <- Inf
    }
    if(length(sxst)==1){
        xst[idx,1] <- rtruncnorm(1, mean=mxst[idx], sd=sxst, a=tau_low, b= tau_high)
      }else{
        xst[idx,1] <- rtruncnorm(1, mean=mxst[idx], sd=sxst[idx], a=tau_low, b= tau_high)
      }
    }
  return(xst)
}
