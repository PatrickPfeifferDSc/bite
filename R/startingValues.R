#' Computes starting values for parameter search
#'
#' Internal. Computes starting values for alpha and beta parameters.
#' Alpha parameters are estimated by a GLM fit, beta parameters by a linear model.
#'
#' @param model Contains model parameters
#' @param x contains data
#' @param y contains outcome (target) information
#' @param Tn panel times dimension
#'
#' @return Returns a list of
#'         \enumerate{
#'         \item starting values for \eqn{\alpha_\nu}
#'         \item starting values for \eqn{\beta}
#'         \item residual variance
#' }
#'

startingValues <- function(model, x, y, Tn){
  alphav <- glm(x ~ ., data = as.data.frame(model$Wx[, 2:model$dx]), family = "binomial"(link="probit"))
  beta0 <- solve((t(model$W) %*% model$W)) %*% (t(model$W) %*% y) # should give INV(A)*B another method would be left matrix division: A\B
  res_var <- sum((y - (model$W %*% beta0)) ^ 2) / (Tn - dim(model$W)[2])
  return(list(alphav = alphav$coefficients, beta0 = beta0, res_var = res_var))
}

