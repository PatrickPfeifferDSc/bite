#' Computes starting values for parameter search
#'
#' Computes starting values for alpha and beta parameters. Alpha parameters
#' are estimated by a GLM fit, beta parameters by a linear model.
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
#' @import stats

start_vals <- function(model, x, y, Tn){
  alphav <- glm(x ~ ., data= as.data.frame(model$Wx[,2:model$dx]), family = "binomial"(link="probit"))
  beta0 <- solve((t(model$W) %*% model$W)) %*% (t(model$W)%*% y) # should give INV(A)*B another method would be left matrix division: A\B
  res_var <- sum((y-(model$W%*%beta0))^2)/(Tn-dim(model$W)[2])
  iheck <- 0

  # if iheck==1,
  # y=data.y;
  # x=data.x;
  # % Heckman starting values
  # Tmax=max(data.Tvec);
  #
  # mux=model.Wx*alphav;
  # hx=zeros(data.n,1);
  #
  # ind0=(x==0);
  # m0=mux(ind0);
  # hx(ind0)=normpdf(m0)./normcdf(-m0);
  #
  # ind1=(x==1);
  # m1=mux(ind1);
  # hx(ind1)=normpdf(m1)./normcdf(m1);
  #
  # for t=1:Tmax;
  # indt(:,t)=data.Tvec==t;
  # end

  # hy=hx(data.indy);
  # H=repmat(hy,1,Tmax).*indt;
  #
  # ny0=sum(data.indy0);
  # ny1=data.Tn-ny0;
  #
  # H01=sparse(zeros(ny0,Tmax));
  # H10=sparse(zeros(ny1,Tmax));
  #
  # Wh=[H(1:ny0,:), H01
  #     H10,H(ny0+1:data.Tn,:)];
  #
  # W=[model.W,Wh];
  # beta0=(W'*W)\(W'*y)

  # Heckman coefs may have not been as good for the given data -> revisit for general data purpose ?

  return(list(alphav = alphav$coefficients, beta0 = beta0, res_var = res_var))
}

