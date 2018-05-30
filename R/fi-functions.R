.one_minus_phi_i_gamma_gamma <- function(gamma, lamb, t)
{
  -(2*gamma*lamb + lamb^2 + (2*gamma*lamb + lamb^2)*exp(2*gamma*t + 2*lamb*t) + ((gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2)*t^2 - 4*gamma*lamb - 2*lamb^2)*exp(gamma*t + lamb*t))/(gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2 + (gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2)*exp(2*gamma*t + 2*lamb*t) - 2*(gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2)*exp(gamma*t + lamb*t))
}

.one_minus_phi_i_lambda_gamma <- function(gamma, lamb, t)
{
  -(((gamma^2 + 2*gamma*lamb + lamb^2)*t^2 + 2)*exp(gamma*t + lamb*t) - exp(2*gamma*t + 2*lamb*t) - 1)/(gamma^2 + 2*gamma*lamb + lamb^2 + (gamma^2 + 2*gamma*lamb + lamb^2)*exp(2*gamma*t + 2*lamb*t) - 2*(gamma^2 + 2*gamma*lamb + lamb^2)*exp(gamma*t + lamb*t))
}

.one_minus_phi_i_lambda_lambda <- function(gamma, lamb, t)
{
  -(((gamma^2 + 2*gamma*lamb + lamb^2)*t^2 + 2)*exp(gamma*t + lamb*t) - exp(2*gamma*t + 2*lamb*t) - 1)/(gamma^2 + 2*gamma*lamb + lamb^2 + (gamma^2 + 2*gamma*lamb + lamb^2)*exp(2*gamma*t + 2*lamb*t) - 2*(gamma^2 + 2*gamma*lamb + lamb^2)*exp(gamma*t + lamb*t))
}

.phi_i_gamma_gamma <- function(gamma, lamb, t)
{
  (lamb^2*exp(2*gamma*t + 2*lamb*t) - 2*gamma*lamb - lamb^2 + ((gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*t^2 + 2*gamma*lamb - 2*(gamma^2*lamb + 2*gamma*lamb^2 + lamb^3)*t)*exp(gamma*t + lamb*t))/(gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2 + (gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4)*exp(2*gamma*t + 2*lamb*t) + 2*(gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*exp(gamma*t + lamb*t))
}

.phi_i_lambda_gamma <- function(gamma, lamb, t)
{
  (lamb^2*exp(2*gamma*t + 2*lamb*t) + gamma^2 + ((gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*t^2 - gamma^2 - lamb^2 + (gamma^3 + gamma^2*lamb - gamma*lamb^2 - lamb^3)*t)*exp(gamma*t + lamb*t))/(gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2 + (gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4)*exp(2*gamma*t + 2*lamb*t) + 2*(gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*exp(gamma*t + lamb*t))
}

.phi_i_lambda_lambda <- function(gamma, lamb, t)
{
  (gamma^2 - (gamma^2 + 2*gamma*lamb)*exp(2*gamma*t + 2*lamb*t) + ((gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*t^2 + 2*gamma*lamb + 2*(gamma^3 + 2*gamma^2*lamb + gamma*lamb^2)*t)*exp(gamma*t + lamb*t))/(gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2 + (gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4)*exp(2*gamma*t + 2*lamb*t) + 2*(gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*exp(gamma*t + lamb*t))
}

.phi_s_lambda_lambda <- function(gamma, lamb, t)
{
  -(gamma^2 + 2*gamma*lamb + (gamma^2 + 2*gamma*lamb)*exp(2*gamma*t + 2*lamb*t) + ((gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4)*t^2 - 2*gamma^2 - 4*gamma*lamb)*exp(gamma*t + lamb*t))/(gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4 + (gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4)*exp(2*gamma*t + 2*lamb*t) - 2*(gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4)*exp(gamma*t + lamb*t))
}

.phi_s_lambda_gamma <- function(gamma, lamb, t)
{
  -(((gamma^2 + 2*gamma*lamb + lamb^2)*t^2 + 2)*exp(gamma*t + lamb*t) - exp(2*gamma*t + 2*lamb*t) - 1)/(gamma^2 + 2*gamma*lamb + lamb^2 + (gamma^2 + 2*gamma*lamb + lamb^2)*exp(2*gamma*t + 2*lamb*t) - 2*(gamma^2 + 2*gamma*lamb + lamb^2)*exp(gamma*t + lamb*t))
}

.phi_s_gamma_gamma <- function(gamma, lamb, t)
{
  -(((gamma^2 + 2*gamma*lamb + lamb^2)*t^2 + 2)*exp(gamma*t + lamb*t) - exp(2*gamma*t + 2*lamb*t) - 1)/(gamma^2 + 2*gamma*lamb + lamb^2 + (gamma^2 + 2*gamma*lamb + lamb^2)*exp(2*gamma*t + 2*lamb*t) - 2*(gamma^2 + 2*gamma*lamb + lamb^2)*exp(gamma*t + lamb*t))
}

.one_minus_phi_s_gamma_gamma <- function(gamma, lamb, t)
{
  (lamb^2 - (2*gamma*lamb + lamb^2)*exp(2*gamma*t + 2*lamb*t) + ((gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*t^2 + 2*gamma*lamb + 2*(gamma^2*lamb + 2*gamma*lamb^2 + lamb^3)*t)*exp(gamma*t + lamb*t))/(gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4 + (gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2)*exp(2*gamma*t + 2*lamb*t) + 2*(gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*exp(gamma*t + lamb*t))
}

.one_minus_phi_s_gamma_lambda <- function(gamma, lamb, t)
{
  (gamma^2*exp(2*gamma*t + 2*lamb*t) + lamb^2 + ((gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*t^2 - gamma^2 - lamb^2 - (gamma^3 + gamma^2*lamb - gamma*lamb^2 - lamb^3)*t)*exp(gamma*t + lamb*t))/(gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4 + (gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2)*exp(2*gamma*t + 2*lamb*t) + 2*(gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*exp(gamma*t + lamb*t))
}

.one_minus_phi_s_lambda_lambda <- function(gamma, lamb, t)
{
  (gamma^2*exp(2*gamma*t + 2*lamb*t) - gamma^2 - 2*gamma*lamb + ((gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*t^2 + 2*gamma*lamb - 2*(gamma^3 + 2*gamma^2*lamb + gamma*lamb^2)*t)*exp(gamma*t + lamb*t))/(gamma^2*lamb^2 + 2*gamma*lamb^3 + lamb^4 + (gamma^4 + 2*gamma^3*lamb + gamma^2*lamb^2)*exp(2*gamma*t + 2*lamb*t) + 2*(gamma^3*lamb + 2*gamma^2*lamb^2 + gamma*lamb^3)*exp(gamma*t + lamb*t))
}

calculateFisherInformation <- function(gamma, lambda, delta, n_obs, i0=0.316)
{
  P <- matrix(c(gamma+lambda*exp(-delta*(gamma+lambda)), lambda-lambda*exp(-delta*(gamma+lambda)), gamma-gamma*exp(-delta*(gamma+lambda)), lambda+gamma*exp(-delta*(gamma+lambda)))/(gamma+lambda), byrow = T, nrow=2)
  
  x_0 <- c(1-i0, i0)
  
  top_left <- 0
  top_right <- 0
  bottom_left <- 0
  bottom_right <- 0
  
  phi_s <- (lambda - lambda*exp(-delta*(gamma+lambda)))/(gamma+lambda)
  phi_i <- (lambda + gamma*exp(-delta*(gamma+lambda)))/(gamma+lambda)
  
  for (i in 1:n_obs)
  {
    top_left <- top_left + sum((x_0 %*% (P%^%(i-1)))[,1])*( (1-phi_s)*.one_minus_phi_s_lambda_lambda(gamma,lambda,delta) + phi_s*.phi_s_lambda_lambda(gamma,lambda,delta) ) +
      sum((x_0%*% (P%^%i))[,2])*( (1-phi_i)*.one_minus_phi_i_lambda_lambda(gamma,lambda,delta) + phi_i*.phi_i_lambda_lambda(gamma,lambda,delta) )
    
    top_right <- top_right + sum((x_0 %*% (P%^%(i-1)))[,1])*( (1-phi_s)*.one_minus_phi_s_gamma_lambda(gamma,lambda,delta) + phi_s*.phi_s_lambda_gamma(gamma,lambda,delta) ) +
      sum((x_0%*% (P%^%i))[,2])*( (1-phi_i)*.one_minus_phi_i_lambda_gamma(gamma,lambda,delta) + phi_i*.phi_i_lambda_gamma(gamma,lambda,delta) )

    bottom_left <- bottom_left + sum((x_0 %*% (P%^%(i-1)))[,1])*( (1-phi_s)*.one_minus_phi_s_gamma_lambda(gamma,lambda,delta) + phi_s*.phi_s_lambda_gamma(gamma,lambda,delta) ) +
      sum((x_0%*% (P%^%i))[,2])*( (1-phi_i)*.one_minus_phi_i_lambda_gamma(gamma,lambda,delta) + phi_i*.phi_i_lambda_gamma(gamma,lambda,delta) )

    bottom_right <- bottom_right + sum((x_0 %*% (P%^%(i-1)))[,1])*( (1-phi_s)*.one_minus_phi_s_gamma_gamma(gamma,lambda,delta) + phi_s*.phi_s_gamma_gamma(gamma,lambda,delta) ) +
      sum((x_0%*% (P%^%i))[,2])*( (1-phi_i)*.one_minus_phi_i_gamma_gamma(gamma,lambda,delta) + phi_i*.phi_i_gamma_gamma(gamma,lambda,delta) )
  }
  
  return ((top_left*bottom_right)-(top_right*bottom_left))
}
