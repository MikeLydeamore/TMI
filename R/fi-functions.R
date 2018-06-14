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

.PMatrix <- function(delta, gamma, lambda)
{
  matrix(c(gamma+lambda*exp(-delta*(gamma+lambda)), lambda-lambda*exp(-delta*(gamma+lambda)), gamma-gamma*exp(-delta*(gamma+lambda)), lambda+gamma*exp(-delta*(gamma+lambda)))/(gamma+lambda), byrow = T, nrow=2)
}

.phi_s <- function(delta, gamma, lambda)
{
  (lambda - lambda*exp(-delta*(gamma+lambda)))/(gamma+lambda)
}

.phi_i <- function(delta, gamma, lambda)
{
  (lambda + gamma*exp(-delta*(gamma+lambda)))/(gamma+lambda)
}


#' Calculate Fisher Information Matrix in the two-state model
#' 
#' Calculates the determinant of the Fisher Information matrix in a two-state model
#' for a given set of parameters and time between observations.
#' 
#' @param gamma Rate of recovery
#' @param lambda Force of infection
#' @param delta Vector of times between observations
#' @param i0 Initial probability of being infected (prevalence)
#' 
#' @return The determinant of the Fisher Information
calculateFisherInformation <- function(gamma, lambda, delta, i0=0.316)
{
  
  x_0 <- c(1-i0, i0)
  
  top_left <- 0
  top_right <- 0
  bottom_left <- 0
  bottom_right <- 0
  
  for (i in 1:length(delta))
  {
    phi_s <- .phi_s(delta[i], gamma, lambda)
    phi_i <- .phi_i(delta[i], gamma, lambda)
    
    if (i == 1)
      P_iminusone <- .PMatrix(delta[i], gamma, lambda)
    else
      P_iminusone <- P_iminusone %*% .PMatrix(delta[i], gamma, lambda)
      
    top_left <- top_left + sum((x_0 %*% P_iminusone)[,1])*( (1-phi_s)*.one_minus_phi_s_lambda_lambda(gamma,lambda,delta[i]) + phi_s*.phi_s_lambda_lambda(gamma,lambda,delta[i]) ) +
      sum((x_0%*% P_iminusone)[,2])*( (1-phi_i)*.one_minus_phi_i_lambda_lambda(gamma,lambda,delta[i]) + phi_i*.phi_i_lambda_lambda(gamma,lambda,delta[i]) )
    
    top_right <- top_right + sum((x_0 %*% P_iminusone)[,1])*( (1-phi_s)*.one_minus_phi_s_gamma_lambda(gamma,lambda,delta[i]) + phi_s*.phi_s_lambda_gamma(gamma,lambda,delta[i]) ) +
      sum((x_0%*% P_iminusone)[,2])*( (1-phi_i)*.one_minus_phi_i_lambda_gamma(gamma,lambda,delta[i]) + phi_i*.phi_i_lambda_gamma(gamma,lambda,delta[i]) )

    bottom_left <- bottom_left + sum((x_0 %*% P_iminusone)[,1])*( (1-phi_s)*.one_minus_phi_s_gamma_lambda(gamma,lambda,delta[i]) + phi_s*.phi_s_lambda_gamma(gamma,lambda,delta[i]) ) +
      sum((x_0%*% P_iminusone)[,2])*( (1-phi_i)*.one_minus_phi_i_lambda_gamma(gamma,lambda,delta[i]) + phi_i*.phi_i_lambda_gamma(gamma,lambda,delta[i]) )

    bottom_right <- bottom_right + sum((x_0 %*% P_iminusone)[,1])*( (1-phi_s)*.one_minus_phi_s_gamma_gamma(gamma,lambda,delta[i]) + phi_s*.phi_s_gamma_gamma(gamma,lambda,delta[i]) ) +
      sum((x_0%*% P_iminusone)[,2])*( (1-phi_i)*.one_minus_phi_i_gamma_gamma(gamma,lambda,delta[i]) + phi_i*.phi_i_gamma_gamma(gamma,lambda,delta[i]) )
  }
  
  return ((top_left*bottom_right)-(top_right*bottom_left))
}
