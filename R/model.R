model_code <- 'functions {
  real p_matrix(real t, real s, real gamma, real lambda) {
    if (s == 0) {
      return ( (lambda - lambda*exp(-t*(gamma+lambda)))/(gamma+lambda) );
    } else {
      return ( (lambda + gamma*exp(-t*(gamma+lambda)))/(gamma+lambda) );
    }
  }
}

data {
  int num_individuals;
  int num_obs[num_individuals];
  real times[num_individuals, max(num_obs)];
  int S[num_individuals, max(num_obs)];
}

parameters {
  real<lower=1e-10, upper=1> lambda;
  real<lower=1e-10, upper=1> gamma;
}

transformed parameters {
  real p[num_individuals, max(num_obs)];
  for (j in 1:num_individuals) {
    p[j][1] = 0;
    for (i in 2:num_obs[j]) {
      p[j][i] = p_matrix(times[j][i]-times[j][i-1], S[j][i-1], gamma, lambda);
    }
  }
}

model {
  lambda ~ uniform(1e-10, 1);
  gamma ~ uniform(1e-10, 1);
  for (j in 1:num_individuals)
  {
    for (t in 2:num_obs[j]) {
      if (S[j][t] != -1)
        S[j][t] ~ bernoulli(p[j][t]);
    }
  }
  
  
}'

tmi_model <- stan_model(model_code = model_code)
