/*
  *Simple model for single grade
*/
functions {
  vector CAR_mean(real rho, vector beta, int K, matrix W, vector sum_Ws) {
    vector[K] mean_vec;
    mean_vec = (rho * W * beta)./(rho * sum_Ws + (1 - rho));
    return mean_vec;
  }
  vector CAR_var(real rho, real tau, int K, vector sum_Ws) {
    vector[K] var_vec;
    var_vec = (tau)./(rho * sum_Ws + (1 - rho));
    return var_vec;
  }
  vector center(vector x) {
    return (x - mean(x));
  }
}
data {
  real a; //hyperparameter
  real b; //hyperparameter
  real mu_alpha; //hyperparameter
  real sigma_sq_alpha; //hyperparameter
  int K; //the number of tracts
  int T_train; //the number of training years
  matrix[K,T_train] Y; // num of students in each tract in each year
  matrix[K,K] W; //adjacency matrix
  vector[K] sum_Ws; //sum of W's to avoid recomputation
  row_vector[T_train] Time; //the linear temporal trend covariate
}
parameters {
  real<lower=0> beta; //common intercept
  real alpha; //common slope

  //Impose sum to 0 on phis and deltas
  //See https://groups.google.com/forum/#!topic/stan-users/Ur7YPqUoUYM
  vector[K] phi_raw; //uncentered spatial intercept component
  vector[K] delta_raw; //uncentered spatial slope component

  real<lower=0> tau_sq_int; //variance for ints phi
  real<lower=0> tau_sq_slo; //variance for slopes delta
  real<lower=0, upper=1> rho_int; //correlation control parameter for ints phi
  real<lower=0, upper=1> rho_slo; //correlation control parameter for slopes delta
  real<lower=0> v_sq; //variance for Ys
}
transformed parameters {
  matrix[K,T_train] psi; //linear model mean
  vector[K] phi; //scaled and centered phi
  vector[K] delta; //scaled and centered phi

  phi = center(phi_raw);
  delta = center(delta_raw);
  psi = rep_matrix((beta + phi),T_train) + (alpha + delta)*Time;
}
model {
  tau_sq_int ~ inv_gamma(a,b); //Hyperprior
  tau_sq_slo ~ inv_gamma(a,b); //Hyperprior
  v_sq ~ inv_gamma(a,b); //Hyperprior
  rho_int ~ uniform(0,1); //Hyperprior
  rho_slo ~ uniform(0,1); //Hyperprior
  alpha ~ normal(mu_alpha,sigma_sq_alpha); //Hyperprior

  //CAR (Gaussin Markov Random Field) Priors
  //CAR_mean(rho_0,beta_0, K, W, sum_Ws), CAR_var(rho_0,tau_sq_0, K, sum_Ws)
  phi ~ normal(CAR_mean(rho_int,phi, K, W, sum_Ws), CAR_var(rho_int,tau_sq_int, K, sum_Ws));
  delta ~ normal(CAR_mean(rho_slo,delta, K, W, sum_Ws), CAR_var(rho_slo,tau_sq_slo, K, sum_Ws));
  for (k in 1:K) {
    Y[k] ~ normal(psi[k],v_sq); //Final emission distribution
  }
}
