/*
  *Linear model for all grades
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

  //Indices
  int K; //the number of tracts
  int G; //the number of grades
  int T_train; //the number of training years

  //Actual data
  matrix[G,T_train] Y[K]; // num of students in each grade, each year, each tract

  matrix[K,K] W; //adjacency matrix
  vector[K] sum_Ws; //sum of W's to avoid recomputation
  matrix[G,T_train] Time; //Time matrix (with constant rows)
  matrix[G,T_train] Grade; //Grade matrix (with constant columns)
}
parameters {
  real<lower=0> alpha_0; //common intercept across all tracts
  real alpha_g; //common grade coefficient across all tracts
  real alpha_t; //common time coefficient across all tracts

  //Impose sum to 0 on phis and deltas
  //See https://groups.google.com/forum/#!topic/stan-users/Ur7YPqUoUYM

  vector[K] beta_0_raw; //uncentered intercept component
  vector[K] beta_g_raw; //uncentered grade coefficient component
  vector[K] beta_t_raw; //uncentered time coefficient component

  real<lower=0> tau_sq_0; //variance for intercepts
  real<lower=0> tau_sq_g; //variance for grades coefficient
  real<lower=0> tau_sq_t; //variance for time coefficient

  real<lower=0, upper=1> rho_0; //corr control param for intercepts
  real<lower=0, upper=1> rho_g; //corr control param for grade
  real<lower=0, upper=1> rho_t; //corr control param for time
  real<lower=0> v_sq; //variance for Ys
}
transformed parameters {
  vector[K] beta_0; //scaled and centered beta_0_raw
  vector[K] beta_g; //scaled and centered beta_g_raw
  vector[K] beta_t; //scaled and centered beta_t_raw
  matrix[G,T_train] mu[K]; //Mean for emission distribution

  beta_0 = center(beta_0_raw);
  beta_g = center(beta_g_raw);
  beta_t = center(beta_t_raw);

  for (k in 1:K){
    //Remember .* is elementwise operation
    mu[k] = rep_matrix((alpha_0 + beta_0[k]),G , T_train) + (alpha_g + beta_g[k]) * Grade + (alpha_t + beta_t[k]) * Time;
  }
}
model {
  v_sq ~ inv_gamma(a,b); //Hyperprior

  tau_sq_0 ~ inv_gamma(a,b); //Hyperprior
  tau_sq_g ~ inv_gamma(a,b); //Hyperprior
  tau_sq_t ~ inv_gamma(a,b); //Hyperprior

  rho_0 ~ uniform(0,1); //Hyperprior
  rho_g ~ uniform(0,1); //Hyperprior
  rho_t ~ uniform(0,1); //Hyperprior

  //CAR priors
  beta_0 ~ normal(CAR_mean(rho_0,beta_0, K, W, sum_Ws), CAR_var(rho_0,tau_sq_0, K, sum_Ws));
  beta_g ~ normal(CAR_mean(rho_g,beta_g, K, W, sum_Ws), CAR_var(rho_g,tau_sq_g, K, sum_Ws));
  beta_t ~ normal(CAR_mean(rho_t,beta_t, K, W, sum_Ws), CAR_var(rho_t,tau_sq_t, K, sum_Ws));

  for (k in 1:K) {
    for (g in 1:G){
      Y[k,g] ~ normal(mu[k,g],v_sq); //Final emission distribution
    }
  }
}
