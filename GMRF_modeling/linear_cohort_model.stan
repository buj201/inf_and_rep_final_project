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
}
data {
  //hyperparameter for GMRF hyperpriors
  real a;
  real b;

  //Indices
  int K; //the number of tracts
  int G; //the number of grades
  int T_train; //the number of training years

  //Actual data
  matrix[G,T_train] Y[K]; // num of students in each grade, each year, each tract

  matrix[K,K] W; //adjacency matrix
  vector[K] sum_Ws; //sum of W's to avoid recomputation
}
parameters {
  vector[K] rate; //Growth rate for tract[k]

  real<lower=0> tau_sq; //variance for rates GMRF prior
  real<lower=0, upper=1> rho; //spatial correlation control parameter
  real<lower=0> v_sq; //variance for change in Ys
}
model {
  tau_sq ~ inv_gamma(a,b); //Hyperprior
  v_sq ~ inv_gamma(a,b); //Hyperprior
  rho ~ uniform(0,1); //Hyperprior

  //CAR (Gaussin Markov Random Field) Prior

  rate ~ normal(CAR_mean(rho,rate, K, W, sum_Ws), CAR_var(rho,tau_sq, K, sum_Ws));
  //Note we could learn rates from g in 2:G and t in 2:T_train,
  //But this approach lets us get a range of errors:
  for (k in 1:K) {
    for (g in 2:G){
      for (t in g:(g+2)){
        Y[k,g,t] ~ normal(Y[k,g-1, t-1] + rate[k], v_sq);
      }
    }
  }
}
