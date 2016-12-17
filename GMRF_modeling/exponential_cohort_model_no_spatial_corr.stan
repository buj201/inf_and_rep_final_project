data {
  //hyperparameter for variance in Y
  real a;
  real b;

  //Indices
  int K; //the number of tracts
  int G; //the number of grades
  int T_train; //the number of training years

  //Actual data
  matrix[G,T_train] Y[K]; // num of students in each grade, each year, each tract

}
parameters {
  vector[K] rate; //Growth rate for tract[k]
  real<lower=0> v_sq; //variance for change in Ys
}
model {
  v_sq ~ inv_gamma(a,b); //Hyperprior

  //Flat prior covers all possible posterior rate space
  //I.e. double in size every year to lose all students every year.
  rate ~ normal(0,0.01);

  //Note we could learn rates from g in 2:G and t in 2:T_train,
  //But this approach lets us get a range of errors:
  for (k in 1:K) {
    for (g in 2:G){
      for (t in g:(g+2)){
        Y[k,g,t] ~ normal(Y[k,g-1, t-1]*(1 + rate[k]), v_sq);
      }
    }
  }
}
