source('get_data.R')

#Now learn same model in STAN

library(rstan)
library(coda)
library(matrixStats)

test_flag = FALSE

if (test_flag){
  num_iter = 100 #set to 10000
  warmup = 100
  num_chains = 1#set to 4
  file_string = "./stan_models/lin_model_stan_fit_test.RData"
} else {
  num_iter = 10000
  warmup = num_iter/5
  num_chains = 4 #set to 4
  file_string = "./stan_models/lin_model_stan_fit.RData"
}

options(mc.cores = parallel::detectCores())

a = 1
b = 0.01
mu_alpha = 0
sigma_sq_alpha = 1000
K = dim(adjacency_matrix)[1]
W = adj
sum_Ws = rowSums(W)

Y = reshape(train[,c('CT2010','School.Year','Count.of.Students.K')],
            timevar = "School.Year",
            idvar = 'CT2010',
            direction = "wide")
Y = Y[ , !(names(Y) %in% 'CT2010')]
rownames(Y) = train$CT2010[1:dim(Y)[1]]

T_train = dim(Y)[2]
Time = seq(1, T_train)
Time = (Time - mean(Time))/T_train

init_phis = rep(0,K)
init_deltas = rep(0,K)
for (i in 1:K){
  data_list = list(Time = Time, Count = as.numeric(Y[i,]))
  uncorrelated_params = lm(formula = Count ~ Time, data = data_list)$coefficients
  uncorrelated_params = as.numeric(uncorrelated_params)
  init_phis[i] = uncorrelated_params[1]
  init_deltas[i] = uncorrelated_params[2]
}
init_alpha = mean(init_deltas)
init_deltas = init_deltas - init_alpha
init_beta = mean(init_phis)
init_phis = init_phis - init_beta

init_psi = ((matrix(init_beta, K, T_train) + matrix(init_phis,nrow=K,ncol=T_train))
            + (matrix(init_alpha, K, T_train) +
                 matrix(Time, nrow=K, ncol=T_train, byrow=TRUE)*matrix(init_deltas,nrow=K, ncol=T_train)))

smart_init = list(alpha = init_alpha,
                  beta = init_beta,
                  phi_raw = init_phis,
                  delta_raw = init_deltas,
                  psi = init_psi,
                  phi_mean = init_phis,
                  delta_mean = init_deltas,
                  phi = init_phis,
                  delta=init_deltas)

data = list(a=a, b=b, mu_alpha=mu_alpha, sigma_sq_alpha=sigma_sq_alpha, K=K, W=W, sum_Ws=sum_Ws, Y=Y, T_train=T_train, Time=Time)
pars = c("beta","phi","alpha",'delta','tau_sq_int','tau_sq_slo','rho_int','rho_slo','v_sq', 'psi')

#Define initialization
#Note this follows the approach taken in the packaged gibbs sampler
#See lines 112-127 of https://github.com/cran/CARBayesST/blob/master/R/gaussian.CARlinear.R

stan_model = stan(file="lin_model.stan",
                  data = data,
                  pars = pars,
                  warmup = warmup,
                  iter = num_iter,
                  init = rep(list(smart_init), num_chains),
                  refresh = 10,
                  chains = num_chains)

#save fit stan model:
save(stan_model, file = file_string)

library(Metrics)

test_split = split(test, test$CT2010)
for (i in 1:length(test_split)){
  test_split[[i]] = test_split[[i]][,'Count.of.Students.K']
}

test_Time = ((T_train+1):(T_train+2) - mean(1:T_train))/T_train

coeffs_for_model = c('beta','phi','alpha','delta')
post_mean = get_posterior_mean(stan_model, pars=coeffs_for_model)[,5] #All chains

alpha_coef = post_mean['alpha']
beta_coef = post_mean['beta']

make_bayesian_predict = function(i, test_data_list){
  vec_time = test_data_list['vec_time'][[1]]
  phi_coef = post_mean[sprintf('phi[%d]',i)]
  delta_coef = post_mean[sprintf('delta[%d]',i)]
  return((beta_coef + phi_coef) +
           (alpha_coef + delta_coef)*vec_time )
}

bayes_se_train = list()
baseline_se_train = list()
bayes_se_test = list()
baseline_se_test = list()

for (i in 1:K){
  train_data_list = list(vec_time = Time,
                         vec_counts = as.numeric(Y[i,]))
  test_data_list = list(vec_time = test_Time,
                        vec_counts = as.numeric(test_split[[i]]))
  baseline_model_K = lm(formula = vec_counts ~ vec_time, data = train_data_list)

  baseline_prediction_train = as.numeric(predict(baseline_model_K, train_data_list))
  baseline_prediction_test = as.numeric(predict(baseline_model_K, test_data_list))
  bayes_prediction_train = make_bayesian_predict(i,train_data_list)
  bayes_prediction_test = make_bayesian_predict(i,test_data_list)

  bayes_se_train[[i]] = se(as.vector(Y[i,]), bayes_prediction_train)
  bayes_se_test[[i]] = se(as.vector(test_split[[i]]), bayes_prediction_test)
  baseline_se_train[[i]] = se(as.vector(Y[i,]), baseline_prediction_train)
  baseline_se_test[[i]] = se(as.vector(test_split[[i]]), baseline_prediction_test)
}

se_to_RMSE = function(x) sqrt(sum(as.numeric(lapply(x, sum)))/(length(x)*length(x[[1]])))

results = rbind(c(se_to_RMSE(bayes_se_train), se_to_RMSE(bayes_se_test)),
                c(se_to_RMSE(baseline_se_train),se_to_RMSE(baseline_se_test)))
rownames(results) = c('Bayes',"Baseline")
colnames(results) = c("Train",'Test')
write.table(results, './stan_models/lin_model_grade_K_scores.csv',sep=',')
