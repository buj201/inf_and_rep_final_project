setwd("~/Desktop/Other/MS_Courses/Inf_and_rep/project/nycSchoolPredictions/GMRF_modeling")
source('get_data.R')

library(rstan)
options(mc.cores = parallel::detectCores())

test_flag = FALSE

if (test_flag){
  num_iter = 200
  warmup = 100
  num_chains = 1
  file_string = "./stan_models/multivar_model_quadratic_stan_fit_test.RData"
} else {
  num_iter = 1500
  warmup = 500
  num_chains = 4
  file_string = "./stan_models/multivar_model_quadratic_stan_fit.RData"
}

K = dim(adj)[1]
G = 6
T_train = 8

train_split = split(train, train$CT2010)
for (i in 1:length(train_split)){
  train_split[[i]] = train_split[[i]][,!(names(train_split[[i]]) %in% c('CT2010','School.Year'))]
  rownames(train_split[[i]]) = c('20012002','20022003','20032004','20042005','20052006','20062007','20072008','20082009')
  train_split[[i]] = t(train_split[[i]])
  train_split[[i]] = train_split[[i]][c(6,1,2,3,4,5),]
  train_split[[i]] = matrix(train_split[[i]],nrow(train_split[[i]]), ncol(train_split[[i]]))
}

#See https://groups.google.com/forum/#!topic/stan-users/-6s9PvVsyug
Y_matrices = array(NA_real_, dim = c(K,G,T_train))
for (i in 1:K){
  Y_matrices[i,,] = train_split[[i]]
}


a = 1
b = 0.01
mu_alpha = 0
sigma_sq_alpha = 1000
Y = Y_matrices
W = adj
sum_Ws = rowSums(W)
Time = matrix((1:T_train - mean(1:T_train))/T_train, G, T_train, byrow = TRUE)
Grade = matrix((1:G - mean(1:G))/G, G, T_train, byrow = FALSE)

init_beta_0 = rep(0,K)
init_beta_g = rep(0,K)
init_beta_t = rep(0,K)
init_beta_gt = rep(0,K)
init_beta_g2 = rep(0,K)
init_beta_t2 = rep(0,K)
for (i in 1:K){
  data_list = list(vec_time = as.vector(Time),
                   vec_grade = as.vector(Grade),
                   vec_counts = as.vector(train_split[[i]]))
  uncorrelated_params = lm(formula = vec_counts ~ vec_grade + vec_time + vec_grade*vec_time + I(vec_grade^2) + I(vec_time^2),
                           data = data_list)
  uncorrelated_params = as.numeric(uncorrelated_params$coefficients)
  init_beta_0[i] = uncorrelated_params[1]
  init_beta_g[i] = uncorrelated_params[2]
  init_beta_t[i] = uncorrelated_params[3]
  init_beta_g2[i] = uncorrelated_params[4]
  init_beta_t2[i] = uncorrelated_params[5]
  init_beta_gt[i] = uncorrelated_params[6]
}
init_alpha_0 = mean(init_beta_0)
init_beta_0 = init_beta_0 - init_alpha_0

init_alpha_g = mean(init_beta_g)
init_beta_g = init_beta_g - init_alpha_g

init_alpha_t = mean(init_beta_t)
init_beta_t = init_beta_t - init_alpha_t

init_alpha_gt = mean(init_beta_gt)
init_beta_gt = init_beta_gt - init_alpha_gt

init_alpha_g2 = mean(init_beta_g2)
init_beta_g2 = init_beta_g2 - init_alpha_g2

init_alpha_t2 = mean(init_beta_t2)
init_beta_t2 = init_beta_t2 - init_alpha_t2

init_mu = list()
for (i in 1:K){
  init_mu[[i]] = (init_alpha_0 + init_beta_0[i]) +
            (init_alpha_g + init_beta_g[i])*Time +
            (init_alpha_g + init_beta_g[i])*Grade +
            (init_alpha_gt + init_beta_gt[i]) * Grade * Time + #* is elementwise in R
            (init_alpha_g2 + init_beta_g2[i]) * Grade * Grade + #* is elementwise in R
            (init_alpha_t2 + init_beta_t2[i]) * Time * Time#* is elementwise in R
}
init_mu_matrices = array(NA_real_, dim = c(K,G,T_train))
for (i in 1:K){
  init_mu_matrices[i,,] = init_mu[[i]]
}


smart_init = list(alpha_0 = init_alpha_0,
                  alpha_g = init_alpha_g,
                  alpha_t = init_alpha_t,
                  alpha_gt = init_alpha_gt,
                  alpha_g2 = init_alpha_g2,
                  alpha_t2 = init_alpha_t2,
                  beta_0_raw = init_beta_0,
                  beta_0 = init_beta_0,
                  beta_g_raw = init_beta_g,
                  beta_g = init_beta_g,
                  beta_t_raw = init_beta_t,
                  beta_t = init_beta_t,
                  beta_gt_raw = init_beta_gt,
                  beta_gt = init_beta_gt,
                  beta_g2_raw = init_beta_g2,
                  beta_g2 = init_beta_g2,
                  beta_t2_raw = init_beta_t2,
                  beta_t2 = init_beta_t2,
                  mu = init_mu_matrices)

data = list(a=a,
            b=b,
            mu_alpha=mu_alpha,
            sigma_sq_alpha=sigma_sq_alpha,
            K=K,
            G=G,
            T_train=T_train,
            Y=Y,
            W=W,
            sum_Ws=sum_Ws,
            Time=Time,
            Grade=Grade)

pars = c("alpha_0","alpha_g","alpha_t",'alpha_gt','alpha_g2','alpha_t2',
         "beta_0","beta_g","beta_t",'beta_gt','beta_g2','beta_t2',
         'tau_sq_0','tau_sq_g','tau_sq_t','tau_sq_gt','tau_sq_t2','tau_sq_g2',
         "rho_0","rho_g","rho_t",'rho_gt',"rho_t2",'rho_g2',
         'v_sq', 'mu')

multi_grade_quadratic_stan_model = stan(file="multi_grade_model_quadratic.stan",
                  data = data,
                  pars = pars,
                  warmup = warmup,
                  iter = num_iter,
                  init = rep(list(smart_init), num_chains),
                  refresh = 10,
                  chains = num_chains)

#save fit stan model:
save(multi_grade_quadratic_stan_model, file = file_string)

#Make predicitions

library(Metrics)

test_split = split(test, test$CT2010)
for (i in 1:length(test_split)){
  test_split[[i]] = test_split[[i]][,!(names(test_split[[i]]) %in% c('CT2010','School.Year'))]
  rownames(test_split[[i]]) = c('20092010','20102011')
  test_split[[i]] = t(test_split[[i]])
  test_split[[i]] = test_split[[i]][c(6,1,2,3,4,5),]
  test_split[[i]] = matrix(test_split[[i]],nrow(test_split[[i]]), ncol(test_split[[i]]))
}

test_Time = matrix(((T_train+1):(T_train+2) - mean(1:T_train))/T_train, G, 2, byrow = TRUE)
test_Grade = matrix((1:G - mean(1:G))/G, G, 2, byrow = FALSE)

coeffs_for_model = c('beta_0','beta_t',
                     'beta_g','beta_gt',
                     'beta_g2','beta_t2',
                     'alpha_0','alpha_t',
                     'alpha_g','alpha_gt',
                     'alpha_g2','alpha_t2')
post_mean = get_posterior_mean(multi_grade_quadratic_stan_model, pars=coeffs_for_model)[,5]#all chains

alpha_0_coef = post_mean['alpha_0']
alpha_t_coef = post_mean['alpha_t']
alpha_g_coef = post_mean['alpha_g']
alpha_gt_coef = post_mean['alpha_gt']
alpha_g2_coef = post_mean['alpha_g2']
alpha_t2_coef = post_mean['alpha_t2']

make_bayesian_predict = function(i, test_data_list){
  vec_time = test_data_list['vec_time'][[1]]
  vec_grade = test_data_list['vec_grade'][[1]]
  beta_0_coef = post_mean[sprintf('beta_0[%d]',i)]
  beta_t_coef = post_mean[sprintf('beta_t[%d]',i)]
  beta_g_coef = post_mean[sprintf('beta_g[%d]',i)]
  beta_gt_coef = post_mean[sprintf('beta_gt[%d]',i)]
  beta_t2_coef = post_mean[sprintf('beta_t2[%d]',i)]
  beta_g2_coef = post_mean[sprintf('beta_g2[%d]',i)]
  return((alpha_0_coef + beta_0_coef) +
         (alpha_t_coef + beta_t_coef)*vec_time +
         (alpha_g_coef + beta_g_coef)*vec_grade +
         (alpha_gt_coef + beta_gt_coef)*vec_time*vec_grade +
         (alpha_g2_coef + beta_g2_coef)*vec_grade*vec_grade +
         (alpha_t2_coef + beta_t2_coef)*vec_time*vec_time
        )
}

bayes_se_train = list()
baseline_se_train = list()
bayes_se_test = list()
baseline_se_test = list()

for (i in 1:K){
  train_data_list = list(vec_time = as.vector(Time),
                   vec_grade = as.vector(Grade),
                   vec_counts = as.vector(train_split[[i]]))
  test_data_list = list(vec_time = as.vector(test_Time),
                        vec_grade = as.vector(test_Grade))
  baseline_model_K = lm(formula = vec_counts ~ vec_grade + vec_time + vec_grade*vec_time + I(vec_grade^2) + I(vec_time^2),
                           data = train_data_list)

  baseline_prediction_train = as.numeric(predict(baseline_model_K, train_data_list))
  baseline_prediction_test = as.numeric(predict(baseline_model_K, test_data_list))
  bayes_prediction_train = make_bayesian_predict(i,train_data_list)
  bayes_prediction_test = make_bayesian_predict(i,test_data_list)

  bayes_se_train[[i]] = se(as.vector(train_split[[i]]), bayes_prediction_train)
  bayes_se_test[[i]] = se(as.vector(test_split[[i]]), bayes_prediction_test)
  baseline_se_train[[i]] = se(as.vector(train_split[[i]]), baseline_prediction_train)
  baseline_se_test[[i]] = se(as.vector(test_split[[i]]), baseline_prediction_test)
}

se_to_RMSE = function(x) sqrt(sum(as.numeric(lapply(x, sum)))/(length(x)*length(x[[1]])))

results = rbind(c(se_to_RMSE(bayes_se_train), se_to_RMSE(bayes_se_test)),
                c(se_to_RMSE(baseline_se_train),se_to_RMSE(baseline_se_test)))
rownames(results) = c('Bayes',"Baseline")
colnames(results) = c("Train",'Test')
write.table(results, './stan_models/multi_grade_model_quadratic_scores.csv',sep=',')
