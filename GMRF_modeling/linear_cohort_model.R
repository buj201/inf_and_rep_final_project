setwd("~/Desktop/Other/MS_Courses/Inf_and_rep/project/nycSchoolPredictions/GMRF_modeling")
source('get_data.R')

library(rstan)
options(mc.cores = parallel::detectCores())

## Set up parameters and data for Stan

test_flag = FALSE

if (test_flag){
  num_iter = 1500
  warmup = 500
  num_chains = 1
} else {
  num_iter = 5000
  warmup = 1000
  num_chains = 4
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

# Learn two models- correlated and uncorrelated

corr_pars = c('rate', 'rho', 'tau_sq','v_sq')
corr_data = list(a=a, b=b, mu_alpha=mu_alpha, sigma_sq_alpha=sigma_sq_alpha, G=G, K=K, W=W, sum_Ws=sum_Ws, Y=Y, T_train=T_train)

diagonal_change_model_corr = stan(file="linear_cohort_model.stan",
                              data = corr_data,
                              pars = corr_pars,
                              warmup = warmup,
                              iter = num_iter,
                              refresh = 10,
                              chains = num_chains)

no_corr_pars = c('rate', 'v_sq')
no_corr_data = list(a=a, b=b, G=G, K=K, Y=Y, T_train=T_train)

diagonal_change_model_no_corr = stan(file="linear_cohort_model_no_spatial_corr.stan",
                                  data = no_corr_data,
                                  pars = no_corr_pars,
                                  warmup = warmup,
                                  iter = num_iter,
                                  refresh = 10,
                                  chains = num_chains)

#save fit stan models
save(diagonal_change_model_corr, file = "./stan_models/linear_change_corr_stan_fit.RData")
save(diagonal_change_model_no_corr, file = "./stan_models/linear_change_no_corr_stan_fit.RData")

# Set up test data and error matrix function

test_split = split(test, test$CT2010)
for (i in 1:length(test_split)){
  test_split[[i]] = test_split[[i]][,!(names(test_split[[i]]) %in% c('CT2010','School.Year'))]
  rownames(test_split[[i]]) = c('20092010','20102011')
  test_split[[i]] = t(test_split[[i]])
  test_split[[i]] = test_split[[i]][c(6,1,2,3,4,5),]
  test_split[[i]] = matrix(test_split[[i]],nrow(test_split[[i]]), ncol(test_split[[i]]))
}

joined_train_test = list()
for (i in 1:length(test_split)){
  joined_train_test[[i]] = cbind(train_split[[i]],test_split[[i]])
}

delta = row(joined_train_test[[1]]) - col(joined_train_test[[1]])

get_all_prediction_errors = function(data_matrix, rate){
  #Function to compute 6x6 prediction error matrix
  num_predict = matrix(0,6,6)
  num_predict[lower.tri(num_predict, diag=TRUE)] = NA

  sum_SE_by_given_predict = matrix(0,6,6)
  sum_SE_by_given_predict[lower.tri(sum_SE_by_given_predict, diag=TRUE)] = NA

  for (r in 1:(nrow(data_matrix))){
    for (c in 1:(ncol(data_matrix))){
      if (delta[r,c] > 1 | delta[r,c] < -2){ #not used for training
        num_ahead = 1
        while ((r + num_ahead) <= nrow(data_matrix) & (c + num_ahead) <= ncol(data_matrix)){
          sum_SE_by_given_predict[r, r + num_ahead] = sum_SE_by_given_predict[r, r + num_ahead] +
            (data_matrix[r+num_ahead, c+num_ahead] -
               (data_matrix[r,c] + (rate)*num_ahead))^2
          num_predict[r, r + num_ahead] = num_predict[r, r + num_ahead] + 1
          num_ahead = num_ahead + 1
        }
      }
    }
  }
  return(list(pred_sum_SE = sum_SE_by_given_predict, num_pred = num_predict))
}

get_error_matrix = function(model){
  rates = get_posterior_mean(model, 'rate')[,5] #All chains

  error_matrices = list()
  for (i in 1:length(joined_train_test)){
    error_matrices[[i]] =  get_all_prediction_errors(joined_train_test[[i]],rates[i])
  }
  error_matrices = lapply(error_matrices, function(x) x[[1]]/x[[2]])

  sum_squared_error_matrix = matrix(0,6,6)
  sum_squared_error_matrix[lower.tri(sum_squared_error_matrix, diag=TRUE)] = NA

  for (i in 1:length(error_matrices)){
    sum_squared_error_matrix = sum_squared_error_matrix + error_matrices[[i]]
  }
  sum_squared_error_matrix = sum_squared_error_matrix / length(error_matrices)
  RMSE_matrix = sqrt(sum_squared_error_matrix)
  rownames(RMSE_matrix) = c('K',"1",'2','3','4','5')
  colnames(RMSE_matrix) = c('K',"1",'2','3','4','5')
  return(RMSE_matrix)
}

corr_RMSE_errors = get_error_matrix(diagonal_change_model_corr)
no_corr_RMSE_errors = get_error_matrix(diagonal_change_model_no_corr)

#Finally let's just get the simplest baseline- just predict no change...
get_no_change_error = function(){
  predict_no_change_error = function(data_matrix){
    #Function to compute 6x6 prediction error matrix
    num_predict = matrix(0,6,6)
    num_predict[lower.tri(num_predict, diag=TRUE)] = NA

    sum_SE_by_given_predict = matrix(0,6,6)
    sum_SE_by_given_predict[lower.tri(sum_SE_by_given_predict, diag=TRUE)] = NA

    for (r in 1:(nrow(data_matrix))){
      for (c in 1:(ncol(data_matrix))){
        if (delta[r,c] > 1 | delta[r,c] < -2){ #not used for training
          num_ahead = 1
          while ((r + num_ahead) <= nrow(data_matrix) & (c + num_ahead) <= ncol(data_matrix)){
            sum_SE_by_given_predict[r, r + num_ahead] = sum_SE_by_given_predict[r, r + num_ahead] +
              (data_matrix[r+num_ahead, c+num_ahead] -
                 data_matrix[r,c])^2
            num_predict[r, r + num_ahead] = num_predict[r, r + num_ahead] + 1
            num_ahead = num_ahead + 1
          }
        }
      }
    }
    return(list(pred_sum_SE = sum_SE_by_given_predict, num_pred = num_predict))
  }
  no_change_error_matrices = list()
  for (i in 1:length(joined_train_test)){
    no_change_error_matrices[[i]] =  predict_no_change_error(joined_train_test[[i]])
  }
  no_change_error_matrices = lapply(no_change_error_matrices, function(x) x[[1]]/x[[2]])

  sum_squared_error_matrix_no_change = matrix(0,6,6)
  sum_squared_error_matrix_no_change[lower.tri(sum_squared_error_matrix_no_change, diag=TRUE)] = NA

  for (i in 1:length(no_change_error_matrices)){
    sum_squared_error_matrix_no_change = sum_squared_error_matrix_no_change + no_change_error_matrices[[i]]
  }
  sum_squared_error_matrix_no_change = sum_squared_error_matrix_no_change / length(no_change_error_matrices)
  RMSE_matrix_no_change = sqrt(sum_squared_error_matrix_no_change)
  rownames(RMSE_matrix_no_change) = c('K',"1",'2','3','4','5')
  colnames(RMSE_matrix_no_change) = c('K',"1",'2','3','4','5')
  return(RMSE_matrix_no_change)
}

RMSE_matrix_no_change = get_no_change_error()


#Save error matrices
write.table(RMSE_matrix_no_change, './stan_models/no_change_RMSE_error.csv',sep=',')
write.table(corr_RMSE_errors, './stan_models/linear_change_corr_RMSE_error.csv',sep=',')
write.table(no_corr_RMSE_errors, './stan_models/linear_change_no_corr_RMSE_error.csv',sep=',')
