source('get_data.R')

library(CARBayesST)

test_flag = TRUE

if (test_flag){
  burnin = 1000
  n.sample=5000
  file_string = "./stan_models/lin_model_CARBayesST_fit_test.RData"
} else {
  burnin = 20000
  n.sample= 60000
  file_string = "./stan_models/lin_model_CARBayesST_fit.RData"
}

# From CRAN vignette- The model can be fitted with the following one-line function call, and we note
# that all data vectors (response, offset and covariates) have to be ordered so that the first K
# data points relate to all spatial units at time 1, the next K data points to all spatial units at
# time 2 and so on.

modelK_spat = ST.CARlinear(formula = Count.of.Students.K ~ 1,
                           family='gaussian',
                           data = train,
                           W = adj,
                           burnin = 1000, #set to 20000
                           n.sample=5000) #set to 60000

hist(modelK_spat$samples$alpha,
     main='Histogram of alpha samples',
     breaks=20,
     xlab='Distribution of MCMC alpha samples')

hist(modelK_spat$samples$beta,
     main='Histogram of beta samples',
     breaks=20,
     xlab='Distribution of MCMC beta samples')

hist(colMeans(modelK_spat$samples$phi),
     main='Histogram of phi samples, Gibbs',
     breaks=20,
     xlab="Posterior means of phi's")

hist(colMeans(modelK_spat$samples$delta),
     main='Histogram of delta samples',
     breaks=20,
     xlab="Posterior means of delta's")

save(modelK_spat, file = "./stan_models/lin_model_CARBayesST_fit_test.RData")