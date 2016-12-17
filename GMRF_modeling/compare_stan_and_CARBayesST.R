source('get_data.R')
library(rstan)
library(coda)
library(matrixStats)
library(CARBayesST)

load("./stan_models/lin_model_stan_fit.RData")

load(file = "./stan_models/lin_model_CARBayesST_fit.RData")

#compare two models using metric from CARBayesST
#Note CARBayesST logLik() returns the following estimate for logLik:
#(from https://github.com/cran/CARBayesST/blob/master/R/gaussian.CARlinear.R)
# deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),N.all), log = TRUE), na.rm=TRUE)
# (-0.5 * deviance.fitted)

#However private correspondence indicated bug in function...

logLik_stan = function(stan_model){
  v_sq_sd = sqrt(median(extract(stan_model, v_sq)))
  v_sq_sd = matrix(v_sq_sd, K, T_train)
  median_means = lapply(stan_model, 'psi'), median)
  median_means_matrix = matrix(unlist(median_means), byrow=FALSE, nrow=K)
  stan_deviance_fitted = 0
  for (i in 1:K){
    stan_deviance_fitted = stan_deviance_fitted -2*sum(dnorm(as.numeric(Y[i,]),mean = median_means_matrix[i,], sd= v_sq_sd[i,], log=TRUE), na.rm=TRUE)
  }
  return(-0.5 * stan_deviance_fitted)
}

logLik_stan(stan_model)

model_pars = c('beta','alpha','phi','delta','nu2')
gibbs_model_medians = lapply(modelK_spat$samples[model_pars],colMedians)
logLik_gibbs_from_scratch = function(gibbs_model_medians){
  v_sq_sd = matrix(gibbs_model_medians$nu2, nrow=K, ncol=T_train)
  median_means = (matrix(rep(gibbs_model_medians$beta, K) + gibbs_model_medians$phi,nrow=K, ncol=T_train,byrow=FALSE) +
                    Time*matrix(rep(gibbs_model_medians$alpha, K)+gibbs_model_medians$delta,nrow=K, ncol=T_train, byrow=FALSE))
  gibbs_deviance_fitted = 0
  for (i in 1:K){
    gibbs_deviance_fitted = gibbs_deviance_fitted -2*sum(dnorm(as.numeric(Y[i,]),mean = median_means[i,], sd= v_sq_sd[i,], log=TRUE), na.rm=TRUE)
  }
  return(-0.5 * gibbs_deviance_fitted)
}

logLik_gibbs_from_scratch(gibbs_model_medians)
logLik(modelK_spat) #Note there is a bug in this function...

moving_avg = function(x) {
  post_warmup = x[-2000:-1]
  return(cumsum(post_warmup)/seq(along=post_warmup))
}
stan_moving_avgs = lapply(stan_model@sim$samples[[1]], moving_avg)

plot((post_means_modelK_spat$beta - stan_moving_avgs$beta)/post_means_modelK_spat$beta,
     type='l',
     main = 'Convergence of stan HMC sampling for beta',
     xlab='Number of MCMC samples',
     ylab='% error ((Gibbs value vs moving MCMC mean)')

post_beta= As.mcmc.list(stan_model,pars="beta")
plot(post_beta)

plot((post_means_modelK_spat$alpha - stan_moving_avgs$alpha)/(post_means_modelK_spat$alpha),
     type='l',
     main = 'Convergence of stan HMC sampling for alpha',
     xlab='Number of MCMC samples',
     ylab='% error ((Gibbs value vs moving MCMC mean)')

post_alpha= As.mcmc.list(stan_model,pars="alpha")
plot(post_alpha, main="Sampling behavior for alpha")

post_means_stan_model = get_posterior_mean(stan_model, pars=pars)
cor(post_means_modelK_spat$phi, post_means_stan_model[grep('phi', rownames(post_means_stan_model))])
cor(post_means_modelK_spat$delta, post_means_stan_model[grep('delta', rownames(post_means_stan_model))])

delta_MAs = as.data.frame(stan_moving_avgs[grep('delta', names(stan_moving_avgs))])
phi_MAs = as.data.frame(stan_moving_avgs[grep('phi', names(stan_moving_avgs))])
plot(apply(delta_MAs, 1, function(x) {cor(x, post_means_modelK_spat$delta)}),
     type='l',
     xlab='Number of MCMC samples',
     ylab='Correlation(MCMC deltas, Gibbs deltas)',
     main='Convergence of deltas (posterior means)')

hist(as.numeric(phi_MAs[2000,]),
     main='Histogram of phi samples, HMC',
     breaks=20,
     xlab="Posterior means of phi's")

hist(as.numeric(delta_MAs[2000,]),
     main='Histogram of delta samples, HMC',
     breaks=20,
     xlab="Posterior means of deltas's")

plot(apply(phi_MAs, 1, function(x) {cor(x, post_means_modelK_spat$phi)}),
     type='l',
     xlab='Number of MCMC samples',
     ylab='Correlation(HMC phis, Gibbs phis)',
     main='Convergence of phis (posterior means)')

stan_rhat(stan_model)
stan_diag(stan_model, 
          information = c("stepsize"), 
          chain = 0)