library(rstan)
library(loo) # compute waic and looic
library(coda)
library(bayesplot)
library(bayesutils)
library(xtable)




# --------------------------------- Model 1 --------------------------------- #
stanmod = "
data {
  int<lower=1> n; // number of observations
  int<lower=1> m; // number of covariates
  vector[n] y; // response
  matrix[n,m] X; // covariates - includes beta0
  vector[n] log_population; // population offset
  vector[m] mu0; // prior mean for betas
  vector[m] sigmasq0; // prior variance for betas
  real<lower=0> v; // sample variance of y
}
parameters {
  real<lower=0> sigmasq;
  vector[m] beta;
}
transformed parameters{
  vector[n] mu;  // mean of responses
  vector[n] logmu;  // mean of responses
  vector[n] logsigmasq;  // log likelihood of data
  mu = X * beta;
  logmu = exp(mu + sigmasq/2);
  logsigmasq = exp(2 * mu + sigmasq) * (exp(sigmasq) - 1);
}
model {
  // specify priors
  for(i in 1:m){
    beta[i] ~ normal(mu0[i], sqrt(sigmasq0[i]));
  }
  sigmasq ~ inv_gamma(0.01, 0.01);

  // data distribution
  for(i in 1:n){
    y[i] ~ lognormal(logmu[i], sqrt(logsigmasq[i]));
  }
}
generated quantities {
  real Rbsq; // goodness-of-fit
  vector[m] exp_beta; // transform beta
  vector[n] log_lik;  // log likelihood of data
  
  Rbsq = 1 - sigmasq/v;
  for (i in 1:m) exp_beta[i] = exp(beta[i]);
  for (i in 1:n){
    log_lik[i] = lognormal_lpdf(y[i]| logmu[i], sqrt(logsigmasq[i]));
  } 
}
"

# --------------------------------- Model 2 --------------------------------- #
sl_stanmod = "
data {
  int<lower=1> n; // number of observations
  int<lower=1> m; // number of covariates
  int<lower=1> q; // number of interactions
  
  vector[n] y; // response
  matrix[n,m] X; // covariates - includes beta0
  vector[q,2] I; // which interactions indices - max two per interaction
  vector[m] mu_alpha; // prior mean for alphas
  vector[m] sigmasq_alpha; // prior variance for alphas
  vector[m] mu_beta; // prior mean for betas
  vector[m] sigmasq_beta; // prior variance for betas
  real<lower=0> v; // sample variance of y
}
parameters {
  real<lower=0> sigmasq;
  vector[m] beta;
  vector[q] alpha;
}
transformed parameters{
  vector[n] mu;  // mean of responses
  vector[q] interactions;

  for(i in 1:q){
    interactions[i] = X[,I[i,1]] * X[,I[i,2]];
  }
  mu = exp(X * beta + alpha * interactions + sigmasq/2);
}
model {
  // specify priors
  for(i in 1:m){
    beta[i] ~ normal(mu_beta[i], sqrt(sigmasq_beta[i]));
  }
  for(i in 1:q){
    alpha[i] ~ normal(mu_alpha[i], sqrt(sigmasq_alpha[i]));
  }
  sigmasq ~ inv_gamma(0.01, 0.01);

  // data distribution
  for(i in 1:n){
    y[i] ~ lognormal(mu[i], sqrt(sigmasq));
  }
}
generated quantities {
  real Rbsq; // goodness-of-fit
  vector[m] exp_beta; // transform beta
  vector[n] log_lik;  // log likelihood of data
  
  Rbsq = 1 - sigmasq/v;
  for (i in 1:m) exp_beta[i] = exp(beta[i]);
  for (i in 1:n) log_lik[i] = lognormal_lpdf(y[i]| mu[i], sqrt(sigmasq));
}
"





