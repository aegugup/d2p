library(rstan)
library(loo) # compute waic and looic
library(coda)
library(bayesplot)
library(bayesutils)


source("dataPrep.R")
total_pfas_df = get_pfc()
acs_df = get_acs(total_pfas_df$COUNTY)
total_pfas_df = cbind(total_pfas_df, acs_df)

demos = c("Female",
          "35 years and over",
          "Black or African American",
          "American Indian and Alaska Native",
          "Asian",
          "Native Hawaiian and Other Pacific Islander",
          "Some other race",
          "Hispanic or Latino (of any race)")

pfas_vars = c("AREA_KM2", "IOI_COUNT", "NUM_IOI_PER_KM2", 
              "BOOL_DISCHARGE", "NUM_DISCHARGE", "AVG_AVG_CONC_mgL" ,                         
              "BOOL_FED_SITES", "NUM_FED_SITES",                             
              "BOOL_BIOSOLIDS", "NUM_BIOSOLIDS","NUM_BIOSOLIDS_PER_KM2")

counties = unique(total_pfas_df$COUNTY)

for (i in 1:length(counties)){
  cty = (total_pfas_df$COUNTY == counties[i]) + 0
  total_pfas_df[,counties[i]] = cty
}

# --------------------------------- Model 1 --------------------------------- #
sl_stanmod = "
data {
  int<lower=1> n; // number of observations
  int<lower=1> m; // number of covariates
  vector[n] y; // response
  matrix[n,m] X; // covariates - includes beta0
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
  mu = exp(X * beta + sigmasq/2);
}
model {
  // specify priors
  for(i in 1:m){
    beta[i] ~ normal(mu0[i], sqrt(sigmasq0[i]));
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

pfas_vars = c("AREA_KM2", "IOI_COUNT", "NUM_IOI_PER_KM2", 
              "BOOL_DISCHARGE", "NUM_DISCHARGE", "AVG_AVG_CONC_mgL" ,                         
              "BOOL_FED_SITES", "NUM_FED_SITES",                             
              "BOOL_BIOSOLIDS", "NUM_BIOSOLIDS","NUM_BIOSOLIDS_PER_KM2")

mod1_vars = c(demos, "NUM_IOI_PER_KM2")
mod1_x = cbind(x0 = rep(1, each = length(total_pfas_df$CONCENTRATION_ngL)),
               total_pfas_df[,mod1_vars])

mod2_vars = c(demos, "NUM_IOI_PER_KM2", "NUM_DISCHARGE")
mod2_x = cbind(x0 = rep(1, each = length(total_pfas_df$CONCENTRATION_ngL)),
               total_pfas_df[,mod2_vars])

mod3_vars = c(demos, "NUM_IOI_PER_KM2", "NUM_DISCHARGE", "NUM_FED_SITES")
mod3_x = cbind(x0 = rep(1, each = length(total_pfas_df$CONCENTRATION_ngL)),
               total_pfas_df[,mod3_vars])

mod4_vars = c(demos, "NUM_IOI_PER_KM2", "NUM_DISCHARGE", "NUM_FED_SITES", 
              "NUM_BIOSOLIDS_PER_KM2")
mod4_x = cbind(x0 = rep(1, each = length(total_pfas_df$CONCENTRATION_ngL)),
               total_pfas_df[,mod4_vars])


all_mod_x = list(mod1_x, mod2_x, mod3_x, mod4_x)


get_fit = function(i_mod, mod_x){
  
  mod_data = list(n = length(total_pfas_df$CONCENTRATION_ngL),
                  m = ncol(mod_x),
                  y = total_pfas_df$CONCENTRATION_ngL,
                  X = mod_x,
                  log_population = log(total_pfas_df$`Total population`),
                  mu0 = numeric(length = ncol(mod_x)),
                  sigmasq0 = numeric(length = ncol(mod_x)) + 1,
                  v = var(total_pfas_df$CONCENTRATION_ngL))
  
  mod = stan_model(model_code = stanmod)
  # draw samples from the model
  fit = sampling(mod, 
                 data = mod_data,
                 iter = 1e5, 
                 cores = 4)
  
  # extract log likelihoods
  ll = extract_log_lik(fit, merge_chains = FALSE)
  
  # compute relative efficiency of log likelihoods
  rel_eff = exp(relative_eff(ll))
  
  # compute looic for each model
  (mod_looic = loo(ll, r_eff = rel_eff))
  save(mod_looic, file=sprintf("mod%d_looic.rda",i_mod))
  
  # compute waic on both fitted models
  (mod_waic = waic(ll))
  save(mod_waic, file = sprintf("mod%d_waic.rda",i_mod))
  
  # posterior means, sd, credible intervals, gelman-rubin
  pars = c("beta", "exp_beta")
  mod_summary = summary(fit, 
                        pars = pars, 
                        probs = c(0.025, 0.975), 
                        use_cache = FALSE)$summary
  
  print(mod_summary)
  
  save(mod_summary, file = sprintf("mod%d_summary.rda", i_mod))
}

mod_data = list(n = length(total_pfas_df$CONCENTRATION_ngL),
                m = ncol(mod1_x),
                y = total_pfas_df$CONCENTRATION_ngL,
                log_population = log(total_pfas_df$`Total population`),
                X = mod1_x,
                mu0 = numeric(length = ncol(mod1_x)),
                sigmasq0 = numeric(length = ncol(mod1_x)) + 1,
                v = var(total_pfas_df$CONCENTRATION_ngL))

mod = stan_model(model_code = stanmod)
# draw samples from the model
fit = sampling(mod, 
               data = mod_data,
               iter = 1e5,
               # control = list(adapt_delta=0.99),
               cores = 4)

pars = c("beta", "exp_beta")
mod_summary = summary(fit, 
                      pars = pars, 
                      probs = c(0.025, 0.975), 
                      use_cache = FALSE)$summary
mod_summary

pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]",
         "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]")

mcmc_parcoord(as.array(fit), np = nuts_params(fit), pars = pars)

mcmc_pairs(as.array(fit), np = nuts_params(fit), pars = c("beta[1]", "beta[2]", "beta[3]"),
           off_diag_args = list(size = 0.75))


# extract log likelihoods
ll = extract_log_lik(fit, merge_chains = FALSE)

# compute relative efficiency of log likelihoods
rel_eff = exp(relative_eff(ll))

# compute looic for each model
(mod_looic = loo(ll, r_eff = rel_eff))

# compute waic on both fitted models
(mod_waic = waic(ll))

for (i in 2:length(all_mod_x)){
  get_fit(i, all_mod_x[[i]])
}

# Heidelberg & Welch: If the halfwidth test fails, extend the chain(s)
coda::heidel.diag(fit)
# trace plots of results
mcmc_trace(fit, pars = pars)

# plot of acf of chains
stan_ac(fit, "beta")
stan_ac(fit, "sigmasq")

# ---------------- Model Checking
# length, mean, and sd of data
fit = as.data.frame(fit)
i_mu = which(!is.na(str_locate(colnames(fit),"mu")[c(1:ncol(fit)),1]))
i_sigma = 1
B = 1e5
y = log(total_pfas_df$CONCENTRATION_ngL)
n = length(y)
m = mean(y)
s = sd(y)

# store nyrep samples of yrep from from posterior predictive distribution
nyrep = 1e4
yrep = matrix(0, nrow = nyrep, ncol = n)
loglik_yrep = matrix(0, nrow = nyrep, ncol = n)
for (i in seq_len(nyrep)) {
  for (j in seq_len(n)){
    mu = fit[i,i_mu[j]]
    sigmasq = fit[i, i_sigma]
    # yrep[i, j] = rlnorm(1, meanlog = mu, sdlog = sqrt(sigmasq))
    # loglik_yrep[i,j] = dlnorm(yrep[i,j],
    #                           meanlog = mu,
    #                           sdlog = sqrt(sigmasq),
    #                           log = TRUE)
    yrep[i, j] = rnorm(1, mean = exp(mu), sd = sqrt(sigmasq))
    loglik_yrep[i,j] = dnorm(yrep[i,j],
                             mean = exp(mu),
                             sd = sqrt(sigmasq),
                             log = TRUE)
  }
}

plot(density(y))
lines(density(yrep[4,]))
# compare minimum of observed data to minimum from
# replicated samples

# minimum of replicated sample
mins = apply(yrep, 1, min)
# estimated p-value
(sum(mins <= min(y + 1)/(length(mins) + 1)))

# histogram comparing T(y) and T(yrep)
ppc_stat(y, yrep, stat = "min")

# look at asymmetry of distribution
# by comparing order statistics to
# samples of posterior mean
d_sim = d_obs =  matrix(0, nrow = nyrep, ncol = n)
sort_y = sort(y)
orders = quantile(y, probs = c(0.10, 0.90))
order_10 = which(sort_y >= orders[1])[1]
order_90 = which(sort_y >= orders[2])[1]
for (i in 1:nrow(yrep)) {
  sort_yrep = sort(yrep[i, ])
  for (j in 1:ncol(yrep)){
    mu = exp(fit[i,i_mu[j]])
    d_sim[i, j] = abs(sort_yrep[order_90] - mu) -
      abs(sort_yrep[order_10] - mu)
    d_obs[i, j] = abs(sort_y[order_90] - mu) -
      abs(sort_y[order_10] - mu)
  }
}
# estimated posterior predictive p-value
(sum(d_sim[,1] >= d_obs[,1]) + 1)/(nrow(d_sim) + 1)

# compare observed and simulated discrepancy measures
plot(d_sim[,1] ~ d_obs[,1],
     xlab = "observed discrepancy",
     ylab = "simulated discrepancy")
abline(0, 1)

# compare histograms of y and yrep
ppc_hist(y, yrep[sample(nyrep, 8), ])
# compare boxplots of y and yrep
ppc_boxplot(y, yrep[sample(nyrep, 8), ])
# compare densities of y and yrep
ppc_dens_overlay(y, yrep[sample(nyrep, 20), ])
# compare ecdfs of y and yrep
ppc_ecdf_overlay(y, yrep[sample(nyrep, 20) , ])
# compare histograms of y - yrep
ppc_error_hist(y, yrep[sample(nyrep, 9) , ])
# compare scatterplots of y vs yrep
ppc_scatter(y, yrep[sample(nyrep, 9) , ])
# compare scatterplots of y vs y - yrep
ppc_error_scatter(y, yrep[sample(nyrep, 9) , ])

# marginal predictive checks
# comparison of observed data and 90% predictive cpis
ppc_intervals(y, yrep)
ppc_ribbon(y, yrep)


# compute relative effective MCMC sample size
# divided by the total sample size for likelihood
r_eff = relative_eff(exp(loglik_yrep),
                     chain_id = rep(1, nyrep))
# compute leave-one-out information
loo_info = loo(loglik_yrep, r_eff = r_eff, save_psis = TRUE)

# construct leave-one-out prediction intervals
ppc_loo_intervals(y, yrep,
                  psis_object = loo_info$psis_object)

# construct leave-one-out quantiles
ppc_loo_pit_qq(y, yrep,
               lw = weights(loo_info$psis_object))

# ppo vs y
PPO = colMeans(exp(loglik_yrep))
plot(PPO ~ y, ylab = "PPO")
dy = density(y)
# scale density to match scale of PPO
dy$y = dy$y/max(dy$y)*(max(PPO) - min(PPO)) + min(PPO)
lines(dy)


mean = exp(fit[i,i_mu[1]])
sigmasq = fit[i, i_sigma]
location <- log(mean^2 / sqrt(sigmasq + mean^2))
shape <- sqrt(log(1 + (sigmasq / mean^2)))
set.seed(123)
samples1 = rlnorm(1000, meanlog = location, sdlog = shape)
samples2 = rnorm(1000, mean = mean, sd = sqrt(sigmasq))

