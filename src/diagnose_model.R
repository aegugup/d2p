source("stan_models.R")
source("get_model_fit.R")


# Can do specific diagnostics for fits here prior to running the loop.
fit = get_fit(all_mod_x[[2]])
# posterior means, sd, credible intervals, gelman-rubin
pars = c("beta", "exp_beta", "sigmasq")
mod_summary = summary(fit, 
                      pars = pars, 
                      probs = c(0.025, 0.975), 
                      use_cache = FALSE)$summary
xtable(mod_summary, digits = 4, 
       display=c("s", "f", "e", "e", "f", "f", "d", "f"))
print(mod_summary[,"mean"])
table(t(mod_summary[,"mean"]),)
pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]","beta[6]", 
         "beta[7]", "beta[8]", "beta[9]", "beta[10]","beta[11]",
         "sigmasq")

# Evaluate divergences
mcmc_parcoord(as.array(fit), np = nuts_params(fit), pars = pars)

# Evaluate collinearity and degeneracies: easier to do only a few parameters 
# at a time
mcmc_pairs(as.array(fit), np = nuts_params(fit), 
           pars = c("beta[1]", "beta[2]", "beta[3]"),
           off_diag_args = list(size = 0.75))


# Heidelberg & Welch: If the halfwidth test fails, extend the chain(s)
coda::heidel.diag(fit)

# Trace plots of results
mcmc_trace(fit, pars = pars)

# Plot of autocorrelation of chains
mcmc_acf(fit, pars=pars)


# ---------------- Model Checking
# length, mean, and sd of data
fit = as.data.frame(fit)
i_mu = which(!is.na(str_locate(colnames(fit),"mu")[c(1:ncol(fit)),1]))
i_sigma = 1
i_logsigmasq = which(!is.na(str_locate(colnames(fit),
                                       "logsigmasq")[c(1:ncol(fit)),1]))
B = 1e5
log_y = log(total_pfas_df$CONCENTRATION_ngL)
y = total_pfas_df$CONCENTRATION_ngL
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
    mean = exp(mu + sigmasq/2)
    sd = sqrt(exp(2*mu + sigmasq)*exp(sigmasq - 1))
    yrep[i, j] = rlnorm(1, meanlog = mean, sdlog = sd)
    loglik_yrep[i,j] = dlnorm(yrep[i,j],
                              meanlog = mean,
                              sdlog = sd,
                              log = TRUE)
    
    # yrep[i, j] = rnorm(1, mean = mean, sd = sd)
    # loglik_yrep[i,j] = dnorm(yrep[i,j],
    #                          mean = mean,
    #                          sd = sd,
    #                          log = TRUE)
  }
}
y = log_y
yrep = log(yrep)
plot(density(y))
lines(density(yrep[5,]))

ys <- c(log_y,yrep[7,])
group <- c(rep(TeX("$\\log(y)$"), n), rep(TeX("$\\log(y^{rep})$"), n))
df <- data.frame(ys, group)

ggplot(df, aes(x=ys, fill=group)) + 
  geom_density(alpha=0.6) +
  labs(x="Concentration (ppt)", y="Density", title=title) + 
  theme_bw() 

ggsave("yrep_vs_y_dens.png")
plot(density(yrep[3,]))
# compare minimum of observed data to minimum from
# replicated samples
# minimum of replicated sample
maxs = apply(yrep, 1, max)
# estimated p-value
(sum(maxs >= max(y + 1))/(length(maxs) + 1))
# histogram comparing T(y) and T(yrep)
ppc_stat(y, yrep, stat = "max")
ggsave("mod2_ppc_stat_max.png")

# minimum of replicated sample
mins = apply(yrep, 1, min)
# estimated p-value
(sum(mins <= min(y + 1))/(length(mins) + 1))
ppc_stat(y, yrep, stat = "min")
ggsave("mod2_ppc_stat_min.png")

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
    mu = fit[i,i_mu[j]]
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
ggsave("mod2_dens_overlay.png")
# compare ecdfs of y and yrep
ppc_ecdf_overlay(y, yrep[sample(nyrep, 20) , ])
ggsave("mod2_ecdf_overlay.png")
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

# Sanity check on lognormal
mean = exp(fit[i,i_mu[1]])
sigmasq = fit[i, i_sigma]
location <- log(mean^2 / sqrt(sigmasq + mean^2))
shape <- sqrt(log(1 + (sigmasq / mean^2)))
set.seed(123)
samples1 = rlnorm(1000, meanlog = location, sdlog = shape)
samples2 = rnorm(1000, mean = mean, sd = sqrt(sigmasq))

# Run the model loop and 
run_fits()

# Compare summaries and information criteria
compares = run_compare()
looic_compare = compares[[1]]
waic_compare = compares[[2]]
xtable(looic_compare)
xtable(waic_compare)

xtable(cbind(looic_compare[,"looic"], waic_compare[,"waic"]))