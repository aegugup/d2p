
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

mod5_vars = c(demos, "NUM_IOI_PER_KM2", "NUM_DISCHARGE", "NUM_FED_SITES", 
              "NUM_BIOSOLIDS_PER_KM2", "AVG_AVG_CONC_mgL")
mod5_x = cbind(x0 = rep(1, each = length(total_pfas_df$CONCENTRATION_ngL)),
               total_pfas_df[,mod5_vars])


all_mod_x = list(mod1_x, mod2_x, mod3_x, mod4_x, mod5_x)


get_fit = function(mod_x){
  
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
  
  return(fit)
}
  
save_fit = function(i_mod, mod_x){
  
  fit = get_fit(mod_x)
  
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

run_fits = function(){
  for (i in 1:length(all_mod_x)){
    save_fit(i, all_mod_x[[i]])
  }
}



# ---------------- Compare LOOIC and WAIC

run_compare = function(){
  load("mod1_looic.rda")
  mod1_looic = mod_looic
  load("mod2_looic.rda")
  mod2_looic = mod_looic
  load("mod3_looic.rda")
  mod3_looic = mod_looic
  load("mod4_looic.rda")
  mod4_looic = mod_looic
  load("mod5_looic.rda")
  mod5_looic = mod_looic
  
  load("mod1_waic.rda")
  mod1_waic = mod_waic
  load("mod2_waic.rda")
  mod2_waic = mod_waic
  load("mod3_waic.rda")
  mod3_waic = mod_waic
  load("mod4_waic.rda")
  mod4_waic = mod_waic
  load("mod5_waic.rda")
  mod5_waic = mod_waic
  
  looic_compare = loo_compare(mod1_looic, 
                              mod2_looic,
                              mod3_looic,
                              mod4_looic,
                              mod5_looic)
  looic_compare
  xtable(looic_compare)
  
  
  waic_compare = loo_compare(mod1_waic, 
                             mod2_waic,
                             mod3_waic,
                             mod4_waic,
                             mod5_looic)
  waic_compare
  xtable(waic_compare)
  
  xtable(cbind(looic_compare[,"looic"], waic_compare[,"waic"]))
  
  return(list(looic_compare, waic_compare))
}
