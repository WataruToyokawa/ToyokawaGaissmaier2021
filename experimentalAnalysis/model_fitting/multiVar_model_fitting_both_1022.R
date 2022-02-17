###############################################################################################
##
## Analyses of the behavioural patterns
## - (1) Were group members' behaviours correlated?
##
## 豊川 航 / TOYOKAWA, Wataru (11 August 2020, University of Konstanz, Germany)
###############################################################################################
if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace # note, this might cause some cmdstanr error (not sure why tho)

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
#rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores()) # regacy from rstan?

## Load Functions
source('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/functions.R')
## Load Functions
source('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/parameter_recovery/SL00/functions.R')

# path -- mac mini
this_path <- '~/Documents/riskySocialLearning/data_4ab_1022/analysis' # mac mini
save_path <- '~/Documents/riskySocialLearning/data_4ab_1022/analysis/fitting_results' # mac mini

# path -- mac book air
# this_path <- '~/Documents/riskySocialLearning_mba/data_4ab_1022/analysis' # mac book air
# save_path <- '~/Documents/riskySocialLearning/data_4ab_1022/analysis/fitting_results' # mac book air

# # path -- mac book pro
# this_path <- '~/Documents/RiskySocialLearning_R/data_4ab_1022/analysis' # mac book pro
# save_path <- '~/Documents/RiskySocialLearning_R/data_4ab_1022/analysis/fitting_results' # mac book pro

# path -- drop box
dropbox_path <- "~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/fitting_results/"

## =======================================================================
## MCMC setting
#  =======================================================================
# debug
chains <- 2
parallel_chains <- 2
thin <- 1
iter_warmup <- 20
iter_sampling <- 120
#
chains = 6
parallel_chains = 6
thin = 2
iter_warmup = 1500
iter_sampling = 3000

# parameter estimation
chains <- 6
parallel_chains <- 6
thin <- 4
iter_warmup <- 2000
iter_sampling <- 6000

# behaivoural data
allBehaviour1022_group <- read.csv("~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/results/allBehaviour1022_group.csv")
allBehaviour1022_group$sub <- allBehaviour1022_group$amazonID %>% as.factor() %>% as.numeric()
allBehaviour1022_group$task <- 1
allBehaviour1022_group$task[which(allBehaviour1022_group$riskDistributionId==12)] <- 2

summary_each_individual <- allBehaviour1022_group %>%
	group_by(sub) %>%
	summarise(mean_risky_choice_rate = mean(risky_choice), na.rm=TRUE
		, mean_payoff = mean(payoff/100)
		, task = task[1]
		)

# preparing data fed to Stan
ReinforcementLearningStanData_group = list(
    All = nrow(allBehaviour1022_group),
    Nsub = length(table(allBehaviour1022_group$amazonID)),
    Ncue = 4, # number of options
    Ntrial = 70,
    Ncon = 2,
    sub = allBehaviour1022_group$sub,
    Y = allBehaviour1022_group$choice_num + 1,
    trial = allBehaviour1022_group$round,
    payoff = allBehaviour1022_group$payoff / 100,
    condition = summary_each_individual$task
    )

## compfun() で高速化
library(compiler)
F_calculation_4ab = function(array, Nsub, Ntrial, data)
{
    F <- array
    # F (Simulation data とは t+1 ずれてるから注意！)
    for(i in 1:Nsub) {
        F[i,1,1] = 0; F[i,2,1] = 0; F[i,3,1] = 0; F[i,4,1] = 0;
        for(t in 2:Ntrial) {
            lastChoice = 0
            if( subset(data,sub==i&round==(t-1))$choice_num+1>0 ) {
                lastChoice = subset(data,sub==i&round==(t-1))$choice_num + 1
            }
            F[i,1,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_safe1
            F[i,2,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_safe2
            F[i,3,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_safe3
            F[i,4,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_risky
            # F[i,3,t] = subset(data,sub==i&round==t)$socialFreq2
            if(lastChoice>0){
                F[i,lastChoice,t] = F[i,lastChoice,t] - 1
            }
            if(length(which(F[i,,t]>0))==0){
                F[i,,t] <- c(-1,-1,-1,-1)
            }
        }
    }
    return(F)
}

F_calculation_4ab.compiled <- cmpfun(F_calculation_4ab)
F0 = array(rep(NA,nrow(allBehaviour1022_group)), c(ReinforcementLearningStanData_group$Nsub, ReinforcementLearningStanData_group$Ncue, ReinforcementLearningStanData_group$Ntrial))
F = F_calculation_4ab.compiled(F0, ReinforcementLearningStanData_group$Nsub, ReinforcementLearningStanData_group$Ntrial, allBehaviour1022_group)

ReinforcementLearningStanData_group$F = F
## ==========  Data cleaning END


# stan model
stan_file_SL00_multiVar_LKJ <- file.path(this_path, 'model_SL00_multiVar_LKJ_multi.stan')
stanmodel_SL00_multiVar_LKJ <- cmdstan_model(stan_file_SL00_multiVar_LKJ)
stanmodel_SL00_multiVar_LKJ$exe_file()
num_parameter <- 4

## =======================================================================
## SL00_multiVar_LKJ model fitting
#  =======================================================================
fit_SL00_multiVar_LKJ_1022 = stanmodel_SL00_multiVar_LKJ$sample(
# fit_SL00_multiVar_LKJ_1022 = stanmodel_SL00_multiVar_LKJ$variational( # Activate if you want to run VB to quickly check the model

  data = ReinforcementLearningStanData_group
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
  #refresh = 200,
  , init = function(chain_id) {
    list(mu_beta = runif(ReinforcementLearningStanData_group$Ncon, -2, 2)
        , mu_alpha = runif(ReinforcementLearningStanData_group$Ncon, -2, 2)
        , mu_soc0 = runif(ReinforcementLearningStanData_group$Ncon, -3, 1)
        , mu_theta = runif(ReinforcementLearningStanData_group$Ncon, -1, 3)
        # , mu_soc_slope = runif(ReinforcementLearningStanData_group$Ncon, -2, 2)
        # , mu_epsilon = runif(ReinforcementLearningStanData_group$Ncon, -2, 2)
        )
  }
  , adapt_delta = 0.9 # default 0.8

  , chains = chains, parallel_chains = parallel_chains, thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 
)


pars <- c('mu_alpha'
    ,'mu_beta'
    # ,'mu_epsilon'
    ,'mu_soc0'
    # ,'mu_soc_slope'
    ,'mu_theta'
    ,'s_alpha'
    ,'s_beta'
    # ,'s_epsilon'
    ,'s_soc0'
    # ,'s_soc_slope'
    ,'s_theta'
    )

# save global parameters
fit_SL00_multiVar_LKJ_1022_globalparameters = fit_SL00_multiVar_LKJ_1022$summary(c(pars,'Rho_ID'))
write.csv(fit_SL00_multiVar_LKJ_1022_globalparameters,
      paste0(dropbox_path, "fit_SL00_multiVar_LKJ_1022_globalparameters.csv"),
      row.names=FALSE)

# Plotting the two fit parameters
fit_SL00_multiVar_LKJ_1022_parameters = data.frame(
  sub = 1:max(ReinforcementLearningStanData_group$sub),
  alpha_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('alpha')$mean,
  alpha_median_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('alpha')$median,
  beta_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('beta')$mean,
  beta_median_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('beta')$median,
  # epsilon_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('epsilon')$mean,
  # epsilon_median_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('epsilon')$median,
  soc_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('soc')$mean,
  soc_median_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('soc')$median,
  # logit_soc_slope_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('logit_soc_slope')$mean,
  # logit_soc_slope_median_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('logit_soc_slope')$median,
  theta_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('theta')$mean,
  theta_median_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('theta')$median,

  log_lik_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_1022$summary('log_lik')$mean,
  waic_SL00_multiVar_LKJ = WAIC_cmdstanr_indv_MCMC(fit_SL00_multiVar_LKJ_1022, nsample=chains*(iter_sampling/thin), nsub=max(ReinforcementLearningStanData_group$sub))$waic
  # waic_SL00_multiVar_LKJ = WAIC_cmdstanr_indv_VB(fit_SL00_multiVar_LKJ_1022)$waic
  )

# Computing PSIS-LOO and checking diagnostics
(loo_SL00_multiVar_LKJ_1022 <- loo(fit_SL00_multiVar_LKJ_1022$draws('log_lik'), save_psis = TRUE))
#loo_AL_0820 <- loo(fit_asocial_learning_mcmc$draws('log_lik'), save_psis = TRUE)

# save fit object
# fit_SL00_multiVar_LKJ_1022$save_object(file = "~/Documents/riskySocialLearning/data_4ab_1022/analysis/fitting_results/fit_SL00_multiVar_LKJ_1022.RDS")
write.csv(fit_SL00_multiVar_LKJ_1022_parameters, paste0(dropbox_path, "fit_SL00_multiVar_LKJ_1022_parameters.csv"), row.names=FALSE)


gc(); gc()











