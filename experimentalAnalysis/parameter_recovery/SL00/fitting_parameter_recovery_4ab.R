###############################################################################################
##
## Parameter recovery check for the task used in 20/08/2020 (two-armed risky bandit)
## Note: This is a code for recovering. So you have to run a simulation before running this code
##
## 豊川 航 / TOYOKAWA, Wataru (11 January 2021, University of Konstanz, Germany)
###############################################################################################
if(FALSE) rstudioapi::restartSession(command = "print('Session Restarted!')") # restart Rstudio session
# use this if you think that r code executed ahead of the following code have 
# used up many memories and/or temporal data
if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace # note, this might cause some cmdstanr error (not sure why tho)

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(cowplot)
library(loo)
# rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores()) # regacy from rstan?

# path -- mac book pro
this_path <- '~/Documents/RiskySocialLearning_R/data_4ab_1022/analysis' # mac book pro
save_path <- '~/Documents/RiskySocialLearning_R/data_4ab_1022/analysis/fitting_results' # mac book pro

# path -- mac mini
this_path <- '~/Documents/riskySocialLearning/data_4ab_1022/analysis' # mac mini
save_path <- '~/Documents/riskySocialLearning/data_4ab_1022/analysis/fitting_results' # mac mini

# path -- drop box
dropbox_path <- "~/Dropbox/wataru/papers/RiskySocialLearning/experiment/parameter_recovery/SL00/"
external_path <- "/Volumes/LaCie/riskySocialLearning/parameter_recovery_4ab"

## Load Functions
source(paste0(dropbox_path , 'functions.R'))

# debug
chains = 4
parallel_chains = 4
thin = 1
iter_warmup = 40
iter_sampling = 80
# debug - 2
chains = 6
parallel_chains = 6
thin = 1
iter_warmup = 1000
iter_sampling = 2000
# # mcmc params
# chains <- 6
# parallel_chains <- 6
# thin <- 2
# iter_warmup <- 2000
# iter_sampling <- 6000

## data
pseudoData_4ab_riskID11_data <- read.csv(paste0(dropbox_path , "pseudoData_4ab_riskID11_data.csv"))
true_global_parameters_4ab_riskID11 <- read.csv(paste0(dropbox_path , "true_global_parameters_4ab_riskID11.csv"))
true_individual_parameters_4ab_riskID11 <- read.csv(paste0(dropbox_path , "true_individual_parameters_4ab_riskID11.csv"))

pseudoData_4ab_riskID12_data <- read.csv(paste0(dropbox_path , "pseudoData_4ab_riskID12_data.csv"))
true_global_parameters_4ab_riskID12 <- read.csv(paste0(dropbox_path , "true_global_parameters_4ab_riskID12.csv"))
true_individual_parameters_4ab_riskID12 <- read.csv(paste0(dropbox_path , "true_individual_parameters_4ab_riskID12.csv"))
pseudoData_4ab_riskID12_data$indiv <- pseudoData_4ab_riskID12_data$indiv + max(pseudoData_4ab_riskID11_data$indiv)
pseudoData_4ab_riskID12_data$group <- pseudoData_4ab_riskID12_data$group + max(pseudoData_4ab_riskID11_data$group)

pseudoData_4ab_data <- pseudoData_4ab_riskID11_data %>% rbind(pseudoData_4ab_riskID12_data)
pseudoData_4ab_data$task <- rep(1:2, each = nrow(pseudoData_4ab_riskID11_data))
true_global_parameters_4ab <- true_global_parameters_4ab_riskID11 %>% rbind(true_global_parameters_4ab_riskID12)
true_global_parameters_4ab$task <- rep(1:2, each= nrow(true_global_parameters_4ab_riskID11))
true_individual_parameters_4ab <- true_individual_parameters_4ab_riskID11 %>% rbind(true_individual_parameters_4ab_riskID12)
true_individual_parameters_4ab$task <- rep(1:2, each= nrow(true_individual_parameters_4ab_riskID11))
true_individual_parameters_4ab$indiv <- 1:nrow(true_individual_parameters_4ab)


# stan model
stan_file_SL00_multiVar_LKJ <- file.path(this_path, 'model_SL00_multiVar_LKJ_multi.stan')
stanmodel_SL00_multiVar_LKJ <- cmdstan_model(stan_file_SL00_multiVar_LKJ)
stanmodel_SL00_multiVar_LKJ$exe_file()
num_parameter <- 4


## ==========  Data cleaning
# preparing data fed to Stan
parameterRecovery_4ab_stan_data = list(
    All = nrow(pseudoData_4ab_data),
    Nsub = length(table(pseudoData_4ab_data$indiv)),
    Ncue = 4, # number of options
    Ntrial = 70,
    Ncon = 2,
    sub = pseudoData_4ab_data$indiv,
    Y = pseudoData_4ab_data$choices,
    trial = pseudoData_4ab_data$trial,
    payoff = pseudoData_4ab_data$payoff,
    condition = true_individual_parameters_4ab$task
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
            if( subset(data,indiv==i&trial==(t-1))$choices>0 ) {
                lastChoice = subset(data,indiv==i&trial==(t-1))$choices
            }
            F[i,1,t] <- data %>% dplyr::filter(indiv==i&trial==t) %>% .$socialFreq_1
            F[i,2,t] <- data %>% dplyr::filter(indiv==i&trial==t) %>% .$socialFreq_2
            F[i,3,t] <- data %>% dplyr::filter(indiv==i&trial==t) %>% .$socialFreq_3
            F[i,4,t] <- data %>% dplyr::filter(indiv==i&trial==t) %>% .$socialFreq_4
            # F[i,3,t] = subset(data,indiv==i&trial==t)$socialFreq2
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
F0 = array(rep(NA,nrow(pseudoData_4ab_data)), c(parameterRecovery_4ab_stan_data$Nsub
  , parameterRecovery_4ab_stan_data$Ncue
  , parameterRecovery_4ab_stan_data$Ntrial))
F = F_calculation_4ab.compiled(F0
  , parameterRecovery_4ab_stan_data$Nsub
  , parameterRecovery_4ab_stan_data$Ntrial
  , pseudoData_4ab_data)

parameterRecovery_4ab_stan_data$F = F
## ==========  Data cleaning END



## =======================================================================
## SL00_multiVar_LKJ model fitting
#  =======================================================================
fit_SL00_multiVar_LKJ_param_recov_4ab = stanmodel_SL00_multiVar_LKJ$sample(
# fit_SL00_multiVar_LKJ_param_recov_4ab = stanmodel_SL00_multiVar_LKJ$variational( # Activate if you want to run VB to quickly check the model

  data = parameterRecovery_4ab_stan_data
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
  #refresh = 200,
  , init = function(chain_id) {
    list(mu_beta = runif(parameterRecovery_4ab_stan_data$Ncon, -2, 2)# because (beta + annealing*t/70) both beta and annealing shouldn't be too large!
        , mu_alpha = runif(parameterRecovery_4ab_stan_data$Ncon, -2, 2)
        , mu_soc0 = runif(parameterRecovery_4ab_stan_data$Ncon, -2, 2)
        , mu_theta = runif(parameterRecovery_4ab_stan_data$Ncon, -1, 3)
        # , mu_soc_slope = runif(parameterRecovery_4ab_stan_data$Ncon, -2, 2)
        # , mu_epsilon = runif(parameterRecovery_4ab_stan_data$Ncon, -2, 2)
        )
  }
  , adapt_delta = 0.9 # default 0.8
  # # == VB method! ==
  # algorithm = "meanfield", #'meanfield' "fullrank"
  # iter = 1000 #(Default)
  # )
  # # == END - VB method! ==

  , chains = chains, parallel_chains = parallel_chains, thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling #debug
  # chains = 4, parallel_chains = 4, thin = 2, iter_warmup = 1000, iter_sampling = 2000
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
fit_SL00_multiVar_LKJ_param_recov_4ab_globalparameters = fit_SL00_multiVar_LKJ_param_recov_4ab$summary(pars %>% append('Rho_ID'))
write.csv(fit_SL00_multiVar_LKJ_param_recov_4ab_globalparameters,
      paste0(dropbox_path , "fit_SL00_multiVar_LKJ_param_recov_4ab_globalparameters.csv"),
      row.names=FALSE)

# Plotting the two fit parameters
fit_SL00_multiVar_LKJ_param_recov_4ab_parameters = data.frame(
  sub = 1:max(parameterRecovery_4ab_stan_data$sub),
  alpha_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('alpha')$mean,
  alpha_q5_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('alpha')$q5,
  alpha_q95_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('alpha')$q95,
  log_beta_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('log_beta')$mean,
  log_beta_q5_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('log_beta')$q5,
  log_beta_q95_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('log_beta')$q95,
  # epsilon_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('epsilon')$mean,
  # epsilon_q5_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('epsilon')$q5,
  # epsilon_q95_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('epsilon')$q95,
  logit_soc0_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('logit_soc0')$mean,
  logit_soc0_q5_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('logit_soc0')$q5,
  logit_soc0_q95_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('logit_soc0')$q95,
  # logit_soc_slope_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('logit_soc_slope')$mean,
  # logit_soc_slope_q5_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('logit_soc_slope')$q5,
  # logit_soc_slope_q95_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('logit_soc_slope')$q95,
  theta_mean_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('theta')$mean,
  theta_q5_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('theta')$q5,
  theta_q95_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('theta')$q95,
  log_lik_SL00_multiVar_LKJ = fit_SL00_multiVar_LKJ_param_recov_4ab$summary('log_lik')$mean,
  waic_SL00_multiVar_LKJ = WAIC_cmdstanr_indv_MCMC(fit_SL00_multiVar_LKJ_param_recov_4ab, nsample=chains*(iter_sampling/thin), nsub=max(parameterRecovery_4ab_stan_data$sub))$waic
  # waic_SL00_multiVar_LKJ = WAIC_cmdstanr_indv_VB(fit_SL00_multiVar_LKJ_param_recov_4ab)$waic
  )

write.csv(fit_SL00_multiVar_LKJ_param_recov_4ab_parameters
	, paste0(dropbox_path , "fit_SL00_multiVar_LKJ_param_recov_4ab_parameters.csv")
	, row.names=FALSE)


# # Changing-social-weight parameters
# fit_SL00_multiVar_LKJ_param_recov_4ab$summary('soc') -> soc_data_SL00_multiVar_LKJ_param_recov_4ab
# soc_data_SL00_multiVar_LKJ_param_recov_4ab$t <- rep(1:70, each = max(parameterRecovery_4ab_stan_data$sub))
# soc_data_SL00_multiVar_LKJ_param_recov_4ab$sub <- rep(1:max(parameterRecovery_4ab_stan_data$sub), 70)
# write.csv(soc_data_SL00_multiVar_LKJ_param_recov_4ab
# 	, paste0(dropbox_path , "soc_data_SL00_multiVar_LKJ_param_recov_4ab.csv")
# 	, row.names=FALSE)

# soc_data_SL00_multiVar_LKJ_param_recov_4ab$task <- rep(1:2, each=nrow(true_individual_parameters_4ab)/2) %>% rep(70)

# ( soc_data_SL00_multiVar_LKJ_param_recov_4ab %>% ggplot() +
#     geom_line(aes(t, mean, group=sub)) +
#     stat_summary(aes(t, mean), fun='median', geom='line', colour='red', size = 1)+
#     labs(title='Estimated dynamics of social weight', x = 'Trials', y = 'Fit social weight')+
#     myTheme_Helvetica()+
#     ylim(c(0,1))+
#     facet_grid(task ~ .) +
#     NULL -> soc_data_SL00_multiVar_LKJ_param_recov_4ab_plot
# )
# pseudoData_4ab_data %>% ggplot() +
#     geom_line(aes(trial, net_sigma, group=indiv)) +
#     stat_summary(aes(trial, net_sigma), fun='median', geom='line', colour='red', size = 1)+
#     labs(title='True dynamics of social weight', x = 'Trials', y = 'True social weight')+
#     myTheme_Helvetica()+
#     ylim(c(0,1))+
#     facet_grid(task ~ .) +
#     NULL -> true_sigma_dynamics

# sigma_true_vs_fit <- plot_grid(true_sigma_dynamics, soc_data_SL00_multiVar_LKJ_param_recov_4ab_plot, labels = c('',''), ncol = 2, align = 'v')
# ggsave(file = paste0(dropbox_path , 'sigma_true_vs_fit.pdf')
# 	, plot = sigma_true_vs_fit
# 	, width=6, height=4.5)


# # Annealing parameters
# fit_SL00_multiVar_LKJ_param_recov_4ab$summary('netBeta') -> netBeta_data_SL00_multiVar_LKJ_param_recov_4ab
# netBeta_data_SL00_multiVar_LKJ_param_recov_4ab$t <- rep(1:70, each = max(parameterRecovery_4ab_stan_data$sub))
# netBeta_data_SL00_multiVar_LKJ_param_recov_4ab$sub <- rep(1:max(parameterRecovery_4ab_stan_data$sub), 70)
# write.csv(netBeta_data_SL00_multiVar_LKJ_param_recov_4ab
# 	, paste0(dropbox_path , "netBeta_data_SL00_multiVar_LKJ_param_recov_4ab.csv")
# 	, row.names=FALSE)

# netBeta_data_SL00_multiVar_LKJ_param_recov_4ab$task <- rep(1:2, each=nrow(true_individual_parameters_4ab)/2) %>% rep(70)

# ( netBeta_data_SL00_multiVar_LKJ_param_recov_4ab %>% ggplot() +
#     geom_line(aes(t, mean, group=sub)) +
#     stat_summary(aes(t, mean), fun='median', geom='line', colour='red', size = 1)+
#     labs(title='Estimated inv. temp.', x = 'Trials', y = 'Fit inverse temp.')+
#     myTheme_Helvetica()+
#     ylim(c(0, 100))+
#     facet_grid(task ~ .) +
#     NULL -> netBeta_data_SL00_multiVar_LKJ_param_recov_4ab_plot
# )
# pseudoData_4ab_data %>% ggplot() +
#     geom_line(aes(trial, net_beta, group=indiv)) +
#     stat_summary(aes(trial, net_beta), fun='median', geom='line', colour='red', size = 1)+
#     labs(title='True inv. temp.', x = 'Trials', y = 'True social weight')+
#     myTheme_Helvetica()+
#     ylim(c(0, 100))+
#     facet_grid(task ~ .) +
#     NULL -> true_beta_dynamics

# # merging figures
# dynamic_parameters_true_vs_fit <- plot_grid(
# 	true_sigma_dynamics
# 	, soc_data_SL00_multiVar_LKJ_param_recov_4ab_plot
# 	, true_beta_dynamics
# 	, netBeta_data_SL00_multiVar_LKJ_param_recov_4ab_plot
# 	, labels = c('','','',''), ncol = 2, align = 'v')

# ggsave(file = paste0(dropbox_path , 'dynamic_parameters_true_vs_fit.pdf')
# 	, plot = dynamic_parameters_true_vs_fit
# 	, width=9, height=9)

gc(); gc()


## =======================================================================
## Individual parameters (showing correlations)
#  =======================================================================

# correlations between indiv parameters
individual_parameters_4ab <- true_individual_parameters_4ab %>% left_join(fit_SL00_multiVar_LKJ_param_recov_4ab_parameters, by = c("indiv" = "sub"))

alpha_plot <- ggplot(individual_parameters_4ab) +
	geom_point(aes(alpha, alpha_mean_SL00_multiVar_LKJ, colour=sqrt((alpha - alpha_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = 'alpha' , x = 'True' , y = 'Fit' ) +
	xlim(c(0,1)) +
	ylim(c(0,1)) +
	myTheme_Helvetica() +
	NULL

log_beta_plot <- ggplot(individual_parameters_4ab) +
	geom_point(aes(log_beta, log_beta_mean_SL00_multiVar_LKJ, colour=sqrt((log_beta - log_beta_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = 'log_beta' , x = 'True' , y = 'Fit' ) +
	myTheme_Helvetica() +
	NULL

# epsilon_plot <- ggplot(individual_parameters_4ab) +
# 	geom_point(aes(epsilon, epsilon_mean_SL00_multiVar_LKJ, colour=sqrt((epsilon - epsilon_mean_SL00_multiVar_LKJ)^2)))+
# 	scale_colour_viridis_c('Difference') +
# 	labs(title = 'epsilon' , x = 'True' , y = 'Fit' ) +
# 	myTheme_Helvetica() +
# 	NULL

soc0_plot <- ggplot(individual_parameters_4ab) +
	geom_point(aes(soc0, logit_soc0_mean_SL00_multiVar_LKJ, colour=sqrt((soc0 - logit_soc0_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = 'soc0' , x = 'True' , y = 'Fit' ) +
	myTheme_Helvetica() +
	NULL

# soc_slope_plot <- ggplot(individual_parameters_4ab) +
# 	geom_point(aes(soc_slope, logit_soc_slope_mean_SL00_multiVar_LKJ, colour=sqrt((soc_slope - logit_soc_slope_mean_SL00_multiVar_LKJ)^2)))+
# 	scale_colour_viridis_c('Difference') +
# 	labs(title = 'soc_slope' , x = 'True' , y = 'Fit' ) +
# 	myTheme_Helvetica() +
# 	NULL

theta_plot <- ggplot(individual_parameters_4ab) +
	geom_point(aes(theta, theta_mean_SL00_multiVar_LKJ, colour=sqrt((theta - theta_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = 'theta' , x = 'True' , y = 'Fit' ) +
	myTheme_Helvetica() +
	NULL

## Correlation coefficients
cor.test(individual_parameters_4ab$alpha, individual_parameters_4ab$alpha_mean_SL00_multiVar_LKJ)
cor.test(individual_parameters_4ab$log_beta, individual_parameters_4ab$log_beta_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_4ab$epsilon, individual_parameters_4ab$epsilon_mean_SL00_multiVar_LKJ)
cor.test(individual_parameters_4ab$soc0, individual_parameters_4ab$logit_soc0_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_4ab$soc_slope, individual_parameters_4ab$logit_soc_slope_mean_SL00_multiVar_LKJ)
cor.test(individual_parameters_4ab$theta, individual_parameters_4ab$theta_mean_SL00_multiVar_LKJ)


# merging figures
indiv_parameters_true_vs_fit <- plot_grid(
	alpha_plot
	, log_beta_plot
	# , epsilon_plot
	, soc0_plot
	# , soc_slope_plot
	, theta_plot
	, labels = c('a','b','c','d','e','f'), ncol = 2, align = 'v')

ggsave(file = paste0(dropbox_path , 'indiv_parameters_true_vs_fit.pdf')
	, plot = indiv_parameters_true_vs_fit
	, width=12, height=6)


## =======================================================================
## Global parameters
#  =======================================================================
true_global_parameters_4ab_long <- true_global_parameters_4ab %>% pivot_longer(-task, names_to = "variable", values_to = "true")

fit_global_parameters_4ab_long <-
  fit_SL00_multiVar_LKJ_param_recov_4ab_globalparameters[1:(2 * num_parameter * 2),] %>%
  separate(variable, into=c('variable', 'task]'), sep='\\[') %>%
  separate('task]', into=c('task', ']'), sep='\\]')

ggplot(data = fit_global_parameters_4ab_long) +
	geom_point(data=true_global_parameters_4ab_long, aes(variable, true), size=2, colour='red')+
	geom_errorbar(aes(variable, ymin=q5, ymax=q95)) +
	geom_point(aes(variable, mean)) +
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle=90))+
	labs(title='Global parameters' , x = 'Parameter name' , y = 'Values')+
  facet_grid(task ~ .) +
	NULL -> global_recovery

ggsave(file = paste0(dropbox_path , 'global_recovery.pdf')
	, plot = global_recovery
	, width=9, height=6)


## =======================================================================
## Individual parameter recovery
#  =======================================================================
## Alpha
subject_order_alpha = order(individual_parameters_4ab$alpha)
alpha_recovery_data = individual_parameters_4ab[subject_order_alpha,]
alpha_recovery_data$subject = 1:nrow(alpha_recovery_data)
alpha_recovery <- ggplot(data=alpha_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=alpha_q5_SL00_multiVar_LKJ, ymax=alpha_q95_SL00_multiVar_LKJ), width=.001) +
	geom_point(aes(subject, alpha_mean_SL00_multiVar_LKJ), colour='blue',size=.5)+
	geom_point(aes(subject, alpha), colour='red',size=.5) +
	myTheme_Helvetica()+
	ylim(c(0,1)) +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',alpha,sep="")))

## log_beta
subject_order_log_beta = order(individual_parameters_4ab$log_beta)
log_beta_recovery_data = individual_parameters_4ab[subject_order_log_beta,]
log_beta_recovery_data$subject = 1:nrow(log_beta_recovery_data)
log_beta_recovery <- ggplot(data=log_beta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=log_beta_q5_SL00_multiVar_LKJ, ymax=log_beta_q95_SL00_multiVar_LKJ), width=.001) +
	geom_point(aes(subject, log_beta_mean_SL00_multiVar_LKJ), colour='blue',size=.5)+
	geom_point(aes(subject, log_beta), colour='red',size=.5) +
	myTheme_Helvetica()+
	#ylim(c(0,1)) +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y='Individual log_beta')

## epsilon
# subject_order_epsilon = order(individual_parameters_4ab$epsilon)
# epsilon_recovery_data = individual_parameters_4ab[subject_order_epsilon,]
# epsilon_recovery_data$subject = 1:nrow(epsilon_recovery_data)
# epsilon_recovery <- ggplot(data=epsilon_recovery_data) +
# 	geom_errorbar(aes(x=subject, ymin=epsilon_q5_SL00_multiVar_LKJ, ymax=epsilon_q95_SL00_multiVar_LKJ), width=.001) +
# 	geom_point(aes(subject, epsilon_mean_SL00_multiVar_LKJ), colour='blue',size=.5)+
# 	geom_point(aes(subject, epsilon), colour='red',size=.5) +
# 	myTheme_Helvetica()+
# 	#ylim(c(0,1)) +
# 	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
# 	labs(x='Agents', y='Individual epsilon')

## soc0
subject_order_soc0 = order(individual_parameters_4ab$soc0)
soc0_recovery_data = individual_parameters_4ab[subject_order_soc0,]
soc0_recovery_data$subject = 1:nrow(soc0_recovery_data)
soc0_recovery <- ggplot(data=soc0_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=logit_soc0_q5_SL00_multiVar_LKJ, ymax=logit_soc0_q95_SL00_multiVar_LKJ), width=.001) +
	geom_point(aes(subject, logit_soc0_mean_SL00_multiVar_LKJ), colour='blue',size=.5)+
	geom_point(aes(subject, soc0), colour='red',size=.5) +
	myTheme_Helvetica()+
	#ylim(c(0,1)) +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y='Individual soc0')

## soc_slope
# subject_order_soc_slope = order(individual_parameters_4ab$soc_slope)
# soc_slope_recovery_data = individual_parameters_4ab[subject_order_soc_slope,]
# soc_slope_recovery_data$subject = 1:nrow(soc_slope_recovery_data)
# soc_slope_recovery <- ggplot(data=soc_slope_recovery_data) +
# 	geom_errorbar(aes(x=subject, ymin=logit_soc_slope_q5_SL00_multiVar_LKJ, ymax=logit_soc_slope_q95_SL00_multiVar_LKJ), width=.001) +
# 	geom_point(aes(subject, logit_soc_slope_mean_SL00_multiVar_LKJ), colour='blue',size=.5)+
# 	geom_point(aes(subject, soc_slope), colour='red',size=.5) +
# 	myTheme_Helvetica()+
# 	#ylim(c(0,1)) +
# 	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
# 	labs(x='Agents', y='Individual soc_slope')

## theta
subject_order_theta = order(individual_parameters_4ab$theta)
theta_recovery_data = individual_parameters_4ab[subject_order_theta,]
theta_recovery_data$subject = 1:nrow(theta_recovery_data)
theta_recovery <- ggplot(data=theta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=theta_q5_SL00_multiVar_LKJ, ymax=theta_q95_SL00_multiVar_LKJ), width=.001) +
	geom_point(aes(subject, theta_mean_SL00_multiVar_LKJ), colour='blue',size=.5)+
	geom_point(aes(subject, theta), colour='red',size=.5) +
	myTheme_Helvetica()+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',theta,sep="")))


# merging figures
indiv_parameters_recovery <- plot_grid(
	alpha_recovery
	, log_beta_recovery
	# , epsilon_recovery
	, soc0_recovery
	# , soc_slope_recovery
	, theta_recovery
	, labels = c('a','b','c','d','e','f'), ncol = 2, align = 'v')

ggsave(file = paste0(dropbox_path , 'indiv_parameters_recovery.pdf')
	, plot = indiv_parameters_recovery
	, width=12, height=6)



# diagnosis
# Computing PSIS-LOO and checking diagnostics
(loo_SL00_multiVar_LKJ_param_recov_4ab <- loo(fit_SL00_multiVar_LKJ_param_recov_4ab$draws('log_lik'), save_psis = TRUE))

# Plotting posterior distributions
bayesplot::mcmc_hist(fit_SL00_multiVar_LKJ_param_recov_4ab$draws(pars))

# Extract posterior draws for later use
# np_SL00_multiVar_LKJ_param_recov_4ab <- nuts_params(fit_SL00_multiVar_LKJ_param_recov_4ab$draws(pars))
# posterior_SL00_multiVar_LKJ_param_recov_4ab <- as.array(fit_SL00_multiVar_LKJ_param_recov_4ab$draws(pars))

# trace plot
color_scheme_set("mix-brightblue-gray")
# bayesplot::mcmc_trace(fit_SL00_multiVar_LKJ_param_recov_4ab$draws(), pars = pars)
mcmc_trace_SL00_multiVar_LKJ_param_recov_4ab <- bayesplot::mcmc_trace(fit_SL00_multiVar_LKJ_param_recov_4ab$draws(pars)
	# , pars = c('mu_alpha', 'mu_beta', 'mu_soc0', 'mu_theta')
	# , np = np_SL00_multiVar_LKJ_param_recov_4ab
	)
ggsave(file = paste0(dropbox_path , 'mcmc_trace_SL00_multiVar_LKJ_param_recov_4ab.pdf')
  , plot = mcmc_trace_SL00_multiVar_LKJ_param_recov_4ab
  , width=12, height=12)

# diagnosis
# diagnostics_df <- as_draws_df(fit_SL00_multiVar_LKJ_param_recov_4ab$sampler_diagnostics())
# print(diagnostics_df)
# fit_SL00_multiVar_LKJ_param_recov_4ab$cmdstan_diagnose()



# save
fit_SL00_multiVar_LKJ_param_recov_4ab$save_object(file = paste0(external_path, 'fit_SL00_multiVar_LKJ_param_recov_4ab.RDS'))
# can be read back in using readRDS
# fit_SL00_multiVar_LKJ_param_recov_4ab <- readRDS(paste0(external_path, 'fit_SL00_multiVar_LKJ_param_recov_4ab.RDS')) 

gc(); gc()
