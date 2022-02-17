###############################################################################################
##
## Model fitting and Bayesian model comparison. I considered the following models
## (1) Asocial RL
## (2) Decision biasing social influence 
## (3) Value shaping social influence
##
## Author: 豊川 航 / TOYOKAWA, Wataru (11 January 2022, University of Konstanz, Germany)
###############################################################################################
if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace 
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(LaplacesDemon)
library(cowplot)
#rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores()) # regacy from rstan?

## Load Functions
# source('functions.R')
source('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/functions.R')
source('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/parameter_recovery/SL00/functions.R')

# path -- drop box
dropbox_path <- "~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/model_fitting/fitting_results/"

# path -- macbook air M1
this_path <- '~/Documents/riskySocialLearning_airM1' # mac mini
save_path <- '~/Documents/riskySocialLearning_airM1/fitting_results' # mac mini

# # path -- mac mini
# this_path <- '~/Documents/riskySocialLearning/data_4ab_1022/analysis' # mac mini
# save_path <- '~/Documents/riskySocialLearning/data_4ab_1022/analysis/fitting_results' # mac mini

# # path -- mac book pro
# this_path <- '~/Documents/RiskySocialLearning_R/data_4ab_1022/analysis' # mac book pro
# save_path <- '~/Documents/RiskySocialLearning_R/data_4ab_1022/analysis/fitting_results' # mac book pro

################################
# Stan models
################################
# Decision biasing social influence model
stan_file_SL00_multiVar_LKJ <- file.path(this_path, 'model_SL00_multiVar_LKJ.stan')
stanmodel_SL00_multiVar_LKJ <- cmdstan_model(stan_file_SL00_multiVar_LKJ)
stanmodel_SL00_multiVar_LKJ$exe_file()
# Asocial RL model
stan_file_AL00_multiVar_LKJ <- file.path(this_path, 'model_AL00_multiVar_LKJ.stan')
stanmodel_AL00_multiVar_LKJ <- cmdstan_model(stan_file_AL00_multiVar_LKJ)
stanmodel_AL00_multiVar_LKJ$exe_file()
# Value shaping social influence model
stan_file_VS00_multiVar_LKJ <- file.path(this_path, 'model_VS00_multiVar_LKJ.stan')
stanmodel_VS00_multiVar_LKJ <- cmdstan_model(stan_file_VS00_multiVar_LKJ)
stanmodel_VS00_multiVar_LKJ$exe_file()

# debug
chains = 2
parallel_chains = 2
thin = 1
iter_warmup = 50
iter_sampling = 150
# # 
chains = 6
parallel_chains = 6
thin = 1
iter_warmup = 1000
iter_sampling = 2000

################################# 
# The 2 armed bandit task
#################################
# Data clearning
behaviour_main_0820 <-
  read.csv("~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_0820/behaviour_main_0820.csv")
behaviour_indiv_0820 <-
  read.csv("~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_0820/behaviour_indiv_0820.csv")
behaviour_main_0820 = rbind(behaviour_main_0820, behaviour_indiv_0820)
allBehaviour0820 = behaviour_main_0820
allBehaviour0820 = allBehaviour0820 %>% dplyr::filter(amazonID != 'INHOUSETEST2')
# make the choice data binary
allBehaviour0820$choice_num = NA
allBehaviour0820$choice_num[which(allBehaviour0820$choice=='sure')] = 0
allBehaviour0820$choice_num[which(allBehaviour0820$choice=='risky')] = 1
allBehaviour0820$choice_num[which(allBehaviour0820$choice=='miss')] = -1

allBehaviour0820$riskDistributionId_factor = 'Condition 5'
allBehaviour0820$riskDistributionId_factor[which(allBehaviour0820$riskDistributionId==6)] = 'Condition 6'
allBehaviour0820$riskDistributionId_factor[which(allBehaviour0820$riskDistributionId==7)] = 'Condition 7'

allBehaviour0820$indivOrGroup_factor = 'Individual'
allBehaviour0820$indivOrGroup_factor[allBehaviour0820$indivOrGroup == 1] = 'Social'

allBehaviour0820$groupSize_category = 'Small'
allBehaviour0820$groupSize_category[which(allBehaviour0820$groupSize==1)] = 'Individual'
allBehaviour0820$groupSize_category[which(allBehaviour0820$groupSize>4)] = 'Large'
allBehaviour0820$groupSize_category = factor(allBehaviour0820$groupSize_category, levels = c('Individual','Small','Large'))

## ==========  Data reading and cleaning
# data for stan fitting - social learning
completedIDs = which(table(allBehaviour0820$amazonID) >= 36) %>% names()
allBehaviour0820_social = allBehaviour0820 %>% dplyr::filter(amazonID %in% completedIDs) %>%
  dplyr::filter(indivOrGroup == 1) # note this is only the social condition
# There was a bug in the data storage (might be due to the server crash??)
# So I needed to eliminate this participant:
allBehaviour0820_social = allBehaviour0820_social %>% dplyr::filter(amazonID != '5eac70db94edd22d57fa00c4')
allBehaviour0820_social$sub = as.numeric(as.factor(allBehaviour0820_social$amazonID))
allBehaviour0820_social$group = as.numeric(as.factor(allBehaviour0820_social$room))

# insert missed trials
for(i in 1:length(table(allBehaviour0820_social$amazonID))) {
    thisSubject <- allBehaviour0820_social %>% dplyr::filter(sub==i)
    for(t in 1:70) {
        if(nrow(thisSubject[which(thisSubject$round==t),])==0) {
            newRow = rep(NA, ncol(allBehaviour0820_social))
            names(newRow) <- names(allBehaviour0820_social)
            allBehaviour0820_social <- bind_rows(allBehaviour0820_social, newRow)
            allBehaviour0820_social$choice_num[nrow(allBehaviour0820_social)] = -1
            allBehaviour0820_social$payoff[nrow(allBehaviour0820_social)] = -1
            allBehaviour0820_social$round[nrow(allBehaviour0820_social)] = t
            allBehaviour0820_social$socialFreq_safe[nrow(allBehaviour0820_social)] = 0
            allBehaviour0820_social$socialFreq_risky[nrow(allBehaviour0820_social)] = 0
            allBehaviour0820_social$amazonID[nrow(allBehaviour0820_social)] = allBehaviour0820_social$amazonID[which(allBehaviour0820_social$sub==i)][1]
            allBehaviour0820_social$sub[nrow(allBehaviour0820_social)] = allBehaviour0820_social$sub[which(allBehaviour0820_social$sub==i)][1]
            allBehaviour0820_social$group[nrow(allBehaviour0820_social)] = allBehaviour0820_social$group[which(allBehaviour0820_social$sub==i)][1]
            #allBehaviour0820_social$taskDifficulty_num[nrow(allBehaviour0820_social)] = allBehaviour0820_social$taskDifficulty_num[which(allBehaviour0820_social$sub==i)][1]
        }
    }
}
allBehaviour0820_social = allBehaviour0820_social[order(allBehaviour0820_social$round),] # Sorting by round
allBehaviour0820_social = allBehaviour0820_social %>% group_by(amazonID) %>% arrange(round, .by_group = TRUE)

# preparing data fed to Stan
ReinforcementLearningStanData_group_0820 = list(
    All = nrow(allBehaviour0820_social),
    Nsub = length(table(allBehaviour0820_social$amazonID)),
    Ncue = 2, # number of options
    Ntrial = 70,
    Ncon = 1,
    sub = allBehaviour0820_social$sub,
    Y = allBehaviour0820_social$choice_num + 1,
    trial = allBehaviour0820_social$round,
    payoff = allBehaviour0820_social$payoff / 100
    )

## compfun() で高速化
library(compiler)
F_calculation = function(array, Nsub, Ntrial, data)
{
    F <- array
    # F (Simulation data とは t+1 ずれてるから注意！)
    for(i in 1:Nsub) {
        F[i,1,1] = 0; F[i,2,1] = 0; # F[i,3,1] = 0;
        for(t in 2:Ntrial) {
            lastChoice = 0
            if( subset(data,sub==i&round==(t-1))$choice_num+1>0 ) {
                lastChoice = subset(data,sub==i&round==(t-1))$choice_num + 1
            }
            F[i,1,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_safe
            F[i,2,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_risky
            #F[i,3,t] = subset(data,sub==i&round==t)$socialFreq2
            if(lastChoice>0){
                F[i,lastChoice,t] = F[i,lastChoice,t] - 1
            }
            # if(length(which(F[i,,t]>0))==0){
            #     F[i,,t] <- c(0,0)
            # }
        }
    }
    return(F)
}

F_calculation.compiled <- cmpfun(F_calculation)
F0 = array(rep(NA,nrow(allBehaviour0820_social)), c(ReinforcementLearningStanData_group_0820$Nsub, ReinforcementLearningStanData_group_0820$Ncue, ReinforcementLearningStanData_group_0820$Ntrial))
F = F_calculation.compiled(F0, ReinforcementLearningStanData_group_0820$Nsub, ReinforcementLearningStanData_group_0820$Ntrial,allBehaviour0820_social)

ReinforcementLearningStanData_group_0820$F = F
## ==========  Data cleaning END

# -- Model fitting--

# -- Value shaping model --
pars <- c('mu_alpha'
    ,'mu_beta'
    ,'mu_soc_vs'
    ,'mu_theta'
    ,'s_alpha'
    ,'s_beta'
    ,'s_soc_vs'
    ,'s_theta'
    )
pars_AL <- c('mu_alpha','mu_beta','s_alpha','s_beta')
# Variational Bayesian fit (for init values used in the MCMC blow)
# fit_VS00_multiVar_LKJ_0820_vb = stanmodel_VS00_multiVar_LKJ$variational(

#   data = ReinforcementLearningStanData_group_0820,
#   # == VB method! ==
#   algorithm = "meanfield", #'meanfield' "fullrank"
#   iter = 1000 #(Default)
#   # == END - VB method! ==
# )
# fit_VS00_multiVar_LKJ_0820_globalparameters_vb = fit_VS00_multiVar_LKJ_0820_vb$summary(pars %>% append(c('Rho_ID')))
# write.csv(fit_VS00_multiVar_LKJ_0820_globalparameters_vb
#   , paste0(save_path , "/fit_VS00_multiVar_LKJ_0820_globalparameters.csv"), row.names=FALSE)


fit_VS00_multiVar_LKJ_0820 = stanmodel_VS00_multiVar_LKJ$sample(

  data = ReinforcementLearningStanData_group_0820
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
  , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_VS00_multiVar_LKJ_0820_globalparameters = fit_VS00_multiVar_LKJ_0820$summary(pars %>% append(c('Rho_ID')))
write.csv(fit_VS00_multiVar_LKJ_0820_globalparameters
  , paste0(dropbox_path , "fit_VS00_multiVar_LKJ_0820_globalparameters.csv"), row.names=FALSE)

# Plotting the two fit parameters
fit_VS00_multiVar_LKJ_0820_parameters = data.frame(
  sub = 1:max(ReinforcementLearningStanData_group_0820$sub),
  alpha_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('alpha')$mean,
  alpha_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('alpha')$median,
  beta_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('beta')$mean,
  beta_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('beta')$median,
  soc_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('soc_vs')$mean,
  soc_vs_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('soc_vs')$median,
  theta_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('theta')$mean,
  theta_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('theta')$median,

  log_lik_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_0820$summary('log_lik')$mean,
  waic_VS00_multiVar_LKJ = WAIC_cmdstanr_indv_MCMC(fit_VS00_multiVar_LKJ_0820, nsample=chains*(iter_sampling/thin), nsub=max(ReinforcementLearningStanData_group_0820$sub))$waic
  # waic_VS00_multiVar_LKJ = WAIC_cmdstanr_indv_VB(fit_VS00_multiVar_LKJ_0820)$waic
  )

write.csv(fit_VS00_multiVar_LKJ_0820_parameters, paste0(dropbox_path, "fit_VS00_multiVar_LKJ_0820_parameters.csv"), row.names=FALSE)

# -- loading the fitting results of other two models --
fit_AL00_multiVar_LKJ_group_0820_parameters = read_csv('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_0820/analysis/fitting_results/fit_AL00_multiVar_LKJ_group_0820_parameters.csv')
fit_SL00_multiVar_LKJ_0820_parameters = read_csv('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_0820/analysis/fitting_results/fit_SL00_multiVar_LKJ_0820_parameters.csv')
fit_VS00_multiVar_LKJ_0820_parameters = read_csv(paste0(dropbox_path, "fit_VS00_multiVar_LKJ_0820_parameters.csv"))
###############################################################
# ---- Bayesian Model Selection -----
# Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ. 
#       Bayesian model selection for group studies. NeuroImage. 2009;46(4):1004–17.
###############################################################
## STEP 1: Using WAIC as an approximation of model log-evidence
#-----------------------------------
WAICtable <- data.frame(
	subject=as.character(fit_VS00_multiVar_LKJ_0820_parameters$sub),
	WAIC_AL00=fit_AL00_multiVar_LKJ_group_0820_parameters$waic_AL00_multiVar_LKJ,
	WAIC_VS00=fit_VS00_multiVar_LKJ_0820_parameters$waic_VS00_multiVar_LKJ,
	WAIC_SL00=fit_SL00_multiVar_LKJ_0820_parameters$waic_SL00_multiVar_LKJ
	)
# STEP 2. Running the algorithm which updates estimates of the Dirichlet density parameter, α
#         iteretively until convergence
#------------------------------------------------------------------------------------------
num_model = 3
D_alpha0 = rep(1, num_model)
D_alpha = D_alpha0
for(i in 1:500){
	u_WAIC_1 = exp(-WAICtable$WAIC_AL00 + digamma(D_alpha[1]) - digamma(sum(D_alpha)))
	u_WAIC_2 = exp(-WAICtable$WAIC_VS00 + digamma(D_alpha[2]) - digamma(sum(D_alpha)))
	u_WAIC_3 = exp(-WAICtable$WAIC_SL00 + digamma(D_alpha[3]) - digamma(sum(D_alpha)))
	D_beta = c(sum(u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    ,sum(u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    ,sum(u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    )
	D_alpha <- D_alpha0 + D_beta
	r_WAIC = c(D_alpha[1]/sum(D_alpha)
    ,D_alpha[2]/sum(D_alpha)
    ,D_alpha[3]/sum(D_alpha)
    )
	if(i%%10==0){print(D_alpha)}
	#if(i%%100000==0){print(round(r,digits=10))}
	#if(i%%2000==0){print(rbind(D_alpha,r))}
}
# STEP 3: Plot (Posterior model probability)
# the expected multinomial parameters r and thus the expected likelihood of obtaining the k-th model, 
# i.e. p(mnk=1|r)=Mult(m; 1,r), for any randomly selected subject
# -------------------------------------------------------------------------------------------
rWAIC_for_fig = data.frame(
	models = c('Asocial','Value shaping','Decision biasing'),
	r_WAIC = r_WAIC)
rWAIC_for_fig$models = factor(c('Asocial','Value shaping','Decision biasing'),levels=c('Asocial','Value shaping','Decision biasing'))

# STEP 4: Plot (Exceedance Probability)
# (1) sample: (r1, r2, r3) ~ Dirichlet(alpha_1, alpha_2, alpha_3); 
# (2) Count cases where r3 > r1 and r3 > r2
# -------------------------------------------------------------------------------------------
total_sample_N = 10000
dir_sample = rdirichlet(total_sample_N, D_alpha)
xp_model3 = (which(apply(dir_sample, MARGIN=1, FUN=which.max) == 3) %>% length()) / total_sample_N
xp_model3

(r_plot_0820 <- ggplot(rWAIC_for_fig, aes(models,r_WAIC)) + geom_bar(stat="identity")+
	labs(x='Model', y='Posterior probability', title='The 1-risky-1-safe task') +
    annotate(geom="text", x=3, y=0.85, label=paste('XP = ', xp_model3), size=6) + 
	coord_flip() +
	ylim(c(0,1))+
    myTheme_Arial()+
	theme(
		plot.title = element_text(size=16, family='Arial', face='bold'),
	)
)

ggsave(file = paste0(dropbox_path, "r_plot_0820.png"), plot = r_plot_0820, width=6, height=4.5)


# STEP 5. 
# Individual level plot
#------------------------------------------------------------------------------------------
g_forBarplot = data.frame(
	subject=as.numeric(rep(WAICtable$sub,3)),
	model=c(rep('Asocial',nrow(WAICtable)),rep('Value shaping',nrow(WAICtable)),rep('Decision biasing',nrow(WAICtable))),
	g_WAIC=c(u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
        , u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
        , u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
	)
WAICtable$g_1 = u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable$g_2 = u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable$g_3 = u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable$bestModel = matrix(c(WAICtable$g_1,WAICtable$g_2,WAICtable$g_3), ncol = 3) %>% apply(MARGIN = 1, FUN = which.max)
WAICtable$bestModel[which(WAICtable$bestModel==1)] <- 'Asocial'
WAICtable$bestModel[which(WAICtable$bestModel==2)] <- 'Value shaping'
WAICtable$bestModel[which(WAICtable$bestModel==3)] <- 'Decision biasing'
g_forBarplot$bestModel = rep(WAICtable$bestModel, 3)
g_forBarplot$model = factor(g_forBarplot$model,levels=c('Asocial','Value shaping','Decision biasing'))

(individual_g_plot_2ab <- ggplot(data=g_forBarplot, aes(model, g_WAIC, fill = bestModel)) + geom_bar(stat="identity") +
	coord_flip() + facet_wrap(~ subject,ncol = 12) +
	scale_y_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1")) +
    scale_fill_viridis_d(name="Most likely\nlearning model",
						breaks=c("Asocial", "Value shaping", "Decision biasing"),
                        labels=c("AL", "VS", "DB"))+
	labs(x='', y='', title='Posterior belief of models for each individual (the 1-risky-1-safe task)') +
	myTheme_Times()+
    theme(panel.border = element_rect(colour = "black", fill=NA))
)

ggsave(file = paste0(dropbox_path, "individual_g_plot_2ab.png"), plot = individual_g_plot_2ab, width=12, height=12)








################################# 
# The 4 armed bandit task
#################################
pars <- c('mu_alpha'
    ,'mu_beta'
    ,'mu_soc_vs'
    ,'mu_theta'
    ,'s_alpha'
    ,'s_beta'
    ,'s_soc_vs'
    ,'s_theta'
    )
pars_AL <- c('mu_alpha','mu_beta','s_alpha','s_beta')
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
        , amazonID = amazonID[1]
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


# stan model (this model allows for multiple 4ab tasks simultaneously)
stan_file_VS00_multiVar_LKJ_multi <- file.path(this_path, 'model_VS00_multiVar_LKJ_multi.stan')
stanmodel_VS00_multiVar_LKJ_multi <- cmdstan_model(stan_file_VS00_multiVar_LKJ_multi)
stanmodel_VS00_multiVar_LKJ_multi$exe_file()

stan_file_AL00_multiVar_LKJ_multi <- file.path(this_path, 'model_AL00_multiVar_LKJ_multi.stan')
stanmodel_AL00_multiVar_LKJ_multi <- cmdstan_model(stan_file_AL00_multiVar_LKJ_multi)
stanmodel_AL00_multiVar_LKJ_multi$exe_file()

num_parameter <- 4

## =======================================================================
## VS00_multiVar_LKJ model fitting
#  =======================================================================
fit_VS00_multiVar_LKJ_1022 = stanmodel_VS00_multiVar_LKJ_multi$sample(

  data = ReinforcementLearningStanData_group
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
  , init = 2
  , adapt_delta = 0.85 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_VS00_multiVar_LKJ_1022_globalparameters = fit_VS00_multiVar_LKJ_1022$summary(pars %>% append(c('Rho_ID')))
write.csv(fit_VS00_multiVar_LKJ_1022_globalparameters
  , paste0(dropbox_path , "fit_VS00_multiVar_LKJ_1022_globalparameters.csv"), row.names=FALSE)

# Plotting the two fit parameters
fit_VS00_multiVar_LKJ_1022_parameters = data.frame(
  sub = 1:max(ReinforcementLearningStanData_group$sub),
  alpha_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('alpha')$mean,
  alpha_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('alpha')$median,
  beta_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('beta')$mean,
  beta_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('beta')$median,
  soc_vs_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('soc_vs')$mean,
  soc_vs_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('soc_vs')$median,
  theta_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('theta')$mean,
  theta_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('theta')$median,

  log_lik_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_1022$summary('log_lik')$mean,
  waic_VS00_multiVar_LKJ = WAIC_cmdstanr_indv_MCMC(fit_VS00_multiVar_LKJ_1022, nsample=chains*(iter_sampling/thin), nsub=max(ReinforcementLearningStanData_group$sub))$waic
  )
write.csv(fit_VS00_multiVar_LKJ_1022_parameters
  , paste0(dropbox_path , "fit_VS00_multiVar_LKJ_1022_parameters.csv"), row.names=FALSE)

## =======================================================================
## AL00_multiVar_LKJ model fitting
#  =======================================================================
fit_AL00_multiVar_LKJ_1022 = stanmodel_AL00_multiVar_LKJ_multi$sample(

  data = ReinforcementLearningStanData_group
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
  , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_AL00_multiVar_LKJ_1022_globalparameters = fit_AL00_multiVar_LKJ_1022$summary(pars_AL %>% append(c('Rho_ID')))
write.csv(fit_AL00_multiVar_LKJ_1022_globalparameters
  , paste0(dropbox_path , "fit_AL00_multiVar_LKJ_1022_globalparameters.csv"), row.names=FALSE)

# Plotting the two fit parameters
fit_AL00_multiVar_LKJ_1022_parameters = data.frame(
  sub = 1:max(ReinforcementLearningStanData_group$sub),
  alpha_mean_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_1022$summary('alpha')$mean,
  alpha_median_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_1022$summary('alpha')$median,
  beta_mean_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_1022$summary('beta')$mean,
  beta_median_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_1022$summary('beta')$median,

  log_lik_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_1022$summary('log_lik')$mean,
  waic_AL00_multiVar_LKJ = WAIC_cmdstanr_indv_MCMC(fit_AL00_multiVar_LKJ_1022, nsample=chains*(iter_sampling/thin), nsub=max(ReinforcementLearningStanData_group$sub))$waic
  )
write.csv(fit_AL00_multiVar_LKJ_1022_parameters
  , paste0(dropbox_path , "fit_AL00_multiVar_LKJ_1022_parameters.csv"), row.names=FALSE)

# --- loading other model's results ---
fit_SL00_multiVar_LKJ_1022_parameters = read_csv("~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/fitting_results/fit_SL00_multiVar_LKJ_1022_parameters.csv")
fit_AL00_multiVar_LKJ_1022_parameters = read_csv("~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/model_fitting/fitting_results/fit_AL00_multiVar_LKJ_1022_parameters.csv")
fit_VS00_multiVar_LKJ_1022_parameters = read_csv("~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/model_fitting/fitting_results/fit_VS00_multiVar_LKJ_1022_parameters.csv")
# fit_AL00_multiVar_LKJ_group_riskID11_parameters = read_csv("~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/fitting_results/fit_AL00_multiVar_LKJ_group_riskID11_parameters.csv")
# fit_AL00_multiVar_LKJ_group_riskID12_parameters = read_csv("~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/fitting_results/fit_AL00_multiVar_LKJ_group_riskID12_parameters.csv")


###############################################################
# ---- Bayesian Model Selection -----
# Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ. 
#       Bayesian model selection for group studies. NeuroImage. 2009;46(4):1004–17.
###############################################################
## STEP 1: Using WAIC as an approximation of model log-evidence
#-----------------------------------
WAICtable_1022 <- data.frame(
	subject=as.character(fit_VS00_multiVar_LKJ_1022_parameters$sub),
    amazonID = summary_each_individual$amazonID,
    taskID = summary_each_individual$task,
	WAIC_AL00=fit_AL00_multiVar_LKJ_1022_parameters$waic_AL00_multiVar_LKJ,
	WAIC_VS00=fit_VS00_multiVar_LKJ_1022_parameters$waic_VS00_multiVar_LKJ,
	WAIC_SL00=fit_SL00_multiVar_LKJ_1022_parameters$waic_SL00_multiVar_LKJ
	)
# STEP 2. Running the algorithm which updates estimates of the Dirichlet density parameter, α
#         iteretively until convergence
#------------------------------------------------------------------------------------------
num_model = 3
D_alpha0 = rep(1, num_model)
D_alpha = D_alpha0
for(i in 1:500){
	u_WAIC_1 = exp(-WAICtable_1022$WAIC_AL00 + digamma(D_alpha[1]) - digamma(sum(D_alpha)))
	u_WAIC_2 = exp(-WAICtable_1022$WAIC_VS00 + digamma(D_alpha[2]) - digamma(sum(D_alpha)))
	u_WAIC_3 = exp(-WAICtable_1022$WAIC_SL00 + digamma(D_alpha[3]) - digamma(sum(D_alpha)))
	D_beta = c(sum(u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    ,sum(u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    ,sum(u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    )
	D_alpha <- D_alpha0 + D_beta
	r_WAIC = c(D_alpha[1]/sum(D_alpha)
    ,D_alpha[2]/sum(D_alpha)
    ,D_alpha[3]/sum(D_alpha)
    )
	if(i%%10==0){print(D_alpha)}
	#if(i%%100000==0){print(round(r,digits=10))}
	#if(i%%2000==0){print(rbind(D_alpha,r))}
}
# STEP 3: Plot (Posterior model probability)
# the expected multinomial parameters r and thus the expected likelihood of obtaining the k-th model, 
# i.e. p(mnk=1|r)=Mult(m; 1,r), for any randomly selected subject
# -------------------------------------------------------------------------------------------
rWAIC_for_fig = data.frame(
	models = c('Asocial','Value shaping','Decision biasing'),
	r_WAIC = r_WAIC)
rWAIC_for_fig$models = factor(c('Asocial','Value shaping','Decision biasing'),levels=c('Asocial','Value shaping','Decision biasing'))

# STEP 4: Plot (Exceedance Probability)
# (1) sample: (r1, r2, r3) ~ Dirichlet(alpha_1, alpha_2, alpha_3); 
# (2) Count cases where r3 > r1 and r3 > r2
# -------------------------------------------------------------------------------------------
total_sample_N = 10000
dir_sample = rdirichlet(total_sample_N, D_alpha)
xp_model3_1022 = (which(apply(dir_sample, MARGIN=1, FUN=which.max) == 3) %>% length()) / total_sample_N
xp_model3_1022

(r_plot_task12 <- ggplot(rWAIC_for_fig, aes(models,r_WAIC)) + geom_bar(stat="identity")+
	labs(x='Model', y='Posterior probability', title='The 4-armed tasks') +
    annotate(geom="text", x=3, y=0.85, label=paste('XP = ', xp_model3_1022), size=6) + 
	coord_flip() +
	ylim(c(0,1))+
    myTheme_Arial()+
	theme(
		plot.title = element_text(size=16, family='Arial', face='bold'),
	)
)
ggsave(file = paste0(dropbox_path, "r_plot_task12.png"), plot = r_plot_task12, width=6, height=4.5)


# STEP 5. 
# Individual level plot
#------------------------------------------------------------------------------------------
g_forBarplot_1022 = data.frame(
	subject=as.numeric(rep(WAICtable_1022$sub,3)),
    amazonID=as.numeric(rep(WAICtable_1022$amazonID,3)),
    taskID=as.numeric(rep(WAICtable_1022$taskID,3)),
	model=c(rep('Asocial',nrow(WAICtable_1022)),rep('Value shaping',nrow(WAICtable_1022)),rep('Decision biasing',nrow(WAICtable_1022))),
	g_WAIC=c(u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
        , u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
        , u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
	)
WAICtable_1022$g_1 = u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable_1022$g_2 = u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable_1022$g_3 = u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable_1022$bestModel = matrix(c(WAICtable_1022$g_1,WAICtable_1022$g_2,WAICtable_1022$g_3), ncol = 3) %>% apply(MARGIN = 1, FUN = which.max)
WAICtable_1022$bestModel[which(WAICtable_1022$bestModel==1)] <- 'Asocial'
WAICtable_1022$bestModel[which(WAICtable_1022$bestModel==2)] <- 'Value shaping'
WAICtable_1022$bestModel[which(WAICtable_1022$bestModel==3)] <- 'Decision biasing'
g_forBarplot_1022$bestModel = rep(WAICtable_1022$bestModel, 3)
g_forBarplot_1022$model = factor(g_forBarplot_1022$model,levels=c('Asocial','Value shaping','Decision biasing'))

(individual_g_plot_task11 <- g_forBarplot_1022 %>% dplyr::filter(taskID == 1) %>%
    ggplot(aes(model, g_WAIC, fill = bestModel)) + geom_bar(stat="identity") +
	coord_flip() + facet_wrap(~ subject,ncol = 12) +
	scale_y_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1")) +
    scale_fill_viridis_d(name="Most likely\nlearning model",
						breaks=c("Asocial", "Value shaping", "Decision biasing"),
                        labels=c("AL", "VS", "DB"))+
	labs(x='', y='', title='Posterior belief of models for each individual (the 1-risky 3-safe tasks)') +
	myTheme_Arial()+
    theme(panel.border = element_rect(colour = "black", fill=NA))
)
(individual_g_plot_task12 <- g_forBarplot_1022 %>% dplyr::filter(taskID == 2) %>%
    ggplot(aes(model, g_WAIC, fill = bestModel)) + geom_bar(stat="identity") +
	coord_flip() + facet_wrap(~ subject,ncol = 12) +
	scale_y_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1")) +
    scale_fill_viridis_d(name="Most likely\nlearning model",
						breaks=c("Asocial", "Value shaping", "Decision biasing"),
                        labels=c("AL", "VS", "DB"))+
	labs(x='', y='', title='Posterior belief of models for each individual (the 2-risky 2-safe tasks)') +
	myTheme_Arial()+
    theme(panel.border = element_rect(colour = "black", fill=NA))
)

ggsave(file = paste0(dropbox_path, "individual_g_plot_task11.png"), plot = individual_g_plot_task11, width=12, height=12)
ggsave(file = paste0(dropbox_path, "individual_g_plot_task12.png"), plot = individual_g_plot_task12, width=12, height=12)










################################# 
# The Negative Risk Premium task
#################################
# data loading and cleaning
# directory
storing_directory <- "~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/R_and_R_01/revision_exp/"
# Behavioural data for both sessions (22 and 23 Sep) are in the "210923" file
behaviour_revision_exp_210923 <- read.csv(paste0(storing_directory,"data/behaviour_revision_exp_210923.csv")) # this
behaviour_revision_exp_indiv_210923 <- read.csv(paste0(storing_directory,"data/behaviour_revision_exp_indiv_210923.csv")) # this

allBehaviour_revision_exp = rbind(behaviour_revision_exp_210923, behaviour_revision_exp_indiv_210923)

# --- Data exclusion ---
# Removing debug id's
allBehaviour_revision_exp = allBehaviour_revision_exp %>% dplyr::filter(amazonID != 'INHOUSETEST3')

# The following subjects participated in the same time and all wrote "No Strategy." on the questionnaire.
# Not sure if this happened truly by chance or by a single subject using multiple IDs
# But I eliminate this just in case
botIDs <- c('QJDkPpRWv1bl49SSAAAM' , 'DAE0G7vfINM1SkR1AAAC' , 'Lt0oVpvT600Le_--AAAN' , 'cdyrWF817Wvaw_IWAAAO')
rooms_contaminated <- c('hooDBoR_session_407' , 'BhzTojW_session_400')

allBehaviour_revision_exp = allBehaviour_revision_exp %>% dplyr::filter( !(room %in% rooms_contaminated) )
# --- Data exclusion END ---

# Summary data
groupSummary_revision_exp = allBehaviour_revision_exp %>%
	group_by(room) %>%
	summarise(groupSize = max(groupSize),
		indivOrGroup = min(indivOrGroup))

# visualisation of the behavioural trajectory
allBehaviour_revision_exp$choice_num = 0
allBehaviour_revision_exp$choice_num[allBehaviour_revision_exp$choice=='risky'] = 1
allBehaviour_revision_exp$choice_num[allBehaviour_revision_exp$choice=='miss'] = -1

allBehaviour_revision_exp$group_category = allBehaviour_revision_exp$room
allBehaviour_revision_exp$group_category[allBehaviour_revision_exp$indivOrGroup==0] <- 'Individual'

allBehaviour_revision_exp$indivOrGroup_factor = 'Individual'
allBehaviour_revision_exp$indivOrGroup_factor[allBehaviour_revision_exp$indivOrGroup == 1] = 'Social'

allBehaviour_revision_exp %>%
	group_by(round, group_category) %>%
	summarise(mean_risky_choice_prob = mean(choice_num, na.rm=TRUE)
		, groupSize = max(groupSize)
		) -> allBehaviour_revision_exp_each_group
allBehaviour_revision_exp_each_group$groupSize_with_n <- paste('n =', allBehaviour_revision_exp_each_group$groupSize)

allBehaviour_revision_exp_each_group %>%
	ggplot(aes(round, mean_risky_choice_prob)) +
	geom_line(aes(group=group_category), alpha = 1/2, colour='grey20') +
	stat_summary(fun = mean, geom="line", colour='red')+
	# stat_summary(aes(round, 1 - mean_risky_choice_prob), fun = mean, geom="line", colour='blue')+
	geom_hline(yintercept=0.5, linetype='dashed') +
	# geom_segment(data=horizontal_lines, aes(x=1,xend=70,y=yintercept,yend=yintercept), linetype="dashed", size=0.5) +
	facet_grid(groupSize_with_n ~ .) +
	myTheme_Helvetica() +
	scale_y_continuous(breaks=c(0,0.5,1))+
	labs(x = 'Trial'
		, y = 'Proportion of choosing\nthe best risky option'
		, title = 'The 1-risky-1-safe task\n Risk premium = 2/3 \n(N = 118)')+
	NULL -> each_group_behaviour

allBehaviour_revision_exp %>% ggplot() +
	stat_summary(aes(round, choice_num, colour=indivOrGroup, group=indivOrGroup), geom='line', fun=mean) +
	# facet_grid(indivOrGroup ~ .) +
	geom_hline(yintercept = 0.5, linetype = 'dashed')+
	myTheme_Helvetica() +
	NULL


groupSummary_revision_exp %>% dplyr::filter(groupSize>1) %>%
	ggplot() +
	geom_bar(aes(groupSize), stat = "count")+
	labs(x = 'Group size', y = 'Count')+
	myTheme_Helvetica() +
	scale_y_continuous(breaks=c(0,4,8,12))+
	xlim(c(1,7.5))+
	theme(panel.grid.major = element_line(size = 0.5
		, linetype = 'solid', colour='grey40'))+
	NULL -> group_size_distribution

## ==========  Data reading and cleaning
# data for stan fitting - social learning
completedIDs = which(table(allBehaviour_revision_exp$amazonID) >= 36) %>% names()
allBehaviour_revision_exp_social = allBehaviour_revision_exp %>% dplyr::filter(amazonID %in% completedIDs) %>%
  dplyr::filter(indivOrGroup == 1) # note this is only the social condition
# transforming id and room names to numbers
allBehaviour_revision_exp_social$sub = as.numeric(as.factor(allBehaviour_revision_exp_social$amazonID))
allBehaviour_revision_exp_social$group = as.numeric(as.factor(allBehaviour_revision_exp_social$room))

# insert missed trials
for(i in 1:length(table(allBehaviour_revision_exp_social$amazonID))) {
    thisSubject <- allBehaviour_revision_exp_social %>% dplyr::filter(sub==i)
    for(t in 1:70) {
        if(nrow(thisSubject[which(thisSubject$round==t),])==0) {
            newRow = rep(NA, ncol(allBehaviour_revision_exp_social))
            names(newRow) <- names(allBehaviour_revision_exp_social)
            allBehaviour_revision_exp_social <- bind_rows(allBehaviour_revision_exp_social, newRow)
            allBehaviour_revision_exp_social$choice_num[nrow(allBehaviour_revision_exp_social)] = -1
            allBehaviour_revision_exp_social$payoff[nrow(allBehaviour_revision_exp_social)] = -1
            allBehaviour_revision_exp_social$round[nrow(allBehaviour_revision_exp_social)] = t
            allBehaviour_revision_exp_social$socialFreq_safe[nrow(allBehaviour_revision_exp_social)] = 0
            allBehaviour_revision_exp_social$socialFreq_risky[nrow(allBehaviour_revision_exp_social)] = 0
            allBehaviour_revision_exp_social$amazonID[nrow(allBehaviour_revision_exp_social)] = allBehaviour_revision_exp_social$amazonID[which(allBehaviour_revision_exp_social$sub==i)][1]
            allBehaviour_revision_exp_social$sub[nrow(allBehaviour_revision_exp_social)] = allBehaviour_revision_exp_social$sub[which(allBehaviour_revision_exp_social$sub==i)][1]
            allBehaviour_revision_exp_social$group[nrow(allBehaviour_revision_exp_social)] = allBehaviour_revision_exp_social$group[which(allBehaviour_revision_exp_social$sub==i)][1]
            #allBehaviour_revision_exp_social$taskDifficulty_num[nrow(allBehaviour_revision_exp_social)] = allBehaviour_revision_exp_social$taskDifficulty_num[which(allBehaviour_revision_exp_social$sub==i)][1]
        }
    }
}
allBehaviour_revision_exp_social = allBehaviour_revision_exp_social[order(allBehaviour_revision_exp_social$round),] # Sorting by round
allBehaviour_revision_exp_social = allBehaviour_revision_exp_social %>% group_by(amazonID) %>% arrange(round, .by_group = TRUE)

# preparing data fed to Stan
stanData_revision_exp_group = list(
    All = nrow(allBehaviour_revision_exp_social),
    Nsub = length(table(allBehaviour_revision_exp_social$amazonID)),
    Ncue = 2, # number of options
    Ntrial = 70,
    Ncon = 1,
    sub = allBehaviour_revision_exp_social$sub,
    Y = allBehaviour_revision_exp_social$choice_num + 1,
    trial = allBehaviour_revision_exp_social$round,
    payoff = allBehaviour_revision_exp_social$payoff / 100
    )

## compfun() で高速化
library(compiler)
F_calculation = function(array, Nsub, Ntrial, data)
{
    F <- array
    # F (Simulation data とは t+1 ずれてるから注意！)
    for(i in 1:Nsub) {
        F[i,1,1] = 0; F[i,2,1] = 0; # F[i,3,1] = 0;
        for(t in 2:Ntrial) {
            lastChoice = 0
            if( subset(data,sub==i&round==(t-1))$choice_num+1>0 ) {
                lastChoice = subset(data,sub==i&round==(t-1))$choice_num + 1
            }
            F[i,1,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_safe1
            F[i,2,t] <- data %>% dplyr::filter(sub==i&round==t) %>% .$socialFreq_safe2
            #F[i,3,t] = subset(data,sub==i&round==t)$socialFreq2
            if(lastChoice>0){
                F[i,lastChoice,t] = F[i,lastChoice,t] - 1
            }
            if(length(which(F[i,,t]>0))==0){
                F[i,,t] <- c(-1,-1)
            }
        }
    }
    return(F)
}

F_calculation.compiled <- cmpfun(F_calculation)
F0 = array(rep(NA,nrow(allBehaviour_revision_exp_social)), c(stanData_revision_exp_group$Nsub, stanData_revision_exp_group$Ncue, stanData_revision_exp_group$Ntrial))
F = F_calculation.compiled(F0, stanData_revision_exp_group$Nsub, stanData_revision_exp_group$Ntrial,allBehaviour_revision_exp_social)

stanData_revision_exp_group$F = F
## ==========  Data cleaning END


# --- Model fitting VS00 model ---
fit_VS00_multiVar_LKJ_revision_exp = stanmodel_VS00_multiVar_LKJ$sample(

  data = stanData_revision_exp_group
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
  , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_VS00_multiVar_LKJ_revision_exp_globalparameters = fit_VS00_multiVar_LKJ_revision_exp$summary(pars %>% append(c('Rho_ID')))
write.csv(fit_VS00_multiVar_LKJ_revision_exp_globalparameters
  , paste0(dropbox_path , "fit_VS00_multiVar_LKJ_revision_exp_globalparameters.csv"), row.names=FALSE)

# Plotting the two fit parameters
fit_VS00_multiVar_LKJ_revision_exp_parameters = data.frame(
  sub = 1:max(stanData_revision_exp_group$sub),
  alpha_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('alpha')$mean,
  alpha_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('alpha')$median,
  beta_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('beta')$mean,
  beta_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('beta')$median,
  soc_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('soc_vs')$mean,
  soc_vs_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('soc_vs')$median,
  theta_mean_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('theta')$mean,
  theta_median_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('theta')$median,

  log_lik_VS00_multiVar_LKJ = fit_VS00_multiVar_LKJ_revision_exp$summary('log_lik')$mean,
  waic_VS00_multiVar_LKJ = WAIC_cmdstanr_indv_MCMC(fit_VS00_multiVar_LKJ_revision_exp, nsample=chains*(iter_sampling/thin), nsub=max(stanData_revision_exp_group$sub))$waic
  # waic_VS00_multiVar_LKJ = WAIC_cmdstanr_indv_VB(fit_VS00_multiVar_LKJ_revision_exp)$waic
  )

write.csv(fit_VS00_multiVar_LKJ_revision_exp_parameters, paste0(dropbox_path, "fit_VS00_multiVar_LKJ_revision_exp_parameters.csv"), row.names=FALSE)

# --- Model fitting AL00 model ---
fit_AL00_multiVar_LKJ_revision_exp = stanmodel_AL00_multiVar_LKJ$sample(

  data = stanData_revision_exp_group
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
  , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_AL00_multiVar_LKJ_revision_exp_globalparameters = fit_AL00_multiVar_LKJ_revision_exp$summary(pars_AL %>% append(c('Rho_ID')))
write.csv(fit_AL00_multiVar_LKJ_revision_exp_globalparameters
  , paste0(dropbox_path , "fit_AL00_multiVar_LKJ_revision_exp_globalparameters.csv"), row.names=FALSE)

# Plotting the two fit parameters
fit_AL00_multiVar_LKJ_revision_exp_parameters = data.frame(
  sub = 1:max(stanData_revision_exp_group$sub),
  alpha_mean_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_revision_exp$summary('alpha')$mean,
  alpha_median_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_revision_exp$summary('alpha')$median,
  beta_mean_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_revision_exp$summary('beta')$mean,
  beta_median_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_revision_exp$summary('beta')$median,

  log_lik_AL00_multiVar_LKJ = fit_AL00_multiVar_LKJ_revision_exp$summary('log_lik')$mean,
  waic_AL00_multiVar_LKJ = WAIC_cmdstanr_indv_MCMC(fit_AL00_multiVar_LKJ_revision_exp, nsample=chains*(iter_sampling/thin), nsub=max(stanData_revision_exp_group$sub))$waic
  # waic_AL00_multiVar_LKJ = WAIC_cmdstanr_indv_VB(fit_AL00_multiVar_LKJ_revision_exp)$waic
  )

write.csv(fit_AL00_multiVar_LKJ_revision_exp_parameters, paste0(dropbox_path, "fit_AL00_multiVar_LKJ_revision_exp_parameters.csv"), row.names=FALSE)


# loading fit results of other model
fit_AL00_multiVar_LKJ_revision_exp_parameters = read_csv(paste0(dropbox_path, "fit_AL00_multiVar_LKJ_revision_exp_parameters.csv"))
fit_VS00_multiVar_LKJ_revision_exp_parameters = read_csv(paste0(dropbox_path, "fit_VS00_multiVar_LKJ_revision_exp_parameters.csv"))
fit_SL00_multiVar_LKJ_revision_exp_parameters = read_csv('~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/R_and_R_01/revision_exp/fitting_results/fit_SL00_multiVar_LKJ_revision_exp_parameters.csv')


###############################################################
# ---- Bayesian Model Selection -----
###############################################################
## STEP 1: Using WAIC as an approximation of model log-evidence
#-----------------------------------
WAICtable_revision_exp <- data.frame(
	subject=as.character(fit_VS00_multiVar_LKJ_revision_exp_parameters$sub),
	WAIC_AL00=fit_AL00_multiVar_LKJ_revision_exp_parameters$waic_AL00_multiVar_LKJ,
	WAIC_VS00=fit_VS00_multiVar_LKJ_revision_exp_parameters$waic_VS00_multiVar_LKJ,
	WAIC_SL00=fit_SL00_multiVar_LKJ_revision_exp_parameters$waic_SL00_multiVar_LKJ
	)
# STEP 2. Running the algorithm which updates estimates of the Dirichlet density parameter, α
#         iteretively until convergence
#------------------------------------------------------------------------------------------
num_model = 3
D_alpha0 = rep(1, num_model)
D_alpha = D_alpha0
for(i in 1:500){
	u_WAIC_1 = exp(-WAICtable_revision_exp$WAIC_AL00 + digamma(D_alpha[1]) - digamma(sum(D_alpha)))
	u_WAIC_2 = exp(-WAICtable_revision_exp$WAIC_VS00 + digamma(D_alpha[2]) - digamma(sum(D_alpha)))
	u_WAIC_3 = exp(-WAICtable_revision_exp$WAIC_SL00 + digamma(D_alpha[3]) - digamma(sum(D_alpha)))
	D_beta = c(sum(u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    ,sum(u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    ,sum(u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
    )
	D_alpha <- D_alpha0 + D_beta
	r_WAIC = c(D_alpha[1]/sum(D_alpha)
    ,D_alpha[2]/sum(D_alpha)
    ,D_alpha[3]/sum(D_alpha)
    )
	if(i%%10==0){print(D_alpha)}
	#if(i%%100000==0){print(round(r,digits=10))}
	#if(i%%2000==0){print(rbind(D_alpha,r))}
}
# STEP 3: Plot (Posterior model probability)
# the expected multinomial parameters r and thus the expected likelihood of obtaining the k-th model, 
# i.e. p(mnk=1|r)=Mult(m; 1,r), for any randomly selected subject
# -------------------------------------------------------------------------------------------
rWAIC_for_fig = data.frame(
	models = c('Asocial','Value shaping','Decision biasing'),
	r_WAIC = r_WAIC)
rWAIC_for_fig$models = factor(c('Asocial','Value shaping','Decision biasing'),levels=c('Asocial','Value shaping','Decision biasing'))


# STEP 4: Plot (Exceedance Probability)
# (1) sample: (r1, r2, r3) ~ Dirichlet(alpha_1, alpha_2, alpha_3); 
# (2) Count cases where r3 > r1 and r3 > r2
# -------------------------------------------------------------------------------------------
total_sample_N = 10000
dir_sample = rdirichlet(total_sample_N, D_alpha)
xp_model3_revision_exp = (which(apply(dir_sample, MARGIN=1, FUN=which.max) == 3) %>% length()) / total_sample_N
xp_model3_revision_exp

(r_plot_revision_exp <- ggplot(rWAIC_for_fig, aes(models,r_WAIC)) + geom_bar(stat="identity")+
	labs(x='Model', y='Posterior probability', title='The negative risk premium task') +
    annotate(geom="text", x=3, y=0.85, label=paste('XP = ', floor(100*xp_model3_revision_exp)/100), size=6) + 
	coord_flip() +
	ylim(c(0,1))+
    myTheme_Arial()+
	theme(
		plot.title = element_text(size=16, family='Arial', face='bold'),
	)
)
ggsave(file = paste0(dropbox_path, "r_plot_revision_exp.png"), plot = r_plot_revision_exp, width=6, height=4.5)

# STEP 5. 
# Individual level plot
#------------------------------------------------------------------------------------------
g_forBarplot_revision_exp = data.frame(
	subject=as.numeric(rep(WAICtable_revision_exp$sub,3)),
	model=c(rep('Asocial',nrow(WAICtable_revision_exp)),rep('Value shaping',nrow(WAICtable_revision_exp)),rep('Decision biasing',nrow(WAICtable_revision_exp))),
	g_WAIC=c(u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
        , u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
        , u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3))
	)
WAICtable_revision_exp$g_1 = u_WAIC_1/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable_revision_exp$g_2 = u_WAIC_2/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable_revision_exp$g_3 = u_WAIC_3/(u_WAIC_1+u_WAIC_2+u_WAIC_3)
WAICtable_revision_exp$bestModel = matrix(c(WAICtable_revision_exp$g_1,WAICtable_revision_exp$g_2,WAICtable_revision_exp$g_3), ncol = 3) %>% apply(MARGIN = 1, FUN = which.max)
WAICtable_revision_exp$bestModel[which(WAICtable_revision_exp$bestModel==1)] <- 'Asocial'
WAICtable_revision_exp$bestModel[which(WAICtable_revision_exp$bestModel==2)] <- 'Value shaping'
WAICtable_revision_exp$bestModel[which(WAICtable_revision_exp$bestModel==3)] <- 'Decision biasing'
g_forBarplot_revision_exp$bestModel = rep(WAICtable_revision_exp$bestModel, 3)
g_forBarplot_revision_exp$model = factor(g_forBarplot_revision_exp$model,levels=c('Asocial','Value shaping','Decision biasing'))

(individual_g_plot_NRP <- g_forBarplot_revision_exp %>% 
    ggplot(aes(model, g_WAIC, fill = bestModel)) + geom_bar(stat="identity") +
	coord_flip() + facet_wrap(~ subject,ncol = 12) +
	scale_y_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1")) +
    scale_fill_viridis_d(name="Most likely\nlearning model",
						breaks=c("Asocial", "Value shaping", "Decision biasing"),
                        labels=c("AL", "VS", "DB"))+
	labs(x='', y='', title='Posterior belief of models for each individual (the NRP tasks)') +
	myTheme_Arial()+
    theme(panel.border = element_rect(colour = "black", fill=NA))
)


ggsave(file = paste0(dropbox_path, "individual_g_plot_NRP.png"), plot = individual_g_plot_NRP, width=12, height=12)



exceedance_probability = data.frame(
    task = c('2ab','4ab','negative')
    , xp_decision_biasing = c(xp_model3, xp_model3_1022, xp_model3_revision_exp)
    )

write.csv(exceedance_probability, paste0(dropbox_path, "exceedance_probability.csv"), row.names=FALSE)

# Merging all three plots
model_probability_plot =  plot_grid(r_plot_0820
    , r_plot_task12
    , r_plot_revision_exp
    # ,  common.legend = TRUE
    # ,  legend = 'right'
    , axis = 'bt'
    , labels = c('a','b','c'), ncol = 3, align = 'v'
    )

ggsave(file = paste0(dropbox_path, "model_probability_plot.png"), plot = model_probability_plot, width=18, height=5)









