###############################################################################################
##
##
## Experimental Plots for the 2nd revision 
##
## Wataru Toyokawa
## What's new?
## (1) Figure 6 now contains only the expetimental results with sample aberage curve 
## (2) A new Figure 7 shows the model prediction curves and the histograms  
###############################################################################################

if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace

# Loading
library(tidyverse)
library(ggpubr)
library(cmdstanr)
library(rstanarm)
library(tidybayes)
library(posterior)
library(bayesplot)
library(loo)
library(foreach)
library(MASS)
library(doParallel)
library(lme4)
library(merTools)
library(cowplot)
#registerDoParallel(detectCores()) # this uses as many core as available
registerDoParallel(detectCores())
# rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores()) # regacy from rstan?

# directory
storing_directory <- "~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/"

## Load Functions
source('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_0820/analysis/functions.R')
source('~/Dropbox/wataru/papers/RiskySocialLearning/four_armed_simulation/functions.R')

# --- Behavioural data -- 
# data loading and cleaning
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
botIDs <- c('QJDkPpRWv1bl49SSAAAM' , 'DAE0G7vfINM1SkR1AAAC' , 'Lt0oVpvT600Le_--AAAN' , 'cdyrWF817Wvaw_IWAAAO', 'INHOUSETEST3')
rooms_contaminated <- c('hooDBoR_session_407' , 'BhzTojW_session_400')

allBehaviour_revision_exp = allBehaviour_revision_exp %>% dplyr::filter( !(room %in% rooms_contaminated) )
# --- Data exclusion END ---

# --- credible interbals ---
quantile_probs = seq(0, 1, 0.1)


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


# --- model fitting and post-hoc simulation results
social_learning_model_validation_revision_exp_data <- read_csv(paste0(storing_directory, "fitting_results/social_learning_model_validation_revision_exp_data.csv"))
fit_AL00_multiVar_LKJ_revision_exp_indiv_globalparameters <- read_csv(paste0(storing_directory, "fitting_results/fit_AL00_multiVar_LKJ_revision_exp_indiv_globalparameters.csv"))
fit_AL00_multiVar_LKJ_revision_exp_indiv_parameters <- read_csv(paste0(storing_directory, "fitting_results/fit_AL00_multiVar_LKJ_revision_exp_indiv_parameters.csv"))
fit_SL00_multiVar_LKJ_revision_exp_globalparameters <- read_csv(paste0(storing_directory, "fitting_results/fit_SL00_multiVar_LKJ_revision_exp_globalparameters.csv"))
fit_SL00_multiVar_LKJ_revision_exp_parameters <- read_csv(paste0(storing_directory, "fitting_results/fit_SL00_multiVar_LKJ_revision_exp_parameters.csv"))

# --- data from the questionnaires
Test_riskyBandit_202109221429 <- read_csv(paste0(storing_directory, "data/Test_riskyBandit_202109221429.csv"))
Test_riskyBandit_202109230855 <- read_csv(paste0(storing_directory, "data/Test_riskyBandit_202109230855.csv"))
Test_riskyBandit_202109231201 <- read_csv(paste0(storing_directory, "data/Test_riskyBandit_202109231201.csv"))

Test_riskyBandit_temp <- rbind(Test_riskyBandit_202109221429, Test_riskyBandit_202109230855) %>% rbind(Test_riskyBandit_202109231201)

# Test_riskyBandit_revision <- Test_riskyBandit_temp %>%
# 	dplyr::filter( amazonID %in% c(fit_parameters_group_SL00_mcmc$amazonID, allBehaviour_revision_exp_indiv_summarised_t35$amazonID) )


# ============================================
# --- Figure 6 -- Figure Supplement 1 ---
# ============================================

allBehaviour_revision_exp$groupSize_category = 'Small'
allBehaviour_revision_exp$groupSize_category[which(allBehaviour_revision_exp$groupSize==1)] = 'Individual'
allBehaviour_revision_exp$groupSize_category[which(allBehaviour_revision_exp$groupSize>4)] = 'Large'
allBehaviour_revision_exp$groupSize_category = factor(allBehaviour_revision_exp$groupSize_category, levels = c('Individual','Small','Large'))

allBehaviour_revision_exp$choice_num2 = 0
allBehaviour_revision_exp$choice_num2[allBehaviour_revision_exp$choice=='risky'] = 1
allBehaviour_revision_exp$choice_num2[allBehaviour_revision_exp$choice=='miss'] = NA

# -- individual condition's data --
completedIDs = which(table(allBehaviour_revision_exp$amazonID) >= 36) %>% names()
allBehaviour_revision_exp_indiv = allBehaviour_revision_exp %>% dplyr::filter(amazonID %in% completedIDs) %>%
  dplyr::filter(indivOrGroup == 0) # note this is only the individual condition
allBehaviour_revision_exp_indiv$sub = as.numeric(as.factor(allBehaviour_revision_exp_indiv$amazonID))
allBehaviour_revision_exp_indiv = allBehaviour_revision_exp_indiv %>% group_by(amazonID) %>% arrange(round, .by_group = TRUE)

# -- summarised data --
allBehaviour_revision_exp_indiv_summarised_t35 = allBehaviour_revision_exp_indiv %>%
	dplyr::filter(round>35) %>%
	group_by(sub) %>%
	summarise(
		risky_choice_count = sum(choice_num2, na.rm = TRUE),
		risky_choice_mean = mean(choice_num2, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1],
		amazonID = amazonID[1]
	)
allBehaviour_revision_exp_indiv_summarised_t35$groupID = allBehaviour_revision_exp_indiv_summarised_t35$room
allBehaviour_revision_exp_indiv_summarised_t35$groupID[which(allBehaviour_revision_exp_indiv_summarised_t35$indivOrGroup_factor=='Individual')] = 'Individual'

# Individual Condition Only
parameterfit_indiv_AL00_revision_exp <- left_join(fit_AL00_multiVar_LKJ_revision_exp_indiv_parameters, allBehaviour_revision_exp_indiv_summarised_t35, by = 'sub')
parameterfit_indiv_AL00_revision_exp$hot_stove_susceptibility <- parameterfit_indiv_AL00_revision_exp$alpha_median_AL00_multiVar_LKJ * (parameterfit_indiv_AL00_revision_exp$beta_median_AL00_multiVar_LKJ + 1)
# if the hot stove effect is too large
parameterfit_indiv_AL00_revision_exp$hot_stove_susceptibility_trancated <- parameterfit_indiv_AL00_revision_exp$hot_stove_susceptibility
parameterfit_indiv_AL00_revision_exp$hot_stove_susceptibility_trancated[which(parameterfit_indiv_AL00_revision_exp$hot_stove_susceptibility > 6)] <- 6

# -- group condition's data --
completedIDs = which(table(allBehaviour_revision_exp$amazonID) >= 36) %>% names()
allBehaviour_revision_exp_social = allBehaviour_revision_exp %>% dplyr::filter(amazonID %in% completedIDs) %>%
  dplyr::filter(indivOrGroup == 1) # note this is only the social condition
# transforming id and room names to numbers
allBehaviour_revision_exp_social$sub = as.numeric(as.factor(allBehaviour_revision_exp_social$amazonID))
allBehaviour_revision_exp_social$group = as.numeric(as.factor(allBehaviour_revision_exp_social$room))

# summarised data
allBehaviour_revision_exp_social_summarised_t35 = allBehaviour_revision_exp_social %>%
	dplyr::filter(round>35) %>%
	group_by(sub) %>%
	summarise(
		risky_choice_count = sum(choice_num2, na.rm = TRUE),
		risky_choice_mean = mean(choice_num2, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1],
		amazonID = amazonID[1]
	)
allBehaviour_revision_exp_social_summarised_t35$groupID = allBehaviour_revision_exp_social_summarised_t35$room
allBehaviour_revision_exp_social_summarised_t35$groupID[which(allBehaviour_revision_exp_social_summarised_t35$indivOrGroup_factor=='Individual')] = 'Individual'

# Group condition
fit_parameters_revision_SL00_mcmc <- left_join(fit_SL00_multiVar_LKJ_revision_exp_parameters, allBehaviour_revision_exp_social_summarised_t35, by = 'sub')
fit_parameters_revision_SL00_mcmc$hot_stove_susceptibility <- fit_parameters_revision_SL00_mcmc$alpha_median_SL00_multiVar_LKJ * (fit_parameters_revision_SL00_mcmc$beta_median_SL00_multiVar_LKJ + 1)
# if the hot stove effect is too large
fit_parameters_revision_SL00_mcmc$hot_stove_susceptibility_trancated <- fit_parameters_revision_SL00_mcmc$hot_stove_susceptibility
fit_parameters_revision_SL00_mcmc$hot_stove_susceptibility_trancated[which(fit_parameters_revision_SL00_mcmc$hot_stove_susceptibility > 6)] <- 6
fit_parameters_revision_SL00_mcmc$soc_mean <- fit_parameters_revision_SL00_mcmc$soc_mean_SL00_multiVar_LKJ


# overall means
social_learning_model_validation_revision_exp_data$raw_proportionRiskyChoice_b2 = social_learning_model_validation_revision_exp_data$proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw()
social_learning_model_validation_revision_exp_summary <-
	social_learning_model_validation_revision_exp_data %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[2] %>% convert_alphaRaw_to_alpha,
        proportionRiskyChoice_b2_upr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[8] %>% convert_alphaRaw_to_alpha,
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% mean(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_revision_exp_summary$proportionRiskyChoice_b2_lower <-
	(social_learning_model_validation_revision_exp_summary$raw_proportionRiskyChoice_b2_mean - social_learning_model_validation_revision_exp_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_revision_exp_summary$n)) %>% convert_alphaRaw_to_alpha
social_learning_model_validation_revision_exp_summary$proportionRiskyChoice_b2_upper <-
	(social_learning_model_validation_revision_exp_summary$raw_proportionRiskyChoice_b2_mean + social_learning_model_validation_revision_exp_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_revision_exp_summary$n)) %>% convert_alphaRaw_to_alpha

social_learning_model_validation_revision_exp_summary$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_revision_exp_summary$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha


# modest social learners' means
sigma_quantiles_revision_exp = social_learning_model_validation_revision_exp_data %>% filter(condition_dummy==1) %>% .$soc_mean %>% quantile(probs = quantile_probs)
social_learning_model_validation_revision_exp_summary_sigma_level_1 <-
	social_learning_model_validation_revision_exp_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_revision_exp[1] & soc_mean < sigma_quantiles_revision_exp[3] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_revision_exp_summary_sigma_level_1$sigma_mean = mean(social_learning_model_validation_revision_exp_summary_sigma_level_1$soc_mean)

social_learning_model_validation_revision_exp_summary_sigma_level_2 <-
	social_learning_model_validation_revision_exp_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_revision_exp[3] & soc_mean < sigma_quantiles_revision_exp[5] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_revision_exp_summary_sigma_level_2$sigma_mean = mean(social_learning_model_validation_revision_exp_summary_sigma_level_2$soc_mean)

social_learning_model_validation_revision_exp_summary_sigma_level_3 <-
	social_learning_model_validation_revision_exp_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_revision_exp[5] & soc_mean < sigma_quantiles_revision_exp[7] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_revision_exp_summary_sigma_level_3$sigma_mean = mean(social_learning_model_validation_revision_exp_summary_sigma_level_3$soc_mean)

social_learning_model_validation_revision_exp_summary_sigma_level_4 <-
	social_learning_model_validation_revision_exp_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_revision_exp[7] & soc_mean < sigma_quantiles_revision_exp[9] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_revision_exp_summary_sigma_level_4$sigma_mean = mean(social_learning_model_validation_revision_exp_summary_sigma_level_4$soc_mean)

social_learning_model_validation_revision_exp_summary_sigma_level_5 <-
	social_learning_model_validation_revision_exp_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_revision_exp[9] & soc_mean < sigma_quantiles_revision_exp[11] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_revision_exp_summary_sigma_level_5$sigma_mean = mean(social_learning_model_validation_revision_exp_summary_sigma_level_5$soc_mean)

social_learning_model_validation_revision_exp_summary_sigma_merged = 
    social_learning_model_validation_revision_exp_summary_sigma_level_1 %>%
    rbind(social_learning_model_validation_revision_exp_summary_sigma_level_2) %>%
    rbind(social_learning_model_validation_revision_exp_summary_sigma_level_3) %>%
    rbind(social_learning_model_validation_revision_exp_summary_sigma_level_4) %>%
    rbind(social_learning_model_validation_revision_exp_summary_sigma_level_5)



# modest social learners' means
social_learning_model_validation_revision_exp_summary_reallyHighSigma <-
	social_learning_model_validation_revision_exp_data %>%
	dplyr::filter(soc_mean > 4/10 & soc_mean < 10/10 & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[2] %>% convert_alphaRaw_to_alpha,
        proportionRiskyChoice_b2_upr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[8] %>% convert_alphaRaw_to_alpha,
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_revision_exp_summary_reallyHighSigma$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_revision_exp_summary_reallyHighSigma$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha



## ========================================================
#
# Figure for the experimental results
#
# =========================================================
library(ggpubr)
main_exp_dir <- '~/Dropbox/wataru/papers/RiskySocialLearning/draft/analysis_repo/'
social_learning_model_validation_0820_data <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/social_learning_model_validation_0820_data.csv'))
social_learning_model_validation_1022_riskID11_data <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/social_learning_model_validation_1022_riskID11_data.csv'))
social_learning_model_validation_1022_riskID12_data <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/social_learning_model_validation_1022_riskID12_data.csv'))
# ======================================
# 1-risky 1-safe (2-armed) task
# ======================================
fit_AL00_multiVar_LKJ_indiv_0820_globalparameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_0820_globalparameters.csv'))
fit_AL00_multiVar_LKJ_indiv_0820_parameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_0820_parameters.csv'))
fit_SL00_multiVar_LKJ_0820_globalparameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_SL00_multiVar_LKJ_0820_globalparameters.csv'))
fit_SL00_multiVar_LKJ_0820_parameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_SL00_multiVar_LKJ_0820_parameters.csv'))

behaviour_main_0820 <- read_csv(paste0(main_exp_dir, "experimentalAnalysis/behaviour_main_0820.csv"))
behaviour_indiv_0820 <-read_csv(paste0(main_exp_dir, "experimentalAnalysis/behaviour_indiv_0820.csv"))
allBehaviour0820 <- rbind(behaviour_main_0820, behaviour_indiv_0820)
allBehaviour0820 <- allBehaviour0820 %>%
  dplyr::filter(amazonID != 'INHOUSETEST2') %>%  # eliminating data generated by debug tests
  dplyr::filter(amazonID != '5eac70db94edd22d57fa00c4') # a bug in the data storaging process

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

# individual condition
completedIDs = which(table(allBehaviour0820$amazonID) >= 36) %>% names()
allBehaviour0820_indiv = allBehaviour0820 %>% dplyr::filter(amazonID %in% completedIDs) %>%
  dplyr::filter(indivOrGroup == 0) # note this is only the individual condition
allBehaviour0820_indiv$sub = as.numeric(as.factor(allBehaviour0820_indiv$amazonID))
allBehaviour0820_indiv = allBehaviour0820_indiv %>% group_by(amazonID) %>% arrange(round, .by_group = TRUE)

# summarised data
allBehaviour0820_indiv_summarised_t35 = allBehaviour0820_indiv %>%
	dplyr::filter(round>35) %>%
	group_by(sub) %>%
	summarise(
		risky_choice_count = sum(choice_num, na.rm = TRUE),
		risky_choice_mean = mean(choice_num, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1],
		amazonID = amazonID[1]
	)
allBehaviour0820_indiv_summarised_t35$groupID = allBehaviour0820_indiv_summarised_t35$room
allBehaviour0820_indiv_summarised_t35$groupID[which(allBehaviour0820_indiv_summarised_t35$indivOrGroup_factor=='Individual')] = 'Individual'

# Individual Condition Only
parameterfit_indiv_AL00_0820 <- left_join(fit_AL00_multiVar_LKJ_indiv_0820_parameters, allBehaviour0820_indiv_summarised_t35, by = 'sub')
parameterfit_indiv_AL00_0820$hot_stove_susceptibility <- parameterfit_indiv_AL00_0820$alpha_median_AL00_multiVar_LKJ * (parameterfit_indiv_AL00_0820$beta_median_AL00_multiVar_LKJ + 1)
# if the hot stove effect is too large
parameterfit_indiv_AL00_0820$hot_stove_susceptibility_trancated <- parameterfit_indiv_AL00_0820$hot_stove_susceptibility
parameterfit_indiv_AL00_0820$hot_stove_susceptibility_trancated[which(parameterfit_indiv_AL00_0820$hot_stove_susceptibility > 6)] <- 6


# Group condition
completedIDs = which(table(allBehaviour0820$amazonID) >= 36) %>% names()
allBehaviour0820_social = allBehaviour0820 %>% dplyr::filter(amazonID %in% completedIDs) %>%
  dplyr::filter(indivOrGroup == 1) # note this is only the social condition
allBehaviour0820_social$sub = as.numeric(as.factor(allBehaviour0820_social$amazonID))
allBehaviour0820_social$group = as.numeric(as.factor(allBehaviour0820_social$room))
# summarised data
allBehaviour0820_social_summarised_t35 = allBehaviour0820_social %>%
	dplyr::filter(round>35) %>%
	group_by(sub) %>%
	summarise(
		risky_choice_count = sum(choice_num, na.rm = TRUE),
		risky_choice_mean = mean(choice_num, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1],
		amazonID = amazonID[1]
	)
allBehaviour0820_social_summarised_t35$groupID = allBehaviour0820_social_summarised_t35$room
allBehaviour0820_social_summarised_t35$groupID[which(allBehaviour0820_social_summarised_t35$indivOrGroup_factor=='Individual')] = 'Individual'

# Group condition
fit_parameters_group_SL00_mcmc <- left_join(fit_SL00_multiVar_LKJ_0820_parameters, allBehaviour0820_social_summarised_t35, by = 'sub')
fit_parameters_group_SL00_mcmc$hot_stove_susceptibility <- fit_parameters_group_SL00_mcmc$alpha_median_SL00_multiVar_LKJ * (fit_parameters_group_SL00_mcmc$beta_median_SL00_multiVar_LKJ + 1)
# if the hot stove effect is too large
fit_parameters_group_SL00_mcmc$hot_stove_susceptibility_trancated <- fit_parameters_group_SL00_mcmc$hot_stove_susceptibility
fit_parameters_group_SL00_mcmc$hot_stove_susceptibility_trancated[which(fit_parameters_group_SL00_mcmc$hot_stove_susceptibility > 6)] <- 6

# overall means
social_learning_model_validation_0820_data$raw_proportionRiskyChoice_b2 = social_learning_model_validation_0820_data$proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw()
social_learning_model_validation_0820_summary <-
	social_learning_model_validation_0820_data %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[2] %>% convert_alphaRaw_to_alpha,
        proportionRiskyChoice_b2_upr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[8] %>% convert_alphaRaw_to_alpha,
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
        # # soc_sd = sd(soc_mean),
        # soc_quantile_1 = quantile(soc_mean, probs = quantile_probs)[2],
        # soc_quantile_2 = quantile(soc_mean, probs = quantile_probs)[5],
        # soc_quantile_3 = quantile(soc_mean, probs = quantile_probs)[9],
		n = n()
		)

social_learning_model_validation_0820_summary$proportionRiskyChoice_b2_lower <-
	(social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_mean - social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_0820_summary$n)) %>% convert_alphaRaw_to_alpha
social_learning_model_validation_0820_summary$proportionRiskyChoice_b2_upper <-
	(social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_mean + social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_0820_summary$n)) %>% convert_alphaRaw_to_alpha

social_learning_model_validation_0820_summary$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha


# modest social learners' means
sigma_quantiles_0820 = social_learning_model_validation_0820_data %>% filter(condition_dummy==1) %>% .$soc_mean %>% quantile(probs = quantile_probs)
social_learning_model_validation_0820_summary_sigma_level_1 <-
	social_learning_model_validation_0820_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_0820[1] & soc_mean < sigma_quantiles_0820[3] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_0820_summary_sigma_level_1$sigma_mean = mean(social_learning_model_validation_0820_summary_sigma_level_1$soc_mean)

social_learning_model_validation_0820_summary_sigma_level_2 <-
	social_learning_model_validation_0820_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_0820[3] & soc_mean < sigma_quantiles_0820[5] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_0820_summary_sigma_level_2$sigma_mean = mean(social_learning_model_validation_0820_summary_sigma_level_2$soc_mean)

social_learning_model_validation_0820_summary_sigma_level_3 <-
	social_learning_model_validation_0820_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_0820[5] & soc_mean < sigma_quantiles_0820[7] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_0820_summary_sigma_level_3$sigma_mean = mean(social_learning_model_validation_0820_summary_sigma_level_3$soc_mean)

social_learning_model_validation_0820_summary_sigma_level_4 <-
	social_learning_model_validation_0820_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_0820[7] & soc_mean < sigma_quantiles_0820[9] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_0820_summary_sigma_level_4$sigma_mean = mean(social_learning_model_validation_0820_summary_sigma_level_4$soc_mean)

social_learning_model_validation_0820_summary_sigma_level_5 <-
	social_learning_model_validation_0820_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_0820[9] & soc_mean < sigma_quantiles_0820[11] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_0820_summary_sigma_level_5$sigma_mean = mean(social_learning_model_validation_0820_summary_sigma_level_5$soc_mean)

social_learning_model_validation_0820_summary_sigma_merged = 
    social_learning_model_validation_0820_summary_sigma_level_1 %>%
    rbind(social_learning_model_validation_0820_summary_sigma_level_2) %>%
    rbind(social_learning_model_validation_0820_summary_sigma_level_3) %>%
    rbind(social_learning_model_validation_0820_summary_sigma_level_4) %>%
    rbind(social_learning_model_validation_0820_summary_sigma_level_5)




social_learning_model_validation_0820_summary_reallyHighSigma <-
	social_learning_model_validation_0820_data %>%
	dplyr::filter(soc_mean > 3/10 & soc_mean < 6/10 & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[2] %>% convert_alphaRaw_to_alpha,
        proportionRiskyChoice_b2_upr = quantile(raw_proportionRiskyChoice_b2, probs = quantile_probs)[8] %>% convert_alphaRaw_to_alpha,
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_0820_summary_reallyHighSigma$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_0820_summary_reallyHighSigma$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha

fit_parameters_group_SL00_mcmc$soc_mean <- fit_parameters_group_SL00_mcmc$soc_mean_SL00_multiVar_LKJ




# ======================================
# 1-risky 3-safe (4-armed) task
# ======================================

# fit result -- global parameters
fit_SL00_multiVar_LKJ_1022_globalparameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_globalparameters.csv'))
fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters.csv'))

## behavioural data summary
allBehaviour1022_group <- read_csv(paste0(main_exp_dir, "experimentalAnalysis/allBehaviour1022_group.csv"))
allBehaviour1022_group_riskID11 <- allBehaviour1022_group %>% dplyr::filter(riskDistributionId == 11)

fit_SL00_multiVar_LKJ_1022_parameters <- read_csv(paste0(main_exp_dir, "experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_parameters.csv"))

allBehaviour1022_group_riskID11_summarised_t35 <- allBehaviour1022_group_riskID11 %>%
	dplyr::filter(round>35) %>%
	group_by(amazonID, sub) %>%
	summarise(
		risky_choice_count = sum(best_risky_choice, na.rm = TRUE),
		risky_choice_mean = mean(best_risky_choice, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1]
	)

allBehaviour1022_group_riskID11_summarised_t35$groupID <- allBehaviour1022_group_riskID11_summarised_t35$room


allBehaviour1022_indiv <- read_csv(paste0(main_exp_dir, "experimentalAnalysis/allBehaviour1022_indiv.csv"))
allBehaviour1022_indiv_riskID11 <- allBehaviour1022_indiv %>%
  filter(riskDistributionId_factor=='Con: 0')
# The 1-risky-3-safe task was labeled "Con: 0" or "11" originally
# And in the analysis code, things like "riskID11" means the 1-risky-3-safe task
# while "riskID12" means the 2-risky-2-safe task
allBehaviour1022_indiv_riskID11$sub_old <- allBehaviour1022_indiv_riskID11$sub
allBehaviour1022_indiv_riskID11$sub <- allBehaviour1022_indiv_riskID11$amazonID %>% as.factor() %>% as.numeric()

allBehaviour1022_indiv_riskID11_summarised_t35 <- allBehaviour1022_indiv_riskID11 %>%
	dplyr::filter(round>35) %>%
	group_by(amazonID, sub) %>%
	summarise(
		risky_choice_count = sum(best_risky_choice, na.rm = TRUE),
		risky_choice_mean = mean(best_risky_choice, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1]
	)
allBehaviour1022_indiv_riskID11_summarised_t35$groupID = allBehaviour1022_indiv_riskID11_summarised_t35$room
allBehaviour1022_indiv_riskID11_summarised_t35$groupID[which(allBehaviour1022_indiv_riskID11_summarised_t35$indivOrGroup_factor=='Individual')] = 'Individual'

# Individual fits
fit_AL00_multiVar_LKJ_indiv_riskID11_parameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID11_parameters.csv'))

# Merging the behavioural data with the fit parameters
# hot stove effect - individual
fit_AL_indiv_riskID11_parameters <- right_join(fit_AL00_multiVar_LKJ_indiv_riskID11_parameters, allBehaviour1022_indiv_riskID11_summarised_t35, by = 'sub')
fit_AL_indiv_riskID11_parameters$hot_stove_susceptibility <- fit_AL_indiv_riskID11_parameters$alpha_median_AL00_multiVar_LKJ * (1+ fit_AL_indiv_riskID11_parameters$beta_median_AL00_multiVar_LKJ)
fit_AL_indiv_riskID11_parameters$hot_stove_susceptibility_trancated <- fit_AL_indiv_riskID11_parameters$hot_stove_susceptibility
fit_AL_indiv_riskID11_parameters$hot_stove_susceptibility_trancated[which(fit_AL_indiv_riskID11_parameters$hot_stove_susceptibility > 6)] <- 6

# hot stove effect - group
fit_SL00_riskID11_parameters <- right_join(fit_SL00_multiVar_LKJ_1022_parameters, allBehaviour1022_group_riskID11_summarised_t35, by = 'sub')
fit_SL00_riskID11_parameters$hot_stove_susceptibility <- fit_SL00_riskID11_parameters$alpha_mean_SL00_multiVar_LKJ * (1+ fit_SL00_riskID11_parameters$beta_mean_SL00_multiVar_LKJ)
fit_SL00_riskID11_parameters$hot_stove_susceptibility_trancated <- fit_SL00_riskID11_parameters$hot_stove_susceptibility
fit_SL00_riskID11_parameters$hot_stove_susceptibility_trancated[which(fit_SL00_riskID11_parameters$hot_stove_susceptibility > 6)] <- 6

# overall means
social_learning_model_validation_1022_riskID11_summary <-
	social_learning_model_validation_1022_riskID11_data %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[2],
        proportionRiskyChoice_b2_upr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[8],
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% mean(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_1022_riskID11_summary$proportionRiskyChoice_b2_lower <-
	(social_learning_model_validation_1022_riskID11_summary$raw_proportionRiskyChoice_b2_mean - social_learning_model_validation_1022_riskID11_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_1022_riskID11_summary$n)) %>% convert_alphaRaw_to_alpha
social_learning_model_validation_1022_riskID11_summary$proportionRiskyChoice_b2_upper <-
	(social_learning_model_validation_1022_riskID11_summary$raw_proportionRiskyChoice_b2_mean + social_learning_model_validation_1022_riskID11_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_1022_riskID11_summary$n)) %>% convert_alphaRaw_to_alpha

social_learning_model_validation_1022_riskID11_summary$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_1022_riskID11_summary$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha



# modest social learners' means
sigma_quantiles_1022_riskID11 = social_learning_model_validation_1022_riskID11_data %>% filter(condition_dummy==1) %>% .$soc_mean %>% quantile(probs = quantile_probs)
social_learning_model_validation_1022_riskID11_summary_sigma_level_1 <-
	social_learning_model_validation_1022_riskID11_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID11[1] & soc_mean < sigma_quantiles_1022_riskID11[3] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID11_summary_sigma_level_1$sigma_mean = mean(social_learning_model_validation_1022_riskID11_summary_sigma_level_1$soc_mean)

social_learning_model_validation_1022_riskID11_summary_sigma_level_2 <-
	social_learning_model_validation_1022_riskID11_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID11[3] & soc_mean < sigma_quantiles_1022_riskID11[5] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID11_summary_sigma_level_2$sigma_mean = mean(social_learning_model_validation_1022_riskID11_summary_sigma_level_2$soc_mean)

social_learning_model_validation_1022_riskID11_summary_sigma_level_3 <-
	social_learning_model_validation_1022_riskID11_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID11[5] & soc_mean < sigma_quantiles_1022_riskID11[7] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID11_summary_sigma_level_3$sigma_mean = mean(social_learning_model_validation_1022_riskID11_summary_sigma_level_3$soc_mean)

social_learning_model_validation_1022_riskID11_summary_sigma_level_4 <-
	social_learning_model_validation_1022_riskID11_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID11[7] & soc_mean < sigma_quantiles_1022_riskID11[9] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID11_summary_sigma_level_4$sigma_mean = mean(social_learning_model_validation_1022_riskID11_summary_sigma_level_4$soc_mean)

social_learning_model_validation_1022_riskID11_summary_sigma_level_5 <-
	social_learning_model_validation_1022_riskID11_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID11[9] & soc_mean < sigma_quantiles_1022_riskID11[11] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID11_summary_sigma_level_5$sigma_mean = mean(social_learning_model_validation_1022_riskID11_summary_sigma_level_5$soc_mean)

social_learning_model_validation_1022_riskID11_summary_sigma_merged = 
    social_learning_model_validation_1022_riskID11_summary_sigma_level_1 %>%
    rbind(social_learning_model_validation_1022_riskID11_summary_sigma_level_2) %>%
    rbind(social_learning_model_validation_1022_riskID11_summary_sigma_level_3) %>%
    rbind(social_learning_model_validation_1022_riskID11_summary_sigma_level_4) %>%
    rbind(social_learning_model_validation_1022_riskID11_summary_sigma_level_5)



# modest social learners' means
social_learning_model_validation_1022_riskID11_summary_reallyHighSigma <-
	social_learning_model_validation_1022_riskID11_data %>%
	dplyr::filter(soc_mean > 3/10 & soc_mean < 6/10 & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[2],
        proportionRiskyChoice_b2_upr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[8],
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_1022_riskID11_summary_reallyHighSigma$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_1022_riskID11_summary_reallyHighSigma$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha


# ======================================
# 2-risky 2-safe (4-armed) task
# ======================================

# fit result -- global parameters
fit_SL00_multiVar_LKJ_1022_globalparameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_globalparameters.csv'))
fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters.csv'))

## behavioural data summary
allBehaviour1022_group <- read_csv(paste0(main_exp_dir, "experimentalAnalysis/allBehaviour1022_group.csv"))
allBehaviour1022_group_riskID12 <- allBehaviour1022_group %>% dplyr::filter(riskDistributionId == 12)

fit_SL00_multiVar_LKJ_1022_parameters <- read_csv(paste0(main_exp_dir, "experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_parameters.csv"))

allBehaviour1022_group_riskID12_summarised_t35 <- allBehaviour1022_group_riskID12 %>%
	dplyr::filter(round>35) %>%
	group_by(amazonID, sub) %>%
	summarise(
		risky_choice_count = sum(best_risky_choice, na.rm = TRUE),
		risky_choice_mean = mean(best_risky_choice, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1]
	)

allBehaviour1022_group_riskID12_summarised_t35$groupID <- allBehaviour1022_group_riskID12_summarised_t35$room


allBehaviour1022_indiv <- read_csv(paste0(main_exp_dir, "experimentalAnalysis/allBehaviour1022_indiv.csv"))
allBehaviour1022_indiv_riskID12 <- allBehaviour1022_indiv %>% filter(riskDistributionId_factor=='Con: 2') #'Con: 2' means the 2-risky-2-safe task
allBehaviour1022_indiv_riskID12$sub_old <- allBehaviour1022_indiv_riskID12$sub
allBehaviour1022_indiv_riskID12$sub <- allBehaviour1022_indiv_riskID12$amazonID %>% as.factor() %>% as.numeric()

allBehaviour1022_indiv_riskID12_summarised_t35 <- allBehaviour1022_indiv_riskID12 %>%
	dplyr::filter(round>35) %>%
	group_by(amazonID, sub) %>%
	summarise(
		risky_choice_count = sum(best_risky_choice, na.rm = TRUE),
		risky_choice_mean = mean(best_risky_choice, na.rm=TRUE),
		trial_num = n(),
		indivOrGroup_factor = indivOrGroup_factor[1],
		room = room[1]
	)
allBehaviour1022_indiv_riskID12_summarised_t35$groupID = allBehaviour1022_indiv_riskID12_summarised_t35$room
allBehaviour1022_indiv_riskID12_summarised_t35$groupID[which(allBehaviour1022_indiv_riskID12_summarised_t35$indivOrGroup_factor=='Individual')] = 'Individual'

# Individual fits
fit_AL00_multiVar_LKJ_indiv_riskID12_parameters <- read_csv(paste0(main_exp_dir, 'experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID12_parameters.csv'))

# Merging the behavioural data with the fit parameters
# hot stove effect - individual
fit_AL_indiv_riskID12_parameters <- right_join(fit_AL00_multiVar_LKJ_indiv_riskID12_parameters, allBehaviour1022_indiv_riskID12_summarised_t35, by = 'sub')
fit_AL_indiv_riskID12_parameters$hot_stove_susceptibility <- fit_AL_indiv_riskID12_parameters$alpha_median_AL00_multiVar_LKJ * (1+ fit_AL_indiv_riskID12_parameters$beta_median_AL00_multiVar_LKJ)
fit_AL_indiv_riskID12_parameters$hot_stove_susceptibility_trancated <- fit_AL_indiv_riskID12_parameters$hot_stove_susceptibility
fit_AL_indiv_riskID12_parameters$hot_stove_susceptibility_trancated[which(fit_AL_indiv_riskID12_parameters$hot_stove_susceptibility > 6)] <- 6

# hot stove effect - group
fit_SL00_riskID12_parameters <- right_join(fit_SL00_multiVar_LKJ_1022_parameters, allBehaviour1022_group_riskID12_summarised_t35, by = 'sub')
fit_SL00_riskID12_parameters$hot_stove_susceptibility <- fit_SL00_riskID12_parameters$alpha_mean_SL00_multiVar_LKJ * (1+ fit_SL00_riskID12_parameters$beta_mean_SL00_multiVar_LKJ)
fit_SL00_riskID12_parameters$hot_stove_susceptibility_trancated <- fit_SL00_riskID12_parameters$hot_stove_susceptibility
fit_SL00_riskID12_parameters$hot_stove_susceptibility_trancated[which(fit_SL00_riskID12_parameters$hot_stove_susceptibility > 6)] <- 6

# overall means
social_learning_model_validation_1022_riskID12_summary <-
	social_learning_model_validation_1022_riskID12_data %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[2],
        proportionRiskyChoice_b2_upr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[8],
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% mean(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_1022_riskID12_summary$proportionRiskyChoice_b2_lower <-
	(social_learning_model_validation_1022_riskID12_summary$raw_proportionRiskyChoice_b2_mean - social_learning_model_validation_1022_riskID12_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_1022_riskID12_summary$n)) %>% convert_alphaRaw_to_alpha
social_learning_model_validation_1022_riskID12_summary$proportionRiskyChoice_b2_upper <-
	(social_learning_model_validation_1022_riskID12_summary$raw_proportionRiskyChoice_b2_mean + social_learning_model_validation_1022_riskID12_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_1022_riskID12_summary$n)) %>% convert_alphaRaw_to_alpha

social_learning_model_validation_1022_riskID12_summary$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_1022_riskID12_summary$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha


# modest social learners' means
sigma_quantiles_1022_riskID12 = social_learning_model_validation_1022_riskID12_data %>% filter(condition_dummy==1) %>% .$soc_mean %>% quantile(probs = quantile_probs)
social_learning_model_validation_1022_riskID12_summary_sigma_level_1 <-
	social_learning_model_validation_1022_riskID12_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID12[1] & soc_mean < sigma_quantiles_1022_riskID12[3] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID12_summary_sigma_level_1$sigma_mean = mean(social_learning_model_validation_1022_riskID12_summary_sigma_level_1$soc_mean)

social_learning_model_validation_1022_riskID12_summary_sigma_level_2 <-
	social_learning_model_validation_1022_riskID12_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID12[3] & soc_mean < sigma_quantiles_1022_riskID12[5] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID12_summary_sigma_level_2$sigma_mean = mean(social_learning_model_validation_1022_riskID12_summary_sigma_level_2$soc_mean)

social_learning_model_validation_1022_riskID12_summary_sigma_level_3 <-
	social_learning_model_validation_1022_riskID12_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID12[5] & soc_mean < sigma_quantiles_1022_riskID12[7] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID12_summary_sigma_level_3$sigma_mean = mean(social_learning_model_validation_1022_riskID12_summary_sigma_level_3$soc_mean)

social_learning_model_validation_1022_riskID12_summary_sigma_level_4 <-
	social_learning_model_validation_1022_riskID12_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID12[7] & soc_mean < sigma_quantiles_1022_riskID12[9] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID12_summary_sigma_level_4$sigma_mean = mean(social_learning_model_validation_1022_riskID12_summary_sigma_level_4$soc_mean)

social_learning_model_validation_1022_riskID12_summary_sigma_level_5 <-
	social_learning_model_validation_1022_riskID12_data %>%
	dplyr::filter(soc_mean > sigma_quantiles_1022_riskID12[9] & soc_mean < sigma_quantiles_1022_riskID12[11] & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mid = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median() %>% convert_alphaRaw_to_alpha,
        soc_mean = mean(soc_mean),
		n = n()
		)
social_learning_model_validation_1022_riskID12_summary_sigma_level_5$sigma_mean = mean(social_learning_model_validation_1022_riskID12_summary_sigma_level_5$soc_mean)

social_learning_model_validation_1022_riskID12_summary_sigma_merged = 
    social_learning_model_validation_1022_riskID12_summary_sigma_level_1 %>%
    rbind(social_learning_model_validation_1022_riskID12_summary_sigma_level_2) %>%
    rbind(social_learning_model_validation_1022_riskID12_summary_sigma_level_3) %>%
    rbind(social_learning_model_validation_1022_riskID12_summary_sigma_level_4) %>%
    rbind(social_learning_model_validation_1022_riskID12_summary_sigma_level_5)



# modest social learners' means
social_learning_model_validation_1022_riskID12_summary_reallyHighSigma <-
	social_learning_model_validation_1022_riskID12_data %>%
	dplyr::filter(soc_mean > 3/10 & soc_mean < 6/10 & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
        proportionRiskyChoice_b2_lwr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[2],
        proportionRiskyChoice_b2_upr = quantile(proportionRiskyChoice_b2, probs = quantile_probs)[8],
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_1022_riskID12_summary_reallyHighSigma$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_1022_riskID12_summary_reallyHighSigma$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha



###############################################################################################
##
##
## GLMM model fit 
##
## (1) Figure 6 now contains only the expetimental results with sample aberage curve 
## 
###############################################################################################
# stan model
# path -- mac mini
this_path <- '~/Documents/riskySocialLearning_airM1' # mac mini

# logit model 
stan_file_logit_regression <- file.path(this_path, 'logit_regression.stan')
stanmodel_logit_regression <- cmdstan_model(stan_file_logit_regression)
stanmodel_logit_regression$exe_file()

# debug
chains = 2
parallel_chains = 2
thin = 1
iter_warmup = 50
iter_sampling = 150
# # 
chains = 4
parallel_chains = 4
thin = 1
iter_warmup = 500
iter_sampling = 1000


# 0820 data
glmm_data_0820 = data.frame(
    amazonID = 1:length(c(parameterfit_indiv_AL00_0820$amazonID, fit_parameters_group_SL00_mcmc$amazonID)) %>% as.factor()
    , room = c(parameterfit_indiv_AL00_0820$room, fit_parameters_group_SL00_mcmc$room) %>% as.factor()
    , trial_num = c(parameterfit_indiv_AL00_0820$trial_num, fit_parameters_group_SL00_mcmc$trial_num)
    , risky_choice_count = c(parameterfit_indiv_AL00_0820$risky_choice_count, fit_parameters_group_SL00_mcmc$risky_choice_count)
    , risky_choice_mean = c(parameterfit_indiv_AL00_0820$risky_choice_mean, fit_parameters_group_SL00_mcmc$risky_choice_mean)
    , hot_stove_susceptibility_trancated = c(parameterfit_indiv_AL00_0820$hot_stove_susceptibility_trancated, fit_parameters_group_SL00_mcmc$hot_stove_susceptibility_trancated)
    , soc_median = c(rep(0, nrow(parameterfit_indiv_AL00_0820)), fit_parameters_group_SL00_mcmc$soc_median_SL00_multiVar_LKJ)
    , theta_median = c(rep(0, nrow(parameterfit_indiv_AL00_0820)), fit_parameters_group_SL00_mcmc$theta_median_SL00_multiVar_LKJ)
    , soc_mean = c(rep(0, nrow(parameterfit_indiv_AL00_0820)), fit_parameters_group_SL00_mcmc$soc_mean_SL00_multiVar_LKJ)
    , theta_mean = c(rep(0, nrow(parameterfit_indiv_AL00_0820)), fit_parameters_group_SL00_mcmc$theta_mean_SL00_multiVar_LKJ)
    , group_dummy = c(rep(0, nrow(parameterfit_indiv_AL00_0820)), rep(1, nrow(fit_parameters_group_SL00_mcmc)))
)

glmm_data_0820_stan = list(N = nrow(glmm_data_0820)
    , N_indiv = length(table(glmm_data_0820$amazonID))
    , N_group = length(table(glmm_data_0820$room))
    , risky_choice_count = glmm_data_0820$risky_choice_count
    , trial_num = glmm_data_0820$trial_num
    , hot_stove = glmm_data_0820$hot_stove_susceptibility_trancated
    , soc_mean = glmm_data_0820$soc_mean
    , theta_mean = glmm_data_0820$theta_mean
    , group_dummy = glmm_data_0820$group_dummy
    , Indiv = glmm_data_0820$amazonID %>% as.numeric()
    , Group = glmm_data_0820$room %>% as.numeric()
    , hot_stove_new_data = seq(0,6,length.out=201)
    , soc_mean_average = filter(glmm_data_0820, group_dummy==1)$soc_mean %>% mean()
    , theta_mean_average = filter(glmm_data_0820, group_dummy==1)$theta_mean %>% mean()
)

fit_logit_regression_0820 = stanmodel_logit_regression$sample(

  data = glmm_data_0820_stan
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
#   , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_logit_regression_0820_parameters = fit_logit_regression_0820$summary('beta')
write.csv(fit_logit_regression_0820_parameters
  , paste0(storing_directory , "fit_logit_regression_0820_parameters.csv"), row.names=FALSE)

fit_logit_regression_0820_prediction = fit_logit_regression_0820$summary(c('p_indiv','p_group'))
fit_logit_regression_0820_prediction$condition = rep(c('Individual','Group'),each=201)
fit_logit_regression_0820_prediction$hot_stove_susceptibility_trancated = rep(seq(0,6,length.out=201),2)
write.csv(fit_logit_regression_0820_prediction
  , paste0(storing_directory , "fit_logit_regression_0820_prediction.csv"), row.names=FALSE)



# #############################################################
# The followings are old analysis using glmer from lme4 package
# #############################################################

# pars = c("(Intercept)"
# , "hot_stove_susceptibility_trancated"
# , "group_dummy"
# , "soc_mean"
# , "theta_mean"
# , 'hot_stove_susceptibility_trancated:group_dummy'
# , "group_dummy:soc_mean:theta_mean"
# # , 'group_dummy:theta_mean'
# )

# # stan glmer
# poisson_fit_0820_stan = #glmer(cbind(risky_choice_count, trial_num-risky_choice_count) ~  
#     stan_glmer(cbind(risky_choice_count, trial_num - risky_choice_count) ~ 
#     hot_stove_susceptibility_trancated 
#     + group_dummy
#     + hot_stove_susceptibility_trancated * group_dummy 
#     + group_dummy * soc_mean 
#     + group_dummy * theta_mean
#     # + group_dummy * soc_mean  * theta_mean
#     + (1| amazonID ) + (1| room )
#     , data = glmm_data_0820, family = binomial(link = "logit"))
# summary(poisson_fit_0820_stan, pars = pars, digits = 2, prob=c(.025, .5, .975))
# bayesplot::color_scheme_set("brightblue") # 'viridis'
# plot(poisson_fit_0820_stan)


# glmm_data_0820 %>% 
#     group_by(group_dummy) %>%
#     data_grid(hot_stove_susceptibility_trancated = seq_range(hot_stove_susceptibility_trancated, n = 101)) %>%
#     add_epred_draws(risky_choice_count) %>%
#     ggplot(aes(x = hot_stove_susceptibility_trancated, y = risky_choice_count/trial_num, color = ordered(group_dummy))) +
#     stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), alpha = 0.25) +
#     geom_point(data = glmm_data_0820) +
#     scale_fill_brewer(palette = "Greys")

# # posterior prediction
# nd_indiv_0820 <- data.frame(hot_stove_susceptibility_trancated = seq(0, 6, , length.out = 200)
#     , trial_num = rep(35, 200)
#     , group_dummy = rep(0, 200)
#     , soc_mean = rep(0, 200)
#     , theta_mean = rep(0, 200)
#     , room = sample(glmm_data_0820$room, 200, replace = TRUE)
#     , amazonID = sample(glmm_data_0820$amazonID, 200, replace = TRUE)
#     )
# # PPD_indiv_0820 <- posterior_predict(poisson_fit_0820_stan, newdata = nd_indiv_0820)

# pred_0820 = tidybayes::add_predicted_draws(poisson_fit_0820_stan
#     , newdata = nd_indiv_0820
#     , allow_new_levels = T
# )

# glmm_data_0820 %>%
#   filter(group_dummy == 0) %>% 
#   add_fitted_draws(poisson_fit_0820_stan, n = 100) %>%
#   ggplot(aes(x = hot_stove_susceptibility_trancated, y = net)) +
#     geom_line(
#       aes(y = .value, group = paste(runner, .draw), color = runner),
#       alpha = 0.1) +
#     geom_point(aes(color = runner))


# # lme4 
# poisson_fit_0820 = glmer(cbind(risky_choice_count, trial_num-risky_choice_count) ~  
#     hot_stove_susceptibility_trancated 
#     + group_dummy
#     + hot_stove_susceptibility_trancated * group_dummy 
#     + group_dummy * soc_mean 
#     + group_dummy * theta_mean
#     # + group_dummy * soc_mean  * theta_mean
#     + (1| amazonID ) + (1| room )
#     , data = glmm_data_0820, family = binomial(link = "logit"))
# summary(poisson_fit_0820)
# # effect in the normal scale
# b0_poisson_fit_0820 <- fixef(poisson_fit_0820)['(Intercept)']
# b1_poisson_fit_0820 <- fixef(poisson_fit_0820)['hot_stove_susceptibility_trancated'] 
# b2_poisson_fit_0820 <- fixef(poisson_fit_0820)['group_dummy'] 
# b3_poisson_fit_0820 <- fixef(poisson_fit_0820)['soc_mean'] 
# b4_poisson_fit_0820 <- fixef(poisson_fit_0820)['theta_mean'] 
# b5_poisson_fit_0820 <- fixef(poisson_fit_0820)['hot_stove_susceptibility_trancated:group_dummy'] 
# # b5_poisson_fit_0820 <- fixef(poisson_fit_0820)['hot_stove_susceptibility_trancated:theta_mean']
# # b6_poisson_fit_0820 <- fixef(poisson_fit_0820)['soc_mean:theta_mean']
# # b7_poisson_fit_0820 <- fixef(poisson_fit_0820)['hot_stove_susceptibility_trancated:soc_mean:theta_mean']

# # the fit sigmoidal function
# poisson_prediction_0820 = function(hot_stove, isGroup, soc, theta) {
# 	1/(1+exp(-(b0_poisson_fit_0820 
#     + b1_poisson_fit_0820 * hot_stove 
#     + b2_poisson_fit_0820 * isGroup
#     + b3_poisson_fit_0820 * isGroup * soc
#     + b4_poisson_fit_0820 * isGroup * theta
#     + b5_poisson_fit_0820 * hot_stove * isGroup
#     # + b6_poisson_fit_0820 * soc * theta
#     # + b7_poisson_fit_0820 * hot_stove * soc* theta)
#     ))
# }

# # Confidence intervals -- individual
# new_data_0820_individual = data.frame(hot_stove_susceptibility_trancated = seq(0,4.3,0.01)
#     , group_dummy = rep(0, length(seq(0,4.3,0.01)))
#     , soc_mean = rep(0, length(seq(0,4.3,0.01)))
#     , theta_mean = rep(0, length(seq(0,4.3,0.01)))
#     , amazonID = rep(0, length(seq(0,4.3,0.01)))
#     , room = rep(0, length(seq(0,4.3,0.01)))
# )
# PI_0820_individual <- predictInterval(merMod = poisson_fit_0820, newdata = new_data_0820_individual,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_0820_individual <- as.data.frame(PI_0820_individual)
# PI_0820_individual <- PI_0820_individual %>% cbind(new_data_0820_individual)

# # Confidence intervals -- group
# new_data_0820_group = data.frame(hot_stove_susceptibility_trancated = seq(0,5.1,0.01)
#     , group_dummy = rep(1, length(seq(0,5.1,0.01)))
#     , soc_mean = rep(mean(fit_parameters_group_SL00_mcmc$soc_mean_SL00_multiVar_LKJ), length(seq(0,5.1,0.01)))
#     , theta_mean = rep(mean(fit_parameters_group_SL00_mcmc$theta_mean_SL00_multiVar_LKJ), length(seq(0,5.1,0.01)))
#     , amazonID = rep(0, length(seq(0,5.1,0.01)))
#     , room = rep(0, length(seq(0,5.1,0.01)))
#     # , amazonID = sample(nrow(parameterfit_indiv_AL00_0820):nrow(glmm_data_0820), size = length(seq(0,6,0.01)), replace = TRUE)
#     # , room = sample(nrow(parameterfit_indiv_AL00_0820):nrow(glmm_data_0820), size = length(seq(0,6,0.01)), replace = TRUE)
# )
# PI_0820_group <- predictInterval(merMod = poisson_fit_0820, newdata = new_data_0820_group,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_0820_group <- as.data.frame(PI_0820_group)
# PI_0820_group <- PI_0820_group %>% cbind(new_data_0820_group)

# GLMM prediction plot
(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.5,yend=0.5),colour="grey30", size=0.5) +
    geom_ribbon(data=fit_logit_regression_0820_prediction, mapping=aes(hot_stove_susceptibility_trancated, ymin=q5, ymax=q95, fill=condition), alpha=1/4)+
	geom_point(data = parameterfit_indiv_AL00_0820, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
	geom_point(data = fit_parameters_group_SL00_mcmc, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean), shape = 20, colour='orange') +
    geom_line(data=fit_logit_regression_0820_prediction, mapping=aes(hot_stove_susceptibility_trancated, mean, group=condition, color=condition), size=1.0)+
    scale_colour_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"), name='Condition')+
    scale_fill_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"), name='Condition')+
	# scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Arial()+
	xlim(c(0,6.5))+
	labs(
		# x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Positive risk premium\nThe 1-risky-1-safe task \n(N = 168)') +
	theme(legend.position = c(0.65, 0.7))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "red"))+
	NULL -> fig6_a
)
ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig6_a.png', plot = fig6_a, dpi = 600, width = 6, height = 4.5)


# 1120 riskID11 data
glmm_data_riskID11 = data.frame(
    amazonID = 1:length(c(fit_AL_indiv_riskID11_parameters$amazonID, fit_SL00_riskID11_parameters$amazonID)) %>% as.factor()
    , room = c(fit_AL_indiv_riskID11_parameters$room, fit_SL00_riskID11_parameters$room) %>% as.factor() 
    , trial_num = c(fit_AL_indiv_riskID11_parameters$trial_num, fit_SL00_riskID11_parameters$trial_num)
    , risky_choice_count = c(fit_AL_indiv_riskID11_parameters$risky_choice_count, fit_SL00_riskID11_parameters$risky_choice_count)
    , risky_choice_mean = c(fit_AL_indiv_riskID11_parameters$risky_choice_mean, fit_SL00_riskID11_parameters$risky_choice_mean)
    , hot_stove_susceptibility_trancated = c(fit_AL_indiv_riskID11_parameters$hot_stove_susceptibility_trancated, fit_SL00_riskID11_parameters$hot_stove_susceptibility_trancated)
    , soc_median = c(rep(0, nrow(fit_AL_indiv_riskID11_parameters)), fit_SL00_riskID11_parameters$soc_median_SL00_multiVar_LKJ)
    , theta_median = c(rep(0, nrow(fit_AL_indiv_riskID11_parameters)), fit_SL00_riskID11_parameters$theta_median_SL00_multiVar_LKJ)
    , soc_mean = c(rep(0, nrow(fit_AL_indiv_riskID11_parameters)), fit_SL00_riskID11_parameters$soc_mean_SL00_multiVar_LKJ)
    , theta_mean = c(rep(0, nrow(fit_AL_indiv_riskID11_parameters)), fit_SL00_riskID11_parameters$theta_mean_SL00_multiVar_LKJ)
    , group_dummy = c(rep(0, nrow(fit_AL_indiv_riskID11_parameters)), rep(1, nrow(fit_SL00_riskID11_parameters)))
)

glmm_data_riskID11_stan = list(N = nrow(glmm_data_riskID11)
    , N_indiv = length(table(glmm_data_riskID11$amazonID))
    , N_group = length(table(glmm_data_riskID11$room))
    , risky_choice_count = glmm_data_riskID11$risky_choice_count
    , trial_num = glmm_data_riskID11$trial_num
    , hot_stove = glmm_data_riskID11$hot_stove_susceptibility_trancated
    , soc_mean = glmm_data_riskID11$soc_mean
    , theta_mean = glmm_data_riskID11$theta_mean
    , group_dummy = glmm_data_riskID11$group_dummy
    , Indiv = glmm_data_riskID11$amazonID %>% as.numeric()
    , Group = glmm_data_riskID11$room %>% as.numeric()
    , hot_stove_new_data = seq(0,6,length.out=201)
    , soc_mean_average = filter(glmm_data_riskID11, group_dummy==1)$soc_mean %>% mean()
    , theta_mean_average = filter(glmm_data_riskID11, group_dummy==1)$theta_mean %>% mean()
)

fit_logit_regression_riskID11 = stanmodel_logit_regression$sample(

  data = glmm_data_riskID11_stan
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
#   , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_logit_regression_riskID11_parameters = fit_logit_regression_riskID11$summary('beta')
write.csv(fit_logit_regression_riskID11_parameters
  , paste0(storing_directory , "fit_logit_regression_riskID11_parameters.csv"), row.names=FALSE)

fit_logit_regression_riskID11_prediction = fit_logit_regression_riskID11$summary(c('p_indiv','p_group'))
fit_logit_regression_riskID11_prediction$condition = rep(c('Individual','Group'),each=201)
fit_logit_regression_riskID11_prediction$hot_stove_susceptibility_trancated = rep(seq(0,6,length.out=201),2)
write.csv(fit_logit_regression_riskID11_prediction
  , paste0(storing_directory , "fit_logit_regression_riskID11_prediction.csv"), row.names=FALSE)


# poisson_fit_riskID11 = glmer(cbind(risky_choice_count, trial_num-risky_choice_count) ~ 
#     hot_stove_susceptibility_trancated + soc_mean + theta_mean 
#     + hot_stove_susceptibility_trancated * soc_mean 
#     + hot_stove_susceptibility_trancated * theta_mean 
#     + soc_mean * theta_mean 
#     + hot_stove_susceptibility_trancated * soc_mean * theta_mean
#     + (1|amazonID) + (1|room)
#     , data = glmm_data_riskID11, family = "binomial")
# summary(poisson_fit_riskID11)

# # Confidence intervals -- individual
# new_data_riskID11_individual = data.frame(hot_stove_susceptibility_trancated = seq(0,5,0.01)
#     , soc_mean = rep(0, length(seq(0,5,0.01)))
#     , theta_mean = rep(0, length(seq(0,5,0.01)))
#     , amazonID = rep(0, length(seq(0,5,0.01)))
#     , room = rep(0, length(seq(0,5,0.01)))
# )
# PI_riskID11_individual <- predictInterval(merMod = poisson_fit_riskID11, newdata = new_data_riskID11_individual,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_riskID11_individual <- as.data.frame(PI_riskID11_individual)
# PI_riskID11_individual <- PI_riskID11_individual %>% cbind(new_data_riskID11_individual)

# # Confidence intervals -- group
# new_data_riskID11_group = data.frame(hot_stove_susceptibility_trancated = seq(0,6,0.01)
#     , soc_mean = rep(mean(fit_parameters_group_SL00_mcmc$soc_mean_SL00_multiVar_LKJ), length(seq(0,6,0.01)))
#     , theta_mean = rep(mean(fit_parameters_group_SL00_mcmc$theta_mean_SL00_multiVar_LKJ), length(seq(0,6,0.01)))
#     , amazonID = rep(0, length(seq(0,6,0.01)))
#     , room = rep(0, length(seq(0,6,0.01)))
# )
# PI_riskID11_group <- predictInterval(merMod = poisson_fit_riskID11, newdata = new_data_riskID11_group,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_riskID11_group <- as.data.frame(PI_riskID11_group)
# PI_riskID11_group <- PI_riskID11_group %>% cbind(new_data_riskID11_group)

# GLMM prediction plot

(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
    geom_ribbon(data=fit_logit_regression_riskID11_prediction, mapping=aes(hot_stove_susceptibility_trancated, ymin=q5, ymax=q95, fill=condition), alpha=1/4)+
	geom_point(data = fit_AL_indiv_riskID11_parameters, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
	geom_point(data = fit_SL00_riskID11_parameters, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean), colour='orange', shape = 20) +
    geom_line(data=fit_logit_regression_riskID11_prediction, mapping=aes(hot_stove_susceptibility_trancated, mean, group=condition, color=condition), size=1.0)+
    scale_colour_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"),name='Condition')+
    scale_fill_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"),name='Condition')+
	myTheme_Arial()+
	xlim(c(0,6.5))+
	labs(
		# x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Positive risk premium\nThe 1-risky-3-safe task \n(N = 148)') +
	theme(legend.position = c(0.65, 0.7))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "red"))+
	NULL -> fig6_b
)
ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig6_b.png', plot = fig6_b, dpi = 600, width = 6, height = 4.5)


# 1120 riskID12 data
glmm_data_riskID12 = data.frame(
    amazonID = 1:length(c(fit_AL_indiv_riskID12_parameters$amazonID, fit_SL00_riskID12_parameters$amazonID)) %>% as.factor()
    , room = c(fit_AL_indiv_riskID12_parameters$room, fit_SL00_riskID12_parameters$room) %>% as.factor() 
    , trial_num = c(fit_AL_indiv_riskID12_parameters$trial_num, fit_SL00_riskID12_parameters$trial_num)
    , risky_choice_count = c(fit_AL_indiv_riskID12_parameters$risky_choice_count, fit_SL00_riskID12_parameters$risky_choice_count)
    , risky_choice_mean = c(fit_AL_indiv_riskID12_parameters$risky_choice_mean, fit_SL00_riskID12_parameters$risky_choice_mean)
    , hot_stove_susceptibility_trancated = c(fit_AL_indiv_riskID12_parameters$hot_stove_susceptibility_trancated, fit_SL00_riskID12_parameters$hot_stove_susceptibility_trancated)
    , soc_median = c(rep(0, nrow(fit_AL_indiv_riskID12_parameters)), fit_SL00_riskID12_parameters$soc_median_SL00_multiVar_LKJ)
    , theta_median = c(rep(0, nrow(fit_AL_indiv_riskID12_parameters)), fit_SL00_riskID12_parameters$theta_median_SL00_multiVar_LKJ)
    , soc_mean = c(rep(0, nrow(fit_AL_indiv_riskID12_parameters)), fit_SL00_riskID12_parameters$soc_mean_SL00_multiVar_LKJ)
    , theta_mean = c(rep(0, nrow(fit_AL_indiv_riskID12_parameters)), fit_SL00_riskID12_parameters$theta_mean_SL00_multiVar_LKJ)
    , group_dummy = c(rep(0, nrow(fit_AL_indiv_riskID12_parameters)), rep(1, nrow(fit_SL00_riskID12_parameters)))
)


glmm_data_riskID12_stan = list(N = nrow(glmm_data_riskID12)
    , N_indiv = length(table(glmm_data_riskID12$amazonID))
    , N_group = length(table(glmm_data_riskID12$room))
    , risky_choice_count = glmm_data_riskID12$risky_choice_count
    , trial_num = glmm_data_riskID12$trial_num
    , hot_stove = glmm_data_riskID12$hot_stove_susceptibility_trancated
    , soc_mean = glmm_data_riskID12$soc_mean
    , theta_mean = glmm_data_riskID12$theta_mean
    , group_dummy = glmm_data_riskID12$group_dummy
    , Indiv = glmm_data_riskID12$amazonID %>% as.numeric()
    , Group = glmm_data_riskID12$room %>% as.numeric()
    , hot_stove_new_data = seq(0,6,length.out=201)
    , soc_mean_average = filter(glmm_data_riskID12, group_dummy==1)$soc_mean %>% mean()
    , theta_mean_average = filter(glmm_data_riskID12, group_dummy==1)$theta_mean %>% mean()
)

fit_logit_regression_riskID12 = stanmodel_logit_regression$sample(

  data = glmm_data_riskID12_stan
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
#   , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_logit_regression_riskID12_parameters = fit_logit_regression_riskID12$summary('beta')
write.csv(fit_logit_regression_riskID12_parameters
  , paste0(storing_directory , "fit_logit_regression_riskID12_parameters.csv"), row.names=FALSE)

fit_logit_regression_riskID12_prediction = fit_logit_regression_riskID12$summary(c('p_indiv','p_group'))
fit_logit_regression_riskID12_prediction$condition = rep(c('Individual','Group'),each=201)
fit_logit_regression_riskID12_prediction$hot_stove_susceptibility_trancated = rep(seq(0,6,length.out=201),2)
write.csv(fit_logit_regression_riskID12_prediction
  , paste0(storing_directory , "fit_logit_regression_riskID12_prediction.csv"), row.names=FALSE)


# poisson_fit_riskID12 = glmer(cbind(risky_choice_count, trial_num-risky_choice_count) ~ 
#     hot_stove_susceptibility_trancated + soc_mean + theta_mean 
#     + hot_stove_susceptibility_trancated * soc_mean 
#     + hot_stove_susceptibility_trancated * theta_mean 
#     + soc_mean * theta_mean 
#     + hot_stove_susceptibility_trancated * soc_mean * theta_mean
#     + (1|amazonID) + (1|room)
#     , data = glmm_data_riskID12, family = "binomial")
# summary(poisson_fit_riskID12)

# # Confidence intervals -- individual
# new_data_riskID12_individual = data.frame(hot_stove_susceptibility_trancated = seq(0,5,0.01)
#     , soc_mean = rep(0, length(seq(0,5,0.01)))
#     , theta_mean = rep(0, length(seq(0,5,0.01)))
#     , amazonID = rep(0, length(seq(0,5,0.01)))
#     , room = rep(0, length(seq(0,5,0.01)))
# )
# PI_riskID12_individual <- predictInterval(merMod = poisson_fit_riskID12, newdata = new_data_riskID12_individual,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_riskID12_individual <- as.data.frame(PI_riskID12_individual)
# PI_riskID12_individual <- PI_riskID12_individual %>% cbind(new_data_riskID12_individual)

# # Confidence intervals -- group
# new_data_riskID12_group = data.frame(hot_stove_susceptibility_trancated = seq(0,6,0.01)
#     , soc_mean = rep(mean(fit_parameters_group_SL00_mcmc$soc_mean_SL00_multiVar_LKJ), length(seq(0,6,0.01)))
#     , theta_mean = rep(mean(fit_parameters_group_SL00_mcmc$theta_mean_SL00_multiVar_LKJ), length(seq(0,6,0.01)))
#     , amazonID = rep(0, length(seq(0,6,0.01)))
#     , room = rep(0, length(seq(0,6,0.01)))
# )
# PI_riskID12_group <- predictInterval(merMod = poisson_fit_riskID12, newdata = new_data_riskID12_group,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_riskID12_group <- as.data.frame(PI_riskID12_group)
# PI_riskID12_group <- PI_riskID12_group %>% cbind(new_data_riskID12_group)

# GLMM prediction plot

(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
    geom_ribbon(data=fit_logit_regression_riskID12_prediction, mapping=aes(hot_stove_susceptibility_trancated, ymin=q5, ymax=q95, fill=condition), alpha=1/4)+
	geom_point(data = fit_AL_indiv_riskID12_parameters, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
	geom_point(data = fit_SL00_riskID12_parameters, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean), colour='orange', shape = 20) +
    geom_line(data=fit_logit_regression_riskID12_prediction, mapping=aes(hot_stove_susceptibility_trancated, mean, group=condition, color=condition), size=1.0)+
    scale_colour_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"),name='Condition')+
    scale_fill_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"),name='Condition')+
	myTheme_Arial()+
	xlim(c(0,6.5))+
	labs(
		# x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Positive risk premium\nThe 2-risky-2-safe task \n(N = 151)') +
	theme(legend.position = c(0.65, 0.7))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "red"))+
	NULL -> fig6_c
)
ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig6_c.png', plot = fig6_c, dpi = 600, width = 6, height = 4.5)



# 1120 revision_exp data
glmm_data_revision_exp = data.frame(
    amazonID = 1:length(c(parameterfit_indiv_AL00_revision_exp$amazonID, fit_parameters_revision_SL00_mcmc$amazonID)) %>% as.factor()
    , room = c(parameterfit_indiv_AL00_revision_exp$room, fit_parameters_revision_SL00_mcmc$room) %>% as.factor() 
    , trial_num = c(parameterfit_indiv_AL00_revision_exp$trial_num, fit_parameters_revision_SL00_mcmc$trial_num)
    , risky_choice_count = c(parameterfit_indiv_AL00_revision_exp$risky_choice_count, fit_parameters_revision_SL00_mcmc$risky_choice_count)
    , risky_choice_mean = c(parameterfit_indiv_AL00_revision_exp$risky_choice_mean, fit_parameters_revision_SL00_mcmc$risky_choice_mean)
    , hot_stove_susceptibility_trancated = c(parameterfit_indiv_AL00_revision_exp$hot_stove_susceptibility_trancated, fit_parameters_revision_SL00_mcmc$hot_stove_susceptibility_trancated)
    , soc_median = c(rep(0, nrow(parameterfit_indiv_AL00_revision_exp)), fit_parameters_revision_SL00_mcmc$soc_median_SL00_multiVar_LKJ)
    , theta_median = c(rep(0, nrow(parameterfit_indiv_AL00_revision_exp)), fit_parameters_revision_SL00_mcmc$theta_median_SL00_multiVar_LKJ)
    , soc_mean = c(rep(0, nrow(parameterfit_indiv_AL00_revision_exp)), fit_parameters_revision_SL00_mcmc$soc_mean_SL00_multiVar_LKJ)
    , theta_mean = c(rep(0, nrow(parameterfit_indiv_AL00_revision_exp)), fit_parameters_revision_SL00_mcmc$theta_mean_SL00_multiVar_LKJ)
    , group_dummy = c(rep(0, nrow(parameterfit_indiv_AL00_revision_exp)), rep(1, nrow(fit_parameters_revision_SL00_mcmc)))
)


glmm_data_revision_exp_stan = list(N = nrow(glmm_data_revision_exp)
    , N_indiv = length(table(glmm_data_revision_exp$amazonID))
    , N_group = length(table(glmm_data_revision_exp$room))
    , risky_choice_count = glmm_data_revision_exp$risky_choice_count
    , trial_num = glmm_data_revision_exp$trial_num
    , hot_stove = glmm_data_revision_exp$hot_stove_susceptibility_trancated
    , soc_mean = glmm_data_revision_exp$soc_mean
    , theta_mean = glmm_data_revision_exp$theta_mean
    , group_dummy = glmm_data_revision_exp$group_dummy
    , Indiv = glmm_data_revision_exp$amazonID %>% as.numeric()
    , Group = glmm_data_revision_exp$room %>% as.numeric()
    , hot_stove_new_data = seq(0,6,length.out=201)
    , soc_mean_average = filter(glmm_data_revision_exp, group_dummy==1)$soc_mean %>% mean()
    , theta_mean_average = filter(glmm_data_revision_exp, group_dummy==1)$theta_mean %>% mean()
)

fit_logit_regression_revision_exp = stanmodel_logit_regression$sample(

  data = glmm_data_revision_exp_stan
  , seed = 777 #output_dir=".", validate_csv = FALSE, # cmdstanr's bug...
#   , init = 2
  , adapt_delta = 0.80 # maybe 0.9 default 0.8
  , chains = chains, parallel_chains = parallel_chains
  , thin = thin, iter_warmup = iter_warmup, iter_sampling = iter_sampling 

)
fit_logit_regression_revision_exp_parameters = fit_logit_regression_revision_exp$summary('beta')
write.csv(fit_logit_regression_revision_exp_parameters
  , paste0(storing_directory , "fit_logit_regression_revision_exp_parameters.csv"), row.names=FALSE)

fit_logit_regression_revision_exp_prediction = fit_logit_regression_revision_exp$summary(c('p_indiv','p_group'))
fit_logit_regression_revision_exp_prediction$condition = rep(c('Individual','Group'),each=201)
fit_logit_regression_revision_exp_prediction$hot_stove_susceptibility_trancated = rep(seq(0,6,length.out=201),2)
write.csv(fit_logit_regression_revision_exp_prediction
  , paste0(storing_directory , "fit_logit_regression_revision_exp_prediction.csv"), row.names=FALSE)



# poisson_fit_revision_exp = glmer(cbind(risky_choice_count, trial_num-risky_choice_count) ~ 
#     hot_stove_susceptibility_trancated + soc_mean + theta_mean 
#     + hot_stove_susceptibility_trancated * soc_mean 
#     + hot_stove_susceptibility_trancated * theta_mean 
#     + soc_mean * theta_mean 
#     + hot_stove_susceptibility_trancated * soc_mean * theta_mean
#     + (1|amazonID) + (1|room)
#     , data = glmm_data_revision_exp, family = "binomial")
# summary(poisson_fit_revision_exp)

# # Confidence intervals -- individual
# new_data_revision_exp_individual = data.frame(hot_stove_susceptibility_trancated = seq(0,3,0.01)
#     , soc_mean = rep(0, length(seq(0,3,0.01)))
#     , theta_mean = rep(0, length(seq(0,3,0.01)))
#     , amazonID = rep(0, length(seq(0,3,0.01)))
#     , room = rep(0, length(seq(0,3,0.01)))
# )
# PI_revision_exp_individual <- predictInterval(merMod = poisson_fit_revision_exp, newdata = new_data_revision_exp_individual,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_revision_exp_individual <- as.data.frame(PI_revision_exp_individual)
# PI_revision_exp_individual <- PI_revision_exp_individual %>% cbind(new_data_revision_exp_individual)

# # Confidence intervals -- group
# new_data_revision_exp_group = data.frame(hot_stove_susceptibility_trancated = seq(0,6,0.01)
#     , soc_mean = rep(mean(fit_parameters_group_SL00_mcmc$soc_mean_SL00_multiVar_LKJ), length(seq(0,6,0.01)))
#     , theta_mean = rep(mean(fit_parameters_group_SL00_mcmc$theta_mean_SL00_multiVar_LKJ), length(seq(0,6,0.01)))
#     , amazonID = rep(0, length(seq(0,6,0.01)))
#     , room = rep(0, length(seq(0,6,0.01)))
# )
# PI_revision_exp_group <- predictInterval(merMod = poisson_fit_revision_exp, newdata = new_data_revision_exp_group,
#                     level = 0.80, n.sims = 10000,
#                     stat = "median", type="probability",
#                     include.resid.var = TRUE)
# PI_revision_exp_group <- as.data.frame(PI_revision_exp_group)
# PI_revision_exp_group <- PI_revision_exp_group %>% cbind(new_data_revision_exp_group)

# GLMM prediction plot

(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.5,yend=0.5),colour="grey30", size=0.5) +
    geom_ribbon(data=fit_logit_regression_revision_exp_prediction, mapping=aes(hot_stove_susceptibility_trancated, ymin=q5, ymax=q95, fill=condition), alpha=1/4)+
	geom_point(data = parameterfit_indiv_AL00_revision_exp, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
	geom_point(data = fit_parameters_revision_SL00_mcmc, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean), colour='orange', shape = 20) +
    geom_line(data=fit_logit_regression_revision_exp_prediction, mapping=aes(hot_stove_susceptibility_trancated, mean, group=condition, color=condition), size=1.0)+
    scale_colour_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"),name='Condition')+
    scale_fill_manual(values = c("orange", "grey20"), breaks = c("Group", "Individual"),name='Condition')+
	myTheme_Arial()+
	xlim(c(0,6.5))+
    ylim(c(0,1))+
	labs(
		# x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Negative risk premium\nThe 1-risky-1-safe task\n(N = 118)') +
	theme(legend.position = c(0.25, 0.75))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "blue"))+
	NULL -> fig6_d
)
ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig6_d.png', plot = fig6_d, dpi = 600, width = 6, height = 4.5)


figure_6_2nd_rev <- plot_grid(fig6_a, fig6_b, fig6_c, fig6_d
    # ,  common.legend = TRUE
    # ,  legend = 'right'
    , axis = 'bt'
    , labels = c('a','b','c','d'), ncol = 2, align = 'v'
)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/figure_6_2nd_rev.png"
    , plot = figure_6_2nd_rev, dpi = 600, width = 12, height = 8
	# , device = cairo_pdf
	# , device = 'pdf'
	)
# ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/figure_6_2nd_rev.pdf"
#     , plot = figure_6_2nd_rev, dpi = 600, width = 12, height = 8
# 	# , device = cairo_pdf
# 	# , device = 'pdf'
# 	)




###############################################################################################
##
##
## Copmutational model fit
##
## (2) A new Figure 7 shows the model prediction curves and the histograms 
## 
###############################################################################################

(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.5,yend=0.5),colour="grey30", size=0.5) +
	geom_ribbon(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==0 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='grey20', alpha=1/2)+
	geom_ribbon(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==1 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='orange', alpha=1/2)+
    # geom_ribbon(data=social_learning_model_validation_0820_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='purple', alpha=1/2)+
	# geom_point(data = parameterfit_indiv_AL00_0820, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
	# geom_point(data = fit_parameters_group_SL00_mcmc, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean, colour=soc_mean), shape = 20) +
	geom_line(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid), size=1.0, linetype='dashed')+
	# geom_line(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==1), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), size=1.0)+
    geom_line(data=social_learning_model_validation_0820_summary_sigma_merged, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=sigma_mean, color=sigma_mean), linetype='solid', size =1.0) +
	# geom_line(data=social_learning_model_validation_0820_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), linetype = 'dashed', size=1.0)+
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Arial()+
	xlim(c(0,6.5))+
    ylim(c(0,1))+
	labs(
		# x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Positive risk premium\nThe 1-risky-1-safe task') +
	theme(legend.position = c(0.5, 0.77))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "red"))+
	NULL -> fig7_a
)

ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig7_a.png', plot = fig7_a, dpi = 600, width = 6, height = 4.5)

(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
	geom_ribbon(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==0 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='grey20', alpha=1/2)+
	geom_ribbon(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==1 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='orange', alpha=1/2)+
	geom_line(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid), size=1.0, linetype='dashed')+
    geom_line(data=social_learning_model_validation_1022_riskID11_summary_sigma_merged, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=sigma_mean, color=sigma_mean), linetype='solid', size =1.0) +
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Arial()+
	xlim(c(0,6.5))+
    ylim(c(0,1))+
	labs(
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Positive risk premium\nThe 1-risky-3-safe task') +
	theme(legend.position = c(0.65, 0.6))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "red"))+
	NULL -> fig7_b
)

ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig7_b.png', plot = fig7_b, dpi = 600, width = 6, height = 4.5)

(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
	geom_ribbon(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==0 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='grey20', alpha=1/2)+
	geom_ribbon(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==1 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='orange', alpha=1/2)+
	geom_line(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid), size=1.0, linetype='dashed')+
    geom_line(data=social_learning_model_validation_1022_riskID12_summary_sigma_merged, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=sigma_mean, color=sigma_mean), linetype='solid', size =1.0) +
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Arial()+
	xlim(c(0,6.5))+
    ylim(c(0,1))+
	labs(
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Positive risk premium\nThe 2-risky-2-safe task') +
	theme(legend.position = c(0.65, 0.6))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "red"))+
	NULL -> fig7_c
)

ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig7_c.png', plot = fig7_c, dpi = 600, width = 6, height = 4.5)


(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.5,yend=0.5),colour="grey30", size=0.5) +
	geom_ribbon(data=social_learning_model_validation_revision_exp_summary%>%dplyr::filter(condition_dummy==0 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='grey20', alpha=1/2)+
	geom_ribbon(data=social_learning_model_validation_revision_exp_summary%>%dplyr::filter(condition_dummy==1 & hot_stove_susceptibility_rounded<6), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lwr, ymax=proportionRiskyChoice_b2_upr), fill='orange', alpha=1/2)+
	geom_line(data=social_learning_model_validation_revision_exp_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid), size=1.0, linetype='dashed')+
    geom_line(data=social_learning_model_validation_revision_exp_summary_sigma_merged, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=sigma_mean, color=sigma_mean), linetype='solid', size =1.0) +
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Arial()+
	xlim(c(0,6.5))+
    ylim(c(0,1))+
	labs(
		x = expression(paste('Susceptibility to the hot stove effect ', alpha[i], '*(', beta[i], '+1)')),
		y = 'Mean proportion of choosing\nthe risky option',
		title = 'Negative risk premium\nThe 1-risky-1-safe task') +
	theme(legend.position = c(0.65, 0.6))+
	# theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12, face='bold'))+
	theme(legend.text = element_text(size=11))+
    theme(plot.title = element_text(colour = "blue"))+
	NULL -> fig7_d
)

ggsave(file = '~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/fig7_d.png', plot = fig7_d, dpi = 600, width = 6, height = 4.5)



figure_7_2nd_rev <- plot_grid(fig7_a, fig7_b, fig7_c, fig7_d
    # ,  common.legend = TRUE
    # ,  legend = 'right'
    , axis = 'bt'
    , labels = c('a','b','c','d'), ncol = 2, align = 'v'
)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/exp_reanalysis_result/figure_7_2nd_rev.png"
    , plot = figure_7_2nd_rev, dpi = 600, width = 12, height = 8
	# , device = cairo_pdf
	# , device = 'pdf'
	)