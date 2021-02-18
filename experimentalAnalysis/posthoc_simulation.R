###############################################################################################
##
## Experimental results - 22/10/2020 (four-armed bandit task)
##
## Wataru Toyokawa
## 17 Feb 2021
###############################################################################################
if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace

# Loading
library(tidyverse)
library(cowplot)
library(ggpubr)
library(foreach)
library(MASS)
library(doParallel)
library(lme4)
#registerDoParallel(detectCores()) # this uses as many core as available
registerDoParallel(detectCores())

setwd("~/Dropbox/wataru/papers/RiskySocialLearning/experiment/overall_analysis")
source("functions.R")

# number of the simulation run
repetition = 100000 # 1000

# =============================================
# The 1-risky 3-safe task
# =============================================

# fit result -- global parameters
fit_SL00_multiVar_LKJ_1022_globalparameters <- read.csv('fit_SL00_multiVar_LKJ_1022_globalparameters.csv')
fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters <- read.csv('fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters.csv')

## behavioural data summary
allBehaviour1022_group <- read.csv("allBehaviour1022_group.csv")
allBehaviour1022_group_riskID11 <- allBehaviour1022_group %>% dplyr::filter(riskDistributionId == 11)

fit_SL00_multiVar_LKJ_1022_parameters <- read.csv("fit_SL00_multiVar_LKJ_1022_parameters.csv")

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


allBehaviour1022_indiv <- read.csv("allBehaviour1022_indiv.csv")
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
fit_AL00_multiVar_LKJ_indiv_riskID11_parameters <- read.csv('fit_AL00_multiVar_LKJ_indiv_riskID11_parameters.csv')

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


# -- global parameters of the group condition
mu_alpha_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable == 'mu_alpha[1]')]
mu_beta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='mu_beta[1]')]
mu_soc_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable == 'mu_soc0[1]')]
mu_theta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable == 'mu_theta[1]')]
s_alpha_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_alpha[1]')]
s_beta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_beta[1]')]
s_soc_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_soc0[1]')]
s_theta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_theta[1]')]

# -- global parameters for the individual condition
mu_alpha_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$variable=='mu_alpha')]
mu_beta_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$variable=='mu_beta')]
s_alpha_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$variable=='s_alpha')]
s_beta_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters$variable=='s_beta')]

# -- task's setting
# repetition = 100000 # 1000
length_of_the_last_block = 35

## the task: 4-armed bandit problem with the setting 0 (riskID11)
riskPremium = 20/15
rSure = 1
payoff_L1 = 1.5
payoff_H1 = 1.5
payoff_L2 = 1.25
payoff_H2 = 1.25
payoff_L3 = 1.0
payoff_H3 = 1.0
# risky option
rRisky = 0.4
payoff_L4 = 0.50
payoff_H4 = (riskPremium*1.5 - payoff_L4*(1-rRisky))/rRisky
initialExpextation = 0

# simulation's setting
horizon = 70 # = number of trials
numOptions = 4
length_of_the_last_block = 30
# simulated model
model = 'Decision-Biasing'

# -- Simulation
social_learning_model_validation_1022_riskID11 = list()
conditionList = c('Individual', 'Group')
groupSize = 8

s_time = Sys.time()
for(condition in conditionList) {
	social_learning_model_validation_1022_riskID11[[paste("n=", groupSize)]][[paste("condition=", condition)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
		## Initial settings
		choices = matrix(nrow=groupSize, ncol=horizon)
		payoffs = matrix(nrow=groupSize, ncol=horizon)
		performance = matrix(nrow=groupSize, ncol=horizon)
		safeChoiceProb = matrix(nrow=groupSize, ncol=horizon)
		isThisBestOption = matrix(nrow=groupSize, ncol=horizon)
		optimalChoiceProb = matrix(nrow=groupSize, ncol=horizon)
		expectedPerformance = matrix(nrow=groupSize, ncol=horizon)
		Q = array(dim = c(numOptions, horizon, groupSize))
		choiceCounter = array(1, dim = c(numOptions, groupSize))
		netChoiceProb = array(dim = c(numOptions, horizon, groupSize))
		netChoiceProb[,1,] = 1/numOptions
		Q[,1,] = initialExpextation
		socialFrequency = matrix(nrow=numOptions, ncol=horizon)
		socialFrequency[,] = 1e-1
		## Setting individual parameters
		if(condition == 'Group') {
			thisGroup <- sample(1:nrow(fit_SL00_riskID11_parameters), size=groupSize, replace=TRUE)
			# individual parameters are drawn from the fit global parameters
			thisAlpha <- (mu_alpha_group + s_alpha_group * rnorm(groupSize, 0, 1)) %>%
			  convert_alphaRaw_to_alpha()
			thisBeta <- (mu_beta_group + s_beta_group * rnorm(groupSize, 0, 1)) %>% exp()
			thisSigma <- (mu_soc_group + s_soc_group * rnorm(groupSize, 0, 1)) %>%
			  convert_alphaRaw_to_alpha()
			thisTheta <- mu_theta_group + s_theta_group * rnorm(groupSize, 0, 1)
			# # ----- if individual fit parameters are used directly ------
			# thisAlpha <- fit_SL00_riskID11_parameters$alpha_mean[thisGroup]
			# thisBeta <- fit_SL00_riskID11_parameters$beta_mean[thisGroup]
			# thisSigma <- fit_SL00_riskID11_parameters$soc_mean[thisGroup]
			# thisTheta <- fit_SL00_riskID11_parameters$theta_mean[thisGroup]
		} else { # if condition == individual
			thisAlpha <- (mu_alpha_indiv + s_alpha_indiv * rnorm(groupSize, 0, 1)) %>%
			  convert_alphaRaw_to_alpha()
			thisBeta <- (mu_beta_indiv + s_beta_indiv * rnorm(groupSize, 0, 1)) %>% exp()
			thisTheta <- rep(0, groupSize)
			thisSigma <- rep(0, groupSize)
		}
		## running the task
		for(t in 1:horizon){
			# each individual chooses one option based on his/her choice probability
			# choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
			choices[,t] = mapply(function(p1,p2,p3,p4){ sample(1:numOptions, 1, prob=c(p1,p2,p3,p4), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,], netChoiceProb[3,t,], netChoiceProb[4,t,] )
			# each subject earns some money (if lucky)
			# payoffs[,t] = payoffGenerateBinary(groupSize, choices[,t], rSure, rRisky, payoff_sureL, payoff_sureH, payoff_risky1, payoff_risky2)
			payoffs[,t] = payoffGenerate4Arm_unsync(groupSize, choices[,t], rSure, rSure, rSure, rRisky, payoff_H1, payoff_L1, payoff_H2, payoff_L2, payoff_H3, payoff_L3, payoff_H4, payoff_L4)
			# update choiceCounter and learningRate (if the learning rate is an averaging rule in this simulation.)
			updatingPositions = (choices[,t] + numOptions*(1:groupSize-1))
			
			# value updating
			if(t < horizon) {
				if(t == 1) {
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					# In the first trial, all Q values are updated by the first experience
					QQ[,t+1] = QQ[,t] + thisAlpha * (payoffs[,t] - QQ[,t])
					dim(QQ) = c(numOptions, groupSize, horizon)
					Q = aperm(QQ, c(1,3,2))
				} else {
					# Updating Q value based on Rescorla-Wagner model (Weighted return model)
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					#QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + learningRate[updatingPositions] * (payoffs[,t] - QQ[updatingPositions,t])
					QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + thisAlpha * (payoffs[,t] - QQ[updatingPositions,t])
					dim(QQ) = c(numOptions, groupSize, horizon)
					Q = aperm(QQ, c(1,3,2))
				}

				# update socialFrequency
				## Option's frequency
				for(i in 1:numOptions){
					if(length(which(names(table(choices[,t]))==i))>0) {
						socialFrequency[i,t+1] = socialFrequency[i,t+1] + table(choices[,t])[which(names(table(choices[,t]))==i)][1]
					}
				}
				
				###############
				## Softmax choice base solely on Q values
				###############
				Q_exp = ( Q[,t+1,] * rep(thisBeta, each = numOptions) ) %>% apply(2,expCeiling)
				softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
				freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
				## The followings update the choice probability matrix
				###############
				## Softmax -- END
				###############

				if(model=='Decision-Biasing'){
  				netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-thisSigma)) %>% t() +
  				  apply(freqDepenMatrix, 1, multiplyBeta, beta=thisSigma) %>% t()
				}else{
				  netMatrix = softmaxMatrix
				}
				netChoiceProbAperm = aperm(netChoiceProb, c(1,3,2))
				dim(netChoiceProbAperm) = c(numOptions*groupSize, horizon)
				dim(netMatrix) = c(numOptions*groupSize, 1)
				netChoiceProbAperm[,t+1] = netMatrix
				dim(netChoiceProbAperm) = c(numOptions, groupSize, horizon)
				netChoiceProb = aperm(netChoiceProbAperm, c(1,3,2))
			}
		}

		for(i in 1:groupSize) {
			safeChoiceProb[i,] = 1 - netChoiceProb[4,,i]
		}

		safeChoiceProbMean_b1 = safeChoiceProb[,1:length_of_the_last_block] %>% apply(1, mean) #%>% mean()
		safeChoiceProbMean_b2 = safeChoiceProb[,(length_of_the_last_block+1):horizon] %>% apply(1, mean) #%>% mean()

		choiceProb_option1_b1 = netChoiceProb[1, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option1_b2 = netChoiceProb[1,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()
		choiceProb_option2_b1 = netChoiceProb[2, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option2_b2 = netChoiceProb[2,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()
		choiceProb_option3_b1 = netChoiceProb[3, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option3_b2 = netChoiceProb[3,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()
		choiceProb_option4_b1 = netChoiceProb[4, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option4_b2 = netChoiceProb[4,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()

		# Submitting this repetition's result
		condition_dummy <- 0
		if(condition=='Group') condition_dummy <- 1
		print(
			data.frame(condition_dummy = rep(condition_dummy, groupSize)
					, hot_stove_susceptibility = thisAlpha * (thisBeta + 1)
					, alpha_mean = thisAlpha
					, beta_mean = thisBeta
					, soc_mean = thisSigma
					, theta_mean = thisTheta
					, proportionSafeChoice_b1 = safeChoiceProbMean_b1
					, proportionSafeChoice_b2 = safeChoiceProbMean_b2
					, proportionRiskyChoice_b1 = choiceProb_option4_b1
					, proportionRiskyChoice_b2 = choiceProb_option4_b2
				)
		)
	}
gc();gc() # rubbish collection
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
social_learning_model_validation_1022_riskID11_data <- social_learning_model_validation_1022_riskID11[[paste("n=", groupSize)]][[paste("condition=", "Individual")]] %>% data.frame()

social_learning_model_validation_1022_riskID11_data <- social_learning_model_validation_1022_riskID11_data %>% rbind(social_learning_model_validation_1022_riskID11[[paste("n=", groupSize)]][[paste("condition=", "Group")]] %>% data.frame())

social_learning_model_validation_1022_riskID11_data$hot_stove_susceptibility_rounded <- (social_learning_model_validation_1022_riskID11_data$hot_stove_susceptibility * 5) %>% round()/5
social_learning_model_validation_1022_riskID11_data$hot_stove_susceptibility_rounded[which(social_learning_model_validation_1022_riskID11_data$hot_stove_susceptibility_rounded>6)] <- 6

social_learning_model_validation_1022_riskID11_data$soc_mean_category <- 'mild'

# overall means
social_learning_model_validation_1022_riskID11_summary <-
	social_learning_model_validation_1022_riskID11_data %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
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
social_learning_model_validation_1022_riskID11_summary_reallyHighSigma <-
	social_learning_model_validation_1022_riskID11_data %>%
	dplyr::filter(soc_mean > 5/10 & soc_mean < 8/10 & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_1022_riskID11_summary_reallyHighSigma$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_1022_riskID11_summary_reallyHighSigma$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha


(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
	geom_point(data = fit_SL00_riskID11_parameters, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean, colour=soc_mean_SL00_multiVar_LKJ), shape = 20) +
	geom_point(data = fit_AL_indiv_riskID11_parameters, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 18)+ # shape=5: diamond
	geom_line(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid))+
	geom_line(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==1), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)))+
	geom_line(data=social_learning_model_validation_1022_riskID11_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), linetype = 'dashed')+
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Helvetica()+
	xlim(c(0,6.5))+
	labs(
		x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '*(', beta[i], '+1)'))),
		y = 'Mean proportion of choosing\nthe optimal risky option',
		title = 'The 1-risky-3-safe task \n(N = 148)') +
	theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12))+
	theme(legend.text = element_text(size=11))+
	NULL -> fig6_b)

write.csv(social_learning_model_validation_1022_riskID11_data,
      "social_learning_model_validation_1022_riskID11_data.csv", row.names=FALSE)




# =============================================
# The 2-risky 2-safe task
# =============================================

# fit result -- global parameters
fit_SL00_multiVar_LKJ_1022_globalparameters <- read.csv('fit_SL00_multiVar_LKJ_1022_globalparameters.csv')
fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters <- read.csv('fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters.csv')

## behavioural data summary
allBehaviour1022_group <- read.csv("allBehaviour1022_group.csv")
allBehaviour1022_group_riskID12 <- allBehaviour1022_group %>% dplyr::filter(riskDistributionId == 12)

fit_SL00_multiVar_LKJ_1022_parameters <- read.csv("fit_SL00_multiVar_LKJ_1022_parameters.csv")

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


allBehaviour1022_indiv <- read.csv("allBehaviour1022_indiv.csv")
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
fit_AL00_multiVar_LKJ_indiv_riskID12_parameters <- read.csv('fit_AL00_multiVar_LKJ_indiv_riskID12_parameters.csv')

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


# -- global parameters of the group condition
mu_alpha_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable == 'mu_alpha[2]')]
mu_beta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='mu_beta[2]')]
mu_soc_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable == 'mu_soc0[2]')]
mu_theta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable == 'mu_theta[2]')]
s_alpha_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_alpha[2]')]
s_beta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_beta[2]')]
s_soc_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_soc0[2]')]
s_theta_group <- fit_SL00_multiVar_LKJ_1022_globalparameters$mean[which(fit_SL00_multiVar_LKJ_1022_globalparameters$variable=='s_theta[2]')]

# -- global parameters for the individual condition
mu_alpha_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$variable=='mu_alpha')]
mu_beta_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$variable=='mu_beta')]
s_alpha_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$variable=='s_alpha')]
s_beta_indiv <- fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters$variable=='s_beta')]

# -- task's setting
# repetition = 100000 # 1000
length_of_the_last_block = 35

## the task: 4-armed bandit problem with the setting 0 (riskID12)
riskPremium = 20/15
rSure = 1
payoff_L1 = 1.5
payoff_H1 = 1.5
payoff_L2 = 1.25
payoff_H2 = 1.25
# risky option
rRisky = 0.4
payoff_L3 = 0.5 #decoy
payoff_H3 = (1*payoff_L2 - payoff_L3*(1-rRisky))/rRisky
payoff_L4 = 0.5
payoff_H4 = (riskPremium*payoff_L1 - payoff_L4*(1-rRisky))/rRisky
initialExpextation = 0

## the task: 4-armed bandit problem with the setting 0 (riskID12)
# riskPremium = 20/15
# rSure = 1
# payoff_L1 = 1.5
# payoff_H1 = 1.5
# payoff_L2 = 1.25
# payoff_H2 = 1.25
# payoff_L3 = 1.0
# payoff_H3 = 1.0
# # risky option
# rRisky = 0.4
# payoff_L4 = 0.50
# payoff_H4 = (riskPremium*1.5 - payoff_L4*(1-rRisky))/rRisky
# initialExpextation = 0

# simulation's setting
horizon = 70 # = number of trials
numOptions = 4
length_of_the_last_block = 30
# simulated model
model = 'Decision-Biasing'

# -- Simulation
social_learning_model_validation_1022_riskID12 = list()
conditionList = c('Individual', 'Group')
groupSize = 8

s_time = Sys.time()
for(condition in conditionList) {
	social_learning_model_validation_1022_riskID12[[paste("n=", groupSize)]][[paste("condition=", condition)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
		## Initial settings
		choices = matrix(nrow=groupSize, ncol=horizon)
		payoffs = matrix(nrow=groupSize, ncol=horizon)
		performance = matrix(nrow=groupSize, ncol=horizon)
		safeChoiceProb = matrix(nrow=groupSize, ncol=horizon)
		isThisBestOption = matrix(nrow=groupSize, ncol=horizon)
		optimalChoiceProb = matrix(nrow=groupSize, ncol=horizon)
		expectedPerformance = matrix(nrow=groupSize, ncol=horizon)
		Q = array(dim = c(numOptions, horizon, groupSize))
		choiceCounter = array(1, dim = c(numOptions, groupSize))
		netChoiceProb = array(dim = c(numOptions, horizon, groupSize))
		netChoiceProb[,1,] = 1/numOptions
		Q[,1,] = initialExpextation
		socialFrequency = matrix(nrow=numOptions, ncol=horizon)
		socialFrequency[,] = 1e-1
		## Setting individual parameters
		if(condition == 'Group') {
			thisGroup <- sample(1:nrow(fit_SL00_riskID12_parameters), size=groupSize, replace=TRUE)
			# individual parameters are drawn from the fit global parameters
			thisAlpha <- (mu_alpha_group + s_alpha_group * rnorm(groupSize, 0, 1)) %>%
			  convert_alphaRaw_to_alpha()
			thisBeta <- (mu_beta_group + s_beta_group * rnorm(groupSize, 0, 1)) %>% exp()
			thisSigma <- (mu_soc_group + s_soc_group * rnorm(groupSize, 0, 1)) %>%
			  convert_alphaRaw_to_alpha()
			thisTheta <- mu_theta_group + s_theta_group * rnorm(groupSize, 0, 1)
			# # ----- if individual fit parameters are used directly ------
			# thisAlpha <- fit_SL00_riskID12_parameters$alpha_mean[thisGroup]
			# thisBeta <- fit_SL00_riskID12_parameters$beta_mean[thisGroup]
			# thisSigma <- fit_SL00_riskID12_parameters$soc_mean[thisGroup]
			# thisTheta <- fit_SL00_riskID12_parameters$theta_mean[thisGroup]
		} else { # if condition == individual
			thisAlpha <- (mu_alpha_indiv + s_alpha_indiv * rnorm(groupSize, 0, 1)) %>%
			  convert_alphaRaw_to_alpha()
			thisBeta <- (mu_beta_indiv + s_beta_indiv * rnorm(groupSize, 0, 1)) %>% exp()
			thisTheta <- rep(0, groupSize)
			thisSigma <- rep(0, groupSize)
		}
		## running the task
		for(t in 1:horizon){
			# each individual chooses one option based on his/her choice probability
			# choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
			choices[,t] = mapply(function(p1,p2,p3,p4){ sample(1:numOptions, 1, prob=c(p1,p2,p3,p4), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,], netChoiceProb[3,t,], netChoiceProb[4,t,] )
			# each subject earns some money (if lucky)
			# payoffs[,t] = payoffGenerateBinary(groupSize, choices[,t], rSure, rRisky, payoff_sureL, payoff_sureH, payoff_risky1, payoff_risky2)
			payoffs[,t] = payoffGenerate4Arm_unsync(groupSize, choices[,t], rSure, rSure, rSure, rRisky, payoff_H1, payoff_L1, payoff_H2, payoff_L2, payoff_H3, payoff_L3, payoff_H4, payoff_L4)
			# update choiceCounter and learningRate (if the learning rate is an averaging rule in this simulation.)
			updatingPositions = (choices[,t] + numOptions*(1:groupSize-1))
			
			# value updating
			if(t < horizon) {
				if(t == 1) {
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					# In the first trial, all Q values are updated by the first experience
					QQ[,t+1] = QQ[,t] + thisAlpha * (payoffs[,t] - QQ[,t])
					dim(QQ) = c(numOptions, groupSize, horizon)
					Q = aperm(QQ, c(1,3,2))
				} else {
					# Updating Q value based on Rescorla-Wagner model (Weighted return model)
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					#QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + learningRate[updatingPositions] * (payoffs[,t] - QQ[updatingPositions,t])
					QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + thisAlpha * (payoffs[,t] - QQ[updatingPositions,t])
					dim(QQ) = c(numOptions, groupSize, horizon)
					Q = aperm(QQ, c(1,3,2))
				}

				# update socialFrequency
				## Option's frequency
				for(i in 1:numOptions){
					if(length(which(names(table(choices[,t]))==i))>0) {
						socialFrequency[i,t+1] = socialFrequency[i,t+1] + table(choices[,t])[which(names(table(choices[,t]))==i)][1]
					}
				}
				
				###############
				## Softmax choice base solely on Q values
				###############
				Q_exp = ( Q[,t+1,] * rep(thisBeta, each = numOptions) ) %>% apply(2,expCeiling)
				softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
				freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
				## The followings update the choice probability matrix
				###############
				## Softmax -- END
				###############

				if(model=='Decision-Biasing'){
  				netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-thisSigma)) %>% t() +
  				  apply(freqDepenMatrix, 1, multiplyBeta, beta=thisSigma) %>% t()
				}else{
				  netMatrix = softmaxMatrix
				}
				netChoiceProbAperm = aperm(netChoiceProb, c(1,3,2))
				dim(netChoiceProbAperm) = c(numOptions*groupSize, horizon)
				dim(netMatrix) = c(numOptions*groupSize, 1)
				netChoiceProbAperm[,t+1] = netMatrix
				dim(netChoiceProbAperm) = c(numOptions, groupSize, horizon)
				netChoiceProb = aperm(netChoiceProbAperm, c(1,3,2))
			}
		}

		for(i in 1:groupSize) {
			safeChoiceProb[i,] = 1 - netChoiceProb[4,,i]
		}

		safeChoiceProbMean_b1 = safeChoiceProb[,1:length_of_the_last_block] %>% apply(1, mean) #%>% mean()
		safeChoiceProbMean_b2 = safeChoiceProb[,(length_of_the_last_block+1):horizon] %>% apply(1, mean) #%>% mean()

		choiceProb_option1_b1 = netChoiceProb[1, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option1_b2 = netChoiceProb[1,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()
		choiceProb_option2_b1 = netChoiceProb[2, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option2_b2 = netChoiceProb[2,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()
		choiceProb_option3_b1 = netChoiceProb[3, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option3_b2 = netChoiceProb[3,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()
		choiceProb_option4_b1 = netChoiceProb[4, 1:length_of_the_last_block,] %>% apply(2, mean) #%>% mean()
		choiceProb_option4_b2 = netChoiceProb[4,(length_of_the_last_block+1):horizon,] %>% apply(2, mean) #%>% mean()

		# Submitting this repetition's result
		condition_dummy <- 0
		if(condition=='Group') condition_dummy <- 1
		print(
			data.frame(condition_dummy = rep(condition_dummy, groupSize)
					, hot_stove_susceptibility = thisAlpha * (thisBeta + 1)
					, alpha_mean = thisAlpha
					, beta_mean = thisBeta
					, soc_mean = thisSigma
					, theta_mean = thisTheta
					, proportionSafeChoice_b1 = safeChoiceProbMean_b1
					, proportionSafeChoice_b2 = safeChoiceProbMean_b2
					, proportionRiskyChoice_b1 = choiceProb_option4_b1
					, proportionRiskyChoice_b2 = choiceProb_option4_b2
				)
		)
	}
gc();gc() # rubbish collection
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
social_learning_model_validation_1022_riskID12_data <- social_learning_model_validation_1022_riskID12[[paste("n=", groupSize)]][[paste("condition=", "Individual")]] %>% data.frame()

social_learning_model_validation_1022_riskID12_data <- social_learning_model_validation_1022_riskID12_data %>% rbind(social_learning_model_validation_1022_riskID12[[paste("n=", groupSize)]][[paste("condition=", "Group")]] %>% data.frame())

social_learning_model_validation_1022_riskID12_data$hot_stove_susceptibility_rounded <- (social_learning_model_validation_1022_riskID12_data$hot_stove_susceptibility * 5) %>% round()/5
social_learning_model_validation_1022_riskID12_data$hot_stove_susceptibility_rounded[which(social_learning_model_validation_1022_riskID12_data$hot_stove_susceptibility_rounded>6)] <- 6

social_learning_model_validation_1022_riskID12_data$soc_mean_category <- 'mild'

# overall means
social_learning_model_validation_1022_riskID12_summary <-
	social_learning_model_validation_1022_riskID12_data %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
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
social_learning_model_validation_1022_riskID12_summary_reallyHighSigma <-
	social_learning_model_validation_1022_riskID12_data %>%
	dplyr::filter(soc_mean > 5/10 & soc_mean < 8/10 & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_1022_riskID12_summary_reallyHighSigma$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_1022_riskID12_summary_reallyHighSigma$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha


(ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
	geom_point(data = fit_SL00_riskID12_parameters, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean, colour=soc_mean_SL00_multiVar_LKJ), shape = 20) +
	geom_point(data = fit_AL_indiv_riskID12_parameters, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 18)+ # shape=5: diamond
  geom_line(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid))+
	geom_line(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==1), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)))+
	geom_line(data=social_learning_model_validation_1022_riskID12_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), linetype = 'dashed')+
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Helvetica()+
	xlim(c(0,6.5))+
	labs(
		x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '*(', beta[i], '+1)'))),
		y = 'Mean proportion of choosing\nthe optimal risky option',
		title = 'The 2-risky-2-safe task \n(N = 151)') +
	theme(legend.position = c(0.75, 0.7))+
	theme(legend.title = element_text(size=12))+
	theme(legend.text = element_text(size=11))+
	NULL -> fig6_c)

write.csv(social_learning_model_validation_1022_riskID12_data,
      "social_learning_model_validation_1022_riskID12_data.csv", row.names=FALSE)






# ===========================================================
# The 2-armed bandit
# ==========================================================
fit_AL00_multiVar_LKJ_indiv_0820_globalparameters <- read.csv('fit_AL00_multiVar_LKJ_indiv_0820_globalparameters.csv')
fit_AL00_multiVar_LKJ_indiv_0820_parameters <- read.csv('fit_AL00_multiVar_LKJ_indiv_0820_parameters.csv')
fit_SL00_multiVar_LKJ_0820_globalparameters <- read.csv('fit_SL00_multiVar_LKJ_0820_globalparameters.csv')
fit_SL00_multiVar_LKJ_0820_parameters <- read.csv('fit_SL00_multiVar_LKJ_0820_parameters.csv')

behaviour_main_0820 <- read.csv("behaviour_main_0820.csv")
behaviour_indiv_0820 <-read.csv("behaviour_indiv_0820.csv")
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



##############################################################################
# -- Model simulation as to recover the experimental pattern --
# Is the relationship between the hot-stove susceptibility and performance replicable by the fit model?
##############################################################################
## -- simulated model --
model = 'Decision-Biasing'

## -- settings
# repetition = 100000 #1000
horizon = 70 # = number of trials
numOptions = 2
length_of_the_last_block = 35

# -- a two-armed bandit task with 20/08/2020 setting --
riskPremium = 20/15
rSure = 1
payoff_sureL = 1.5 # 1 - rSure
payoff_sureH = 1.5 # rSure
meanSafePayoff = (rSure*payoff_sureH+(1-rSure)*payoff_sureL)
rRisky = 0.3
payoff_risky1 = 0.5
payoff_risky2 = (riskPremium*meanSafePayoff - payoff_risky1*(1-rRisky))/rRisky
initialExpextation = 0

# -- Simulation
social_learning_model_validation_0820 = list()
conditionList = c('Individual', 'Group')
groupSize = 8

# -- global parameters of the individual condition
mu_alpha_indiv <- fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$variable=='mu_alpha')]
mu_beta_indiv <- fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$variable=='mu_beta')]
s_alpha_indiv <- fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$variable=='s_alpha')]
s_beta_indiv <- fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$mean[which(fit_AL00_multiVar_LKJ_indiv_0820_globalparameters$variable=='s_beta')]

# -- global parameters of the group condition
mu_alpha_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='mu_alpha')]
mu_beta_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='mu_beta')]
mu_soc_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='mu_soc0')]
mu_theta_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='mu_theta')]
s_alpha_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='s_alpha')]
s_beta_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='s_beta')]
s_soc_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='s_soc0')]
s_theta_group <- fit_SL00_multiVar_LKJ_0820_globalparameters$mean[which(fit_SL00_multiVar_LKJ_0820_globalparameters$variable=='s_theta')]


s_time = Sys.time()
for(condition in conditionList) {
	social_learning_model_validation_0820[[paste("n=", groupSize)]][[paste("condition=", condition)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
		## Initial settings
		choices = matrix(nrow=groupSize, ncol=horizon)
		payoffs = matrix(nrow=groupSize, ncol=horizon)
		performance = matrix(nrow=groupSize, ncol=horizon)
		safeChoiceProb = matrix(nrow=groupSize, ncol=horizon)
		isThisBestOption = matrix(nrow=groupSize, ncol=horizon)
		optimalChoiceProb = matrix(nrow=groupSize, ncol=horizon)
		expectedPerformance = matrix(nrow=groupSize, ncol=horizon)
		Q = array(dim = c(numOptions, horizon, groupSize))
		choiceCounter = array(1, dim = c(numOptions, groupSize))
		netChoiceProb = array(dim = c(numOptions, horizon, groupSize))
		netChoiceProb[,1,] = 1/numOptions
		Q[,1,] = initialExpextation
		socialFrequency = matrix(nrow=numOptions, ncol=horizon)
		socialFrequency[,] = 1e-1
		## Setting individual parameters
		if(condition == 'Group') {
			thisGroup <- sample(1:nrow(fit_parameters_group_SL00_mcmc), size=groupSize, replace=TRUE)
			thisAlpha <- (mu_alpha_group + s_alpha_group * rnorm(groupSize, 0, 1)) %>% convert_alphaRaw_to_alpha()
			thisBeta <- ( mu_beta_group + s_beta_group * rnorm(groupSize, 0, 1) ) %>% exp()
			thisSigma <- (mu_soc_group + s_soc_group * rnorm(groupSize, 0, 1)) %>% convert_alphaRaw_to_alpha()
			thisTheta <- mu_theta_group + s_theta_group * rnorm(groupSize, 0, 1)
			
		} else {
			thisAlpha <- (mu_alpha_indiv + s_alpha_indiv * rnorm(groupSize, 0, 1)) %>% convert_alphaRaw_to_alpha()
			thisBeta <- ( mu_beta_indiv + s_beta_indiv * rnorm(groupSize, 0, 1) ) %>% exp()
			thisSigma <- rep(0, groupSize)
			thisTheta <- rep(0, groupSize)
		}
		## running the task
		for(t in 1:horizon){
			# each individual chooses one option based on his/her choice probability
			choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
			# each subject earns some money (if lucky)
			#payoffs[,t] = payoffGenerateGaussian(groupSize, choices[,t], mu_sure, mu_risky, sd_sure, sd_risky)
			payoffs[,t] = payoffGenerateBinary(groupSize, choices[,t], rSure, rRisky, payoff_sureL, payoff_sureH, payoff_risky1, payoff_risky2)
			# update choiceCounter and learningRate (if the learning rate is an averaging rule in this simulation.)
			updatingPositions = (choices[,t] + numOptions*(1:groupSize-1))
			
			if(t < horizon) {
				if(t == 1) {
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					# In the first trial, all Q values are updated by the first experience
					QQ[,t+1] = QQ[,t] + thisAlpha * (payoffs[,t] - QQ[,t])
					dim(QQ) = c(numOptions, groupSize, horizon)
					Q = aperm(QQ, c(1,3,2))
				} else {
					# Updating Q value based on Rescorla-Wagner model (Weighted return model)
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					#QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + learningRate[updatingPositions] * (payoffs[,t] - QQ[updatingPositions,t])
					QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + thisAlpha * (payoffs[,t] - QQ[updatingPositions,t])
					dim(QQ) = c(numOptions, groupSize, horizon)
					Q = aperm(QQ, c(1,3,2))
				}
				# update socialFrequency
				## Option 1's frequency
				if(length(which(names(table(choices[,t]))==1))>0) {
					socialFrequency[1,t+1] = socialFrequency[1,t+1] + table(choices[,t])[which(names(table(choices[,t]))==1)][1]
				}
				## Option 2's frequency
				if(length(which(names(table(choices[,t]))==2))>0) {
					socialFrequency[2,t+1] = socialFrequency[2,t+1] + table(choices[,t])[which(names(table(choices[,t]))==2)][1]
				}

				###############
				## Softmax choice base solely on Q values
				###############
				Q_exp = ( Q[,t+1,] * rep(thisBeta, each = numOptions) ) %>% apply(2,expCeiling)
				softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
				freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
				## The followings update the choice probability matrix
				###############
				## Softmax -- END
				###############

				if(model=='Decision-Biasing'){
				netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-thisSigma)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=thisSigma) %>% t()
				}else{
				netMatrix = softmaxMatrix
				}
				netChoiceProbAperm = aperm(netChoiceProb, c(1,3,2))
				dim(netChoiceProbAperm) = c(numOptions*groupSize, horizon)
				dim(netMatrix) = c(numOptions*groupSize, 1)
				netChoiceProbAperm[,t+1] = netMatrix
				dim(netChoiceProbAperm) = c(numOptions, groupSize, horizon)
				netChoiceProb = aperm(netChoiceProbAperm, c(1,3,2))
			}
		}

		for(i in 1:groupSize) {
			safeChoiceProb[i,] = netChoiceProb[1,,i]
		}

		safeChoiceProbMean_b1 = safeChoiceProb[,1:length_of_the_last_block] %>% apply(1, mean) #%>% mean()
		safeChoiceProbMean_b2 = safeChoiceProb[,(length_of_the_last_block+1):horizon] %>% apply(1, mean) #%>% mean()

		risky_choice_count_b1 <- (choices[,1:length_of_the_last_block] - 1) %>% apply(1, sum)
		risky_choice_count_b2 <- (choices[,(length_of_the_last_block+1):horizon] - 1) %>% apply(1, sum)

		# Submitting this repetition's result
		condition_dummy <- 0
		if(condition=='Group') condition_dummy <- 1
		print(
			data.frame(
				condition_dummy = rep(condition_dummy, groupSize),
				hot_stove_susceptibility = thisAlpha * (thisBeta + 1),
				alpha_mean = thisAlpha,
				beta_mean = thisBeta,
				soc_mean = thisSigma,
				theta_mean = thisTheta,
				proportionSafeChoice_b1 = safeChoiceProbMean_b1,
				proportionSafeChoice_b2 = safeChoiceProbMean_b2,
				proportionRiskyChoice_b1 = 1 - safeChoiceProbMean_b1,
				proportionRiskyChoice_b2 = 1 - safeChoiceProbMean_b2,
				risky_choice_count_b1 = risky_choice_count_b1,
				risky_choice_count_b2 = risky_choice_count_b2
				)
		)
	}
gc();gc() # rubbish collection
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
social_learning_model_validation_0820_data <- social_learning_model_validation_0820[[paste("n=", groupSize)]][[paste("condition=", "Individual")]] %>% data.frame()

social_learning_model_validation_0820_data <- social_learning_model_validation_0820_data %>% rbind(social_learning_model_validation_0820[[paste("n=", groupSize)]][[paste("condition=", "Group")]] %>% data.frame())

social_learning_model_validation_0820_data$hot_stove_susceptibility_rounded <- (social_learning_model_validation_0820_data$hot_stove_susceptibility * 5) %>% round()/5
social_learning_model_validation_0820_data$hot_stove_susceptibility_rounded[which(social_learning_model_validation_0820_data$hot_stove_susceptibility_rounded>6)] <- 6

social_learning_model_validation_0820_data$soc_mean_category <- 'mild'

# overall means
social_learning_model_validation_0820_summary <-
	social_learning_model_validation_0820_data %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% mean(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_0820_summary$proportionRiskyChoice_b2_lower <-
	(social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_mean - social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_0820_summary$n)) %>% convert_alphaRaw_to_alpha
social_learning_model_validation_0820_summary$proportionRiskyChoice_b2_upper <-
	(social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_mean + social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_sd / sqrt(social_learning_model_validation_0820_summary$n)) %>% convert_alphaRaw_to_alpha

social_learning_model_validation_0820_summary$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_0820_summary$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha


# modest social learners' means
social_learning_model_validation_0820_summary_reallyHighSigma <-
	social_learning_model_validation_0820_data %>%
	dplyr::filter(soc_mean > 5/10 & soc_mean < 8/10 & hot_stove_susceptibility_rounded < 6) %>%
	group_by(condition_dummy, hot_stove_susceptibility_rounded, soc_mean_category) %>%
	summarise(
		proportionRiskyChoice_b2_mean = mean(proportionRiskyChoice_b2),
		proportionRiskyChoice_b2_sd = sd(proportionRiskyChoice_b2),
		raw_proportionRiskyChoice_b2_mean = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% median(),
		raw_proportionRiskyChoice_b2_sd = proportionRiskyChoice_b2 %>% convert_alpha_to_alphaRaw() %>% sd(),
		soc_mean = mean(soc_mean),
		n = n()
		)

social_learning_model_validation_0820_summary_reallyHighSigma$proportionRiskyChoice_b2_mid <-
	social_learning_model_validation_0820_summary_reallyHighSigma$raw_proportionRiskyChoice_b2_mean %>% convert_alphaRaw_to_alpha

fit_parameters_group_SL00_mcmc$soc_mean <- fit_parameters_group_SL00_mcmc$soc_mean_SL00_multiVar_LKJ

# Plot
(ggplot() +
	#geom_ribbon(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='grey20', alpha=1/2)+
	#geom_ribbon(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==1&soc_mean_category=='mild'), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='orange', alpha=1/2)+
	geom_segment(aes(x=0,xend=6,y=0.5,yend=0.5),colour="grey30", size=0.5) +
	geom_point(data = fit_parameters_group_SL00_mcmc, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean, colour=soc_mean), shape = 20) +
	geom_point(data = parameterfit_indiv_AL00_0820, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 18)+ # shape=5: diamond
	geom_line(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid))+
	geom_line(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==1), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)))+
	geom_line(data=social_learning_model_validation_0820_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), linetype = 'dashed')+
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Helvetica()+
	xlim(c(0,6.5))+
	labs(
		x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '*(', beta[i], '+1)'))),
		y = 'Mean proportion of choosing\nthe optimal risky option',
		title = 'The 1-risky-1-safe task \n(N = 168)') +
	#theme(legend.position = c(0.85, 0.5))+
	theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12))+
	theme(legend.text = element_text(size=11))+
	NULL -> fig6_a
)

write.csv(social_learning_model_validation_0820_data,
      "social_learning_model_validation_0820_data.csv", row.names=FALSE)



## Figure -- experimental result
figure_exp_model_pred <- ggarrange(fig6_a, fig6_b, fig6_c
	# ,  common.legend = TRUE
	# ,  legend = 'right'
	, labels = c('','',''), ncol = 3, align = 'v'
	)

# figure_exp_model_pred <- plot_grid(fig6_a, fig6_b, fig6_c
# 	, labels = c('','',''), ncol = 3, align = 'v')

ggsave(file = "figure_exp_model_pred.pdf", plot = figure_exp_model_pred, dpi = 600, width = 12, height = 4.5
	, device = cairo_pdf
	# , device = 'pdf'
	)

