###############################################################################################
##
## 3-safe 1-risky bandit with exp condition 0 setting
## Parameter Recovery Test
##
## Wataru Toyokawa
##
###############################################################################################
if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace

# Loading
library(tidyverse)
library(cowplot)
library(foreach)
library(MASS)
library(doParallel)
library(lme4)
library(mvtnorm)
#registerDoParallel(detectCores()) # this uses as many core as available
registerDoParallel(detectCores())

## Load Functions
source('~/Dropbox/wataru/papers/RiskySocialLearning/four_armed_simulation/functions.R')

# path -- drop box
dropbox_path <- "~/Dropbox/wataru/papers/RiskySocialLearning/experiment/parameter_recovery/SL00/recovery_param_expfit/"

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


##############################
# -- Global parameters (model SL11) --
##############################
# parameter set exp
mu_alpha <- -1.8
mu_beta <- 1.6
mu_soc0 <- -2.0
mu_theta <- 1.4
s_alpha <- 1.6
s_beta <- 1.1
s_soc0 <- 1.7
s_theta <- 1.5

# # parameter set 1
# mu_alpha <- -1.8
# mu_beta <- 1.5
# # mu_epsilon <- 1.5
# mu_soc0 <- -2
# # mu_soc_slope <- -1
# mu_theta <- 2
# s_alpha <- 1
# s_beta <- 1
# # s_epsilon <- 1
# s_soc0 <- 1
# # s_soc_slope <- 1
# s_theta <- 1

# # parameter set 1
# mu_alpha <- -1.2
# mu_beta <- 1.8
# # mu_epsilon <- 1.5
# mu_soc0 <- -1.5
# # mu_soc_slope <- -1
# mu_theta <- 1.5
# s_alpha <- 1.2
# s_beta <- 1.2
# # s_epsilon <- 1
# s_soc0 <- 1.2
# # s_soc_slope <- 1
# s_theta <- 1.2



# Correlations between individual parameters
# That is, for example, individuals who have a large log_beta tend to have a negative epsilon
# whereas those who have small log_beta might have a positive epsilon (i.e. a negative correlation)
sigma_beta_theta <- matrix(c(s_beta, -0.4, -0.4, s_theta), ncol=2)

# sigma_beta_epsilon <- matrix(c(s_beta, -0.3, -0.3, s_epsilon), ncol=2)
# sigma_soc0_soc_slope <- matrix(c(s_soc0, -0.3, -0.3, s_soc_slope), ncol=2)

##############################
# -- Simulation setting --
##############################
groupSize <- 6
repetition <- 18

# -- Simulation
pseudoData_4ab_riskID11 = list()
conditionList = c('Group')


s_time = Sys.time()
for(condition in conditionList) {
	pseudoData_4ab_riskID11[[paste("n=", groupSize)]][[paste("condition=", condition)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
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
			beta_theta_mvnorm <- rmvnorm(n=groupSize , mean=c(mu_beta, mu_theta) , sigma=sigma_beta_theta)
			# beta_mvnorm <- rmvnorm(n=groupSize , mean=c(mu_beta, mu_epsilon) , sigma=sigma_beta_epsilon)
			# sigma_mvnorm <- rmvnorm(n=groupSize , mean=c(mu_soc0, mu_soc_slope) , sigma=sigma_soc0_soc_slope)
			this_alpha <- (mu_alpha + s_alpha * rnorm(groupSize, 0, 1)) %>% convert_alphaRaw_to_alpha()
			this_log_beta <- beta_theta_mvnorm[,1]
			# this_epsilon <- beta_mvnorm[,2]
			# this_soc0 <- beta_soc0_mvnorm[,2]
			# this_soc_slope <- sigma_mvnorm[,2]
			# this_log_beta <- mu_beta + s_beta * rnorm(groupSize, 0, 1)
			# this_epsilon <- mu_epsilon + s_epsilon * rnorm(groupSize, 0, 1)
			this_soc0 <- mu_soc0 + s_soc0 * rnorm(groupSize, 0, 1)
			# this_soc_slope <- mu_soc_slope + s_soc_slope * rnorm(groupSize, 0, 1)
			# this_theta <- mu_theta + s_theta * rnorm(groupSize, 0, 1)
			this_theta <- beta_theta_mvnorm[,2]
			this_sigma <- this_soc0 %>% convert_alphaRaw_to_alpha()
		} else {
			# Individual condition
		}
		## running the task
		for(t in 1:horizon){
			# == choice and payoff (4-armed bandit)
			choices[,t] = mapply(function(p1,p2,p3,p4){ sample(1:numOptions, 1, prob=c(p1,p2,p3,p4), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,], netChoiceProb[3,t,], netChoiceProb[4,t,] )
			payoffs[,t] = payoffGenerate4Arm_unsync(groupSize, choices[,t], rSure, rSure, rSure, rRisky, payoff_H1, payoff_L1, payoff_H2, payoff_L2, payoff_H3, payoff_L3, payoff_H4, payoff_L4)

			# # == choice and payoff (2-armed bandit)
			# # each individual chooses one option based on his/her choice probability
			# choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
			# payoffs[,t] = payoffGenerateBinary(groupSize, choices[,t], rSure, rRisky, payoff_sureL, payoff_sureH, payoff_risky1, payoff_risky2)

			# update choiceCounter and learningRate (if the learning rate is an averaging rule in this simulation.)
			updatingPositions = (choices[,t] + numOptions*(1:groupSize-1))
			# -- (if the learning rate is an averaging rule in this simulation.) --
			#choiceCounter = aperm(choiceCounter, c(2,1)) # transform the choiceCounter matrix
			#dim(choiceCounter) = c(numOptions*groupSize) # reduce dimension
			#choiceCounter[updatingPositions] = choiceCounter[updatingPositions] + 1
			#learningRate = 1/choiceCounter # Learning rate is now set
			#dim(choiceCounter) = c(numOptions, groupSize)
			# -- END --
			if(t < horizon) {
				if(t == 1) {
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					# In the first trial, all Q values are updated by the first experience
					QQ[,t+1] = QQ[,t] + this_alpha * (payoffs[,t] - QQ[,t])
					dim(QQ) = c(numOptions, groupSize, horizon)
					Q = aperm(QQ, c(1,3,2))
				} else {
					# Updating Q value based on Rescorla-Wagner model (Weighted return model)
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(numOptions*groupSize, horizon)
					#QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + learningRate[updatingPositions] * (payoffs[,t] - QQ[updatingPositions,t])
					QQ[updatingPositions,t+1] = QQ[updatingPositions,t] + this_alpha * (payoffs[,t] - QQ[updatingPositions,t])
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
				Q_exp = ( Q[,t+1,] * rep( exp( this_log_beta ) , each = numOptions) ) %>% apply(2,expCeiling)
				# Q_exp = ( Q[,t+1,] * rep( exp( this_log_beta+this_epsilon*(t+1)/horizon ) , each = numOptions) ) %>% apply(2,expCeiling)
				softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
				freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], this_theta, numOptions)
				## The followings update the choice probability matrix
				###############
				## Softmax -- END
				###############

				# this_sigma = ( this_soc0 + this_soc_slope * (t+1)/horizon ) %>% convert_alphaRaw_to_alpha()
				#soc = 0 # soc = 0 indicates asocial learning (i.e. no social info use)
				##netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-epsilon*numOptions)) %>% t() + epsilon
				#netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
				if(model=='Decision-Biasing'){
					netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-this_sigma)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=this_sigma) %>% t()
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


		print(
			data.frame(
				choices = as.vector(t(choices))
				, payoffs = as.vector(t(payoffs))
				, indiv = rep(1:groupSize, each=horizon) + (rep - 1) * groupSize
				, group = rep
				, trial = rep(1:horizon, groupSize)
				, alpha = rep(this_alpha, each = horizon)
				, log_beta = rep(this_log_beta, each = horizon)
				# , epsilon = rep(this_epsilon, each = horizon)
				, soc0 = rep(this_soc0, each = horizon)
				, sigma = rep(this_sigma, each = horizon)
				# , soc_slope = rep(this_soc_slope, each = horizon)
				, theta = rep(this_theta, each = horizon)
				, socialFreq_1 = socialFrequency[1,]
				, socialFreq_2 = socialFrequency[2,]
				, socialFreq_3 = socialFrequency[3,]
				, socialFreq_4 = socialFrequency[4,]
				)
		)

	}
gc();gc() # rubbish collection
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
pseudoData_4ab_riskID11_data <- pseudoData_4ab_riskID11[[paste("n=", groupSize)]][[paste("condition=", "Group")]] %>% data.frame()

# pseudoData_4ab_riskID11_data$net_beta <- ( pseudoData_4ab_riskID11_data$log_beta + pseudoData_4ab_riskID11_data$epsilon * (pseudoData_4ab_riskID11_data$trial/horizon) ) %>% exp()

# pseudoData_4ab_riskID11_data$net_sigma <- ( pseudoData_4ab_riskID11_data$soc0 + pseudoData_4ab_riskID11_data$soc_slope *(pseudoData_4ab_riskID11_data$trial/horizon) ) %>% convert_alphaRaw_to_alpha


true_global_parameters_4ab_riskID11 <- data.frame(
	mu_alpha=mu_alpha
	, mu_beta=mu_beta
	# , mu_epsilon=mu_epsilon
	, mu_soc0=mu_soc0
	# , mu_soc_slope=mu_soc_slope
	, mu_theta=mu_theta
	, s_alpha=s_alpha
	, s_beta=s_beta
	# , s_epsilon=s_epsilon
	, s_soc0=s_soc0
	# , s_soc_slope=s_soc_slope
	, s_theta=s_theta
	)

true_individual_parameters_4ab_riskID11 <- pseudoData_4ab_riskID11_data %>%
	group_by(indiv) %>%
	summarise(alpha = mean(alpha)
		, log_beta = mean(log_beta)
		# , epsilon = mean(epsilon)
		, soc0 = mean(soc0)
		, sigma = mean(sigma)
		# , soc_slope = mean(soc_slope)
		, theta = mean(theta)
		, riskyChoice = sum(choices - 1)/horizon
		)

# --- save files ---
write.csv(pseudoData_4ab_riskID11_data, paste0(dropbox_path, "pseudoData_4ab_riskID11_data.csv"), row.names=FALSE)
write.csv(true_global_parameters_4ab_riskID11, paste0(dropbox_path, "true_global_parameters_4ab_riskID11.csv"), row.names=FALSE)
write.csv(true_individual_parameters_4ab_riskID11, paste0(dropbox_path, "true_individual_parameters_4ab_riskID11.csv"), row.names=FALSE)





