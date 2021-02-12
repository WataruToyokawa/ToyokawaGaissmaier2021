###############################################################################################
##
## The effect of individual differences in the hot stove susceptibility on
## the social learning in the basic Denrell (2007) task.
##
## 23 November 2020
###############################################################################################
# Aim: Thus far I have shown that social influence can rescue social learners from the hot stove effect
# However, it is still unclear that whether such individuals who have a low susceptibility to the hot stove
# can also benefit from social learning, or the collective rescue emerges only
# in the expense of the detrimental effect for them


if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace

# Loading
library(tidyverse)
library(cowplot)
library(foreach)
library(MASS)
library(doParallel)
registerDoParallel(detectCores())

## Load Functions
source('~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/functions.R')

## -- simulated model --
model = 'Decision-Biasing'

#### simulation 2
# Systematically varying the individual heterogeneity in each learning parameters

# Individual variations in the susceptibility to the hot stove effect -- alpha
# where group members are all susceptible to the hot stove effect
groupSizeList = c(5)
groupSize = groupSizeList[1]
repetition = 20000
conditions = c("Individual", "Group")
## -- Global parameter setting --
sigmaList = 0.3
# sigmaList = seq(0.25, 0.35, 0.1)
# sigmaList = seq(0.05, 0.55, 0.1)
# thetaList = c(1, 2, 4)
thetaList = 2
alphaList = list(
	'1'=rep(0.5, groupSize) ,
	#'2'=c(0.45, 0.50, 0.55, 0.60, 0.65),
	'3'=c(0.3, 0.4, 0.5, 0.6, 0.7),
	#'4'=c(0.25, 0.32, 0.55, 0.78, 0.85),
	#'5'=c(0.2, 0.32, 0.55, 0.78, 0.9),
	# '6'=c(0.15, 0.24, 0.41, 0.8, 0.9),
	'7'=c(0.15, 0.21, 0.24, 0.93, 0.97), # where majority is risk-seeking
	'8'=c(0.05, 0.21, 0.49, 0.85, 0.9)
	)
# alpha_list_indiv = c(0.15, 0.21, 0.24, 0.93, 0.97)
# alpha_list_indiv = c(0.05, 0.21, 0.49, 0.85, 0.9)
alpha_list_indiv = c(0.05, 0.15, 0.21, 0.24, 0.3, 0.4, 0.5, 0.7, 0.85, 0.9, 0.97)
invTemperatureList = c(7)
#hot_stove_susceptibility = (invTemperatureList+1) * alphaList
variation_level_list = 1:length(alphaList)

## -- simulation parameters
horizon = 150 # = number of trials
numOptions = 2

## -- transformed parameters
#alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList)
sigmaRawList = mapply(convert_alpha_to_alphaRaw, sigmaList)

## -- Individual variation
variationAlphaRaw = 0.01 #1.5
variationBeta = 0.01#0.7
variatonAnnealing = 0.01#1.5
variationSigmaRaw = 0.01
variationTheta = 0.01

# -- Setting a two-armed bandit task --
riskPremium = 1.5
mu_sure = 1
mu_risky = mu_sure * riskPremium
sd_sure = 0.05
sd_risky = 1
initialExpextation = 0

# -- Initial Settings --
heterogeneous_groups_denrell2007task_alpha = list()


sigmaFlag = 0
# -- simulation --
s_time = Sys.time()
for (indivOrGroup in conditions) {
	if (sigmaFlag==1 & indivOrGroup=="Individual") next
	for (variation in variation_level_list) {

		if(indivOrGroup=='Individual') {
			variation <- 0
			groupSize <- alpha_list_indiv %>% length()
		}
		else {
			groupSize <- groupSizeList[1]
		}

		heterogeneous_groups_denrell2007task_alpha[[paste("indivOrGroup=",indivOrGroup)]][[paste("variation_level=",variation)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
			set.seed(rep)
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
			if(indivOrGroup=="Individual"){
				alphaRawList <- alpha_list_indiv %>%
					mapply(convert_alpha_to_alphaRaw, .) %>%
					sample() # shuffling individual order to avoid the weird behaviour
			} else {
				alphaRawList <- alphaList[variation][[1]] %>%
					mapply(convert_alpha_to_alphaRaw, .) %>%
					sample()
				# alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList[variation][[1]])
			}
			thisAlpha = ( alphaRawList + variationAlphaRaw * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
			thisBeta = invTemperatureList + variationBeta  * rt(groupSize, df = 14, ncp = 0)
			thisBeta[which(thisBeta<0)] <- 0
			thisHotStoveSusceptibility <- (alphaRawList %>% mapply(convert_alphaRaw_to_alpha,.) * (invTemperatureList + 1))
			# thisHotStoveSusceptibility <- thisAlpha * (thisBeta+1)
			#thisAnnealing = annealing + variatonAnnealing * rt(groupSize, df = 14, ncp = 0)
			if(indivOrGroup == 'Individual'){
				thisSigma <- rep(0, groupSize)
				thisTheta <- rep(0, groupSize)
				theta <- 0
			} else {
				thisSigma = ( runif(groupSize, min=min(sigmaRawList), max=max(sigmaRawList)) + variationSigmaRaw * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
				thisTheta = thetaList + variationTheta * rt(groupSize, df = 14, ncp = 0)
			}
			thisAnnealing = 0

			for(t in 1:horizon){
				# each individual chooses one option based on his/her choice probability
				choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
				# each subject earns some money (if lucky)
				payoffs[,t] = payoffGenerateGaussian(groupSize, choices[,t], mu_sure, mu_risky, sd_sure, sd_risky)
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
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# In the original article, March only considered a proportional choice
					# If you want to implement softmax rule, you should modify this
					#proportionalChoiceMatrix = Q[,t+1,] %>% apply(1, divideVector, denominator = apply(Q[,t+1,],2,sum)) %>% t()

					###############
					## Softmax choice base solely on Q values
					###############
					#Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta = (thisBeta + thisAnnealing * (t+1)/horizon) ) %>% t() %>% apply(2,expCeiling)
					Q_exp = ( Q[,t+1,] * rep((thisBeta + thisAnnealing * (t+1)/horizon), each = numOptions) ) %>% apply(2,expCeiling)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
					## The followings update the choice probability matrix
					###############
					## Softmax -- END
					###############

					#soc = 1/(1+exp(-(soc_raw)))
					#soc = 0 # soc = 0 indicates asocial learning (i.e. no social info use)
					##netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-epsilon*numOptions)) %>% t() + epsilon
					#netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
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

			safeChoiceProbMean = safeChoiceProb[,round(horizon/2+1):horizon] %>% apply(MARGIN=1, mean)
			# safeChoiceProbMean = safeChoiceProb[,round(2*horizon/3+1):horizon] %>% apply(MARGIN=1, mean)

			# Submitting this repetition's result
			print(
				data.frame(
					indivOrGroup = rep(indivOrGroup, groupSize),
					variation_level = variation,
					sub = 1:groupSize,
					alpha = thisAlpha,
					beta = thisBeta,
					hot_stove_susceptibility = thisHotStoveSusceptibility,
					sigma = thisSigma,
					theta = thisTheta,
					#mu_theta = rep(theta, groupSize),
					safeChoiceProbMean = safeChoiceProbMean,
					rep = rep
					)
			)
		}
		if(indivOrGroup=="Individual") sigmaFlag <- 1
		gc();gc() # rubbish collection
	}
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
heterogeneous_groups_denrell2007task_alpha_data = heterogeneous_groups_denrell2007task_alpha[[paste("indivOrGroup=",conditions[1])]][[paste("variation_level=",0)]] %>% data.frame()

for(variation in variation_level_list) {
	#for(sigma in sigmaList) {
		#for (groupSize in groupSizeList) {
			#for (alpha in alphaList) {
				#for (invTemperature in invTemperatureList) {
					#if(variation!=variation_level_list[1]) {
						heterogeneous_groups_denrell2007task_alpha_data = heterogeneous_groups_denrell2007task_alpha_data %>% rbind(data.frame(heterogeneous_groups_denrell2007task_alpha[[paste("indivOrGroup=",'Group')]][[paste("variation_level=",variation)]]))
					#}
				#}
			#}
		#}
	#}
}

#names(heterogeneous_groups_denrell2007task_alpha_data) = c('groupSize', 'learningRate', 'invTemperature', 'copyRate', 'conformityExp', 'proportionSafeChoice')

heterogeneous_groups_denrell2007task_alpha_data$safeChoiceProbMean_raw <- heterogeneous_groups_denrell2007task_alpha_data$safeChoiceProbMean %>% convert_alpha_to_alphaRaw()

heterogeneous_groups_denrell2007task_alpha_summary <- heterogeneous_groups_denrell2007task_alpha_data %>%
	group_by(indivOrGroup, variation_level, hot_stove_susceptibility) %>%
	summarise(
		mean_proportionSafeChoice_direct = mean(safeChoiceProbMean, na.rm = TRUE),
	    #sd_proportionSafeChoice = sd(safeChoiceProbMean, na.rm = TRUE),
	    mean_proportionSafeChoice = mean(safeChoiceProbMean_raw, na.rm = TRUE) %>% convert_alphaRaw_to_alpha(),
	    median_proportionSafeChoice_raw = median(safeChoiceProbMean_raw, na.rm = TRUE),
	    sd_proportionSafeChoice_raw = sd(safeChoiceProbMean_raw, na.rm = TRUE)
		)

heterogeneous_groups_denrell2007task_alpha_summary$hot_stove_susceptibility_slided <- heterogeneous_groups_denrell2007task_alpha_summary$hot_stove_susceptibility + rnorm(nrow(heterogeneous_groups_denrell2007task_alpha_summary), 0, 0.01)

heterogeneous_groups_denrell2007task_alpha_summary_global <-
	heterogeneous_groups_denrell2007task_alpha_summary %>%
	group_by(indivOrGroup, variation_level) %>%
	summarise(
		mean_hot_stove_susceptibility = mean(hot_stove_susceptibility, na.rm=TRUE),
		mean_safeChoiceProb = mean(mean_proportionSafeChoice_direct, na.rm=TRUE)
		# mean_safeChoiceProb = mean(median_proportionSafeChoice_raw, na.rm=TRUE) %>% convert_alphaRaw_to_alpha()
	)

Pr_when_beta = function (X, beta) {
	Z = -beta/(2*(X/(beta+1))-2) * (X - 2)
	return_vector <- 1 / (1 + exp(Z))
	return_vector[which(X >= beta+1)] <- NA
	return(return_vector)
}

(heterogeneous_groups_denrell2007task_alpha_summary %>%
	dplyr::filter(indivOrGroup=='Group') %>%
	ggplot()+
	geom_vline(xintercept=2, linetype='dashed', colour='grey60')+
	geom_hline(yintercept=0.5, linetype='dashed', colour='grey60')+
	stat_function(data=data.frame(X=c(2,9), invTemperature_factor='β = 7', indivOrGroup='Individual'), fun=Pr_when_beta, args=list(beta=7)) +
	# group condition
	geom_point(aes(hot_stove_susceptibility, 1- mean_proportionSafeChoice, colour=as.factor(variation_level)),alpha=3/3)+
	geom_line(aes(hot_stove_susceptibility, 1-mean_proportionSafeChoice, group=variation_level, colour=as.factor(variation_level)), alpha=3/3)+
	geom_point(data=heterogeneous_groups_denrell2007task_alpha_summary_global%>%filter(variation_level!=0), aes(mean_hot_stove_susceptibility, 1-mean_safeChoiceProb, colour=as.factor(variation_level)), size=4, shape=18, alpha = 2/3)+
	# individual condition
	geom_point(data = heterogeneous_groups_denrell2007task_alpha_summary%>%filter(variation_level==0), aes(hot_stove_susceptibility, 1- mean_proportionSafeChoice), colour='grey30',alpha=3/3)+
	geom_line(data = heterogeneous_groups_denrell2007task_alpha_summary%>%filter(variation_level==0), aes(hot_stove_susceptibility, 1-mean_proportionSafeChoice, group=variation_level), colour='grey30', alpha=3/3)+
	# geom_point(data=heterogeneous_groups_denrell2007task_alpha_summary_global%>%filter(variation_level==0), aes(mean_hot_stove_susceptibility, 1-mean_safeChoiceProb), colour='grey30', size=4, shape=18, alpha = 2/3)+
	labs(x=expression(paste("Hot stove susceptibility ",alpha[i],"", (beta[i]+1) ,sep="")),
		y="Proportion of choosing\nthe risky option",
		colour = "Individual\nheterogeneity",
		title=expression(paste("Heterogeneous ", alpha[i], sep=""))
		)+
	ylim(c(0,1))+
	xlim(c(0,8))+
	myTheme_Helvetica()+
	scale_colour_viridis_d(direction=-1)+
	theme(legend.position = 'NaN')+
	theme(plot.title = element_text(vjust = - 10, hjust = 0.7))+
	NULL -> heterogeneous_groups_denrell2007task_alphaEffect_plot)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_alphaEffect_plot.png", plot = heterogeneous_groups_denrell2007task_alphaEffect_plot, dpi = 300, width = 4.5, height = 3)

write.csv(heterogeneous_groups_denrell2007task_alpha_summary,
			"~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_alpha_summary.csv",
			row.names=FALSE)



# Individual variations in the susceptibility to the hot stove effect -- beta
# where group members are all susceptible to the hot stove effect
repetition = 20000
conditions = c("Individual", "Group")
## -- Global parameter setting --
sigmaList = 0.3
# sigmaList = seq(0.25, 0.35, 0.1)
# sigmaList = seq(0.05, 0.55, 0.1)
# thetaList = c(1, 2, 4)
thetaList = 2
# alphaList = list(
# 	'1'=rep(0.55, groupSize),
# 	'2'=c(0.45, 0.50, 0.55, 0.60, 0.65),
# 	'3'=c(0.35, 0.45, 0.55, 0.65, 0.75),
# 	'4'=c(0.2, 0.32, 0.55, 0.78, 0.9)
# 	)
alphaList = 0.50
# invTemperatureList = c(7)
invTemperatureList = list(
	'1'=rep(7, groupSize),
	'2'=c(5.75, 6.375, 7, 7.625, 8.25),
	'3' = c(2, 2.8, 3.6, 11.6, 15),
	'4' = c(1.5, 2, 4.9, 11.6, 15)
	#'3'=c(4.5, 5.75, 7, 8.25, 9.5),
	#'4'=c(2, 4.5, 7, 9.5, 12)
	)
invTemperatureList_indiv = c(1, 1.5, 2, 2.8, 3.6, 4.9, 7, 8.25, 11.6, 15)
#hot_stove_susceptibility = (invTemperatureList+1) * alphaList
variation_level_list = 1:length(invTemperatureList)

## -- simulation parameters
groupSizeList = c(5)
horizon = 150 # = number of trials
numOptions = 2

## -- transformed parameters
#alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList)
sigmaRawList = mapply(convert_alpha_to_alphaRaw, sigmaList)

## -- Individual variation
variationAlphaRaw = 0.01 #1.5
variationBeta = 0.01#0.7
variatonAnnealing = 0.01#1.5
variationSigmaRaw = 0.01
variationTheta = 0.01

# -- Setting a two-armed bandit task --
riskPremium = 1.5
mu_sure = 1
mu_risky = mu_sure * riskPremium
sd_sure = 0.05
sd_risky = 1
initialExpextation = 0

# -- Initial Settings --
heterogeneous_groups_denrell2007task_beta = list()
groupSize = groupSizeList[1]

sigmaFlag = 0
# -- simulation --
s_time = Sys.time()
for (indivOrGroup in conditions) {
	if (sigmaFlag==1 & indivOrGroup=="Individual") next
	for (variation in variation_level_list) {

		if(indivOrGroup=='Individual') {
			variation <- 0
			groupSize <- invTemperatureList_indiv %>% length()
		}
		else {
			groupSize <- groupSizeList[1]
		}

		heterogeneous_groups_denrell2007task_beta[[paste("indivOrGroup=",indivOrGroup)]][[paste("variation_level=",variation)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
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
			alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList[1][[1]])
			thisAlpha = ( alphaRawList + variationAlphaRaw * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
			if(indivOrGroup=="Individual"){
				this_invTemperatureList_indiv <- invTemperatureList_indiv %>% sample()
				thisBeta <- (this_invTemperatureList_indiv + variationBeta  * rt(groupSize, df = 14, ncp = 0))
				thisHotStoveSusceptibility <- (alphaRawList %>% mapply(convert_alphaRaw_to_alpha,.)) * (this_invTemperatureList_indiv + 1)
			} else {
				this_invTemperatureList <- invTemperatureList[variation][[1]] %>% sample()
				thisBeta <- this_invTemperatureList + variationBeta  * rt(groupSize, df = 14, ncp = 0)
				thisHotStoveSusceptibility <- (alphaRawList %>% mapply(convert_alphaRaw_to_alpha,.)) * (this_invTemperatureList + 1)
			}
			thisBeta[which(thisBeta<0)] <- 0
			#thisAnnealing = annealing + variatonAnnealing * rt(groupSize, df = 14, ncp = 0)
			if(indivOrGroup == 'Individual'){
				thisSigma <- rep(0, groupSize)
				thisTheta <- rep(0, groupSize)
				theta <- 0
			} else {
				thisSigma = ( runif(groupSize, min=min(sigmaRawList), max=max(sigmaRawList)) + variationSigmaRaw * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
				thisTheta = thetaList + variationTheta * rt(groupSize, df = 14, ncp = 0)
			}
			thisAnnealing = 0

			for(t in 1:horizon){
				# each individual chooses one option based on his/her choice probability
				choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
				# each subject earns some money (if lucky)
				payoffs[,t] = payoffGenerateGaussian(groupSize, choices[,t], mu_sure, mu_risky, sd_sure, sd_risky)
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
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# In the original article, March only considered a proportional choice
					# If you want to implement softmax rule, you should modify this
					#proportionalChoiceMatrix = Q[,t+1,] %>% apply(1, divideVector, denominator = apply(Q[,t+1,],2,sum)) %>% t()

					###############
					## Softmax choice base solely on Q values
					###############
					#Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta = (thisBeta + thisAnnealing * (t+1)/horizon) ) %>% t() %>% apply(2,expCeiling)
					Q_exp = ( Q[,t+1,] * rep((thisBeta + thisAnnealing * (t+1)/horizon), each = numOptions) ) %>% apply(2,expCeiling)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
					## The followings update the choice probability matrix
					###############
					## Softmax -- END
					###############

					#soc = 1/(1+exp(-(soc_raw)))
					#soc = 0 # soc = 0 indicates asocial learning (i.e. no social info use)
					##netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-epsilon*numOptions)) %>% t() + epsilon
					#netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
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

			safeChoiceProbMean = safeChoiceProb[,round(horizon/2+1):horizon] %>% apply(MARGIN=1, mean)
			# safeChoiceProbMean = safeChoiceProb[,(horizon/2+1):horizon] %>% mean()

			# Submitting this repetition's result
			print(
				data.frame(
					indivOrGroup = rep(indivOrGroup, groupSize),
					variation_level = variation,
					sub = 1:groupSize,
					alpha = thisAlpha,
					beta = thisBeta,
					hot_stove_susceptibility = thisHotStoveSusceptibility,
					sigma = thisSigma,
					theta = thisTheta,
					#mu_theta = rep(theta, groupSize),
					safeChoiceProbMean = safeChoiceProbMean,
					rep = rep
					)
			)
		}
		if(indivOrGroup=="Individual") sigmaFlag <- 1
		gc();gc() # rubbish collection
	}
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
heterogeneous_groups_denrell2007task_beta_data = heterogeneous_groups_denrell2007task_beta[[paste("indivOrGroup=",conditions[1])]][[paste("variation_level=",0)]] %>% data.frame()

for(variation in variation_level_list) {
	#for(sigma in sigmaList) {
		#for (groupSize in groupSizeList) {
			#for (alpha in alphaList) {
				#for (invTemperature in invTemperatureList) {
					#if(variation!=variation_level_list[1]) {
						heterogeneous_groups_denrell2007task_beta_data = heterogeneous_groups_denrell2007task_beta_data %>% rbind(data.frame(heterogeneous_groups_denrell2007task_beta[[paste("indivOrGroup=",'Group')]][[paste("variation_level=",variation)]]))
					#}
				#}
			#}
		#}
	#}
}

#names(heterogeneous_groups_denrell2007task_beta_data) = c('groupSize', 'learningRate', 'invTemperature', 'copyRate', 'conformityExp', 'proportionSafeChoice')

heterogeneous_groups_denrell2007task_beta_data$safeChoiceProbMean_raw <- heterogeneous_groups_denrell2007task_beta_data$safeChoiceProbMean %>% convert_alpha_to_alphaRaw()

heterogeneous_groups_denrell2007task_beta_summary <- heterogeneous_groups_denrell2007task_beta_data %>%
	group_by(indivOrGroup, variation_level, hot_stove_susceptibility) %>%
	summarise(
		#mean_proportionSafeChoice = mean(safeChoiceProbMean, na.rm = TRUE),
	    #sd_proportionSafeChoice = sd(safeChoiceProbMean, na.rm = TRUE),
	    mean_proportionSafeChoice = median(safeChoiceProbMean_raw, na.rm = TRUE) %>% convert_alphaRaw_to_alpha(),
	    median_proportionSafeChoice_raw = median(safeChoiceProbMean_raw, na.rm = TRUE),
	    sd_proportionSafeChoice_raw = sd(safeChoiceProbMean_raw, na.rm = TRUE)
		)

heterogeneous_groups_denrell2007task_beta_summary$hot_stove_susceptibility_slided <- heterogeneous_groups_denrell2007task_beta_summary$hot_stove_susceptibility + rnorm(nrow(heterogeneous_groups_denrell2007task_beta_summary), 0, 0.01)

heterogeneous_groups_denrell2007task_beta_summary_global <-
	heterogeneous_groups_denrell2007task_beta_summary %>%
	group_by(indivOrGroup, variation_level) %>%
	summarise(
		mean_hot_stove_susceptibility = mean(hot_stove_susceptibility, na.rm=TRUE),
		mean_safeChoiceProb = mean(median_proportionSafeChoice_raw, na.rm=TRUE) %>% convert_alphaRaw_to_alpha())

Pr_when_alpha = function (X, alpha) {
	Z = (alpha*(X/alpha - 1)^2)/(2*(2-alpha)) - (X/alpha - 1)/2
	# Z = -beta/(2*(X/(beta+1))-2) * (X - 2)
	return_vector <- 1 / (1 + exp(Z))
	return_vector[which(X < alpha)] <- NA
	return(return_vector)
}

(heterogeneous_groups_denrell2007task_beta_summary %>%
	dplyr::filter(indivOrGroup=='Group') %>%
	ggplot()+
	geom_vline(xintercept=2, linetype='dashed', colour='grey60')+
	geom_hline(yintercept=0.5, linetype='dashed', colour='grey60')+
	stat_function(data=data.frame(X=c(1,8), indivOrGroup='Individual'), fun=Pr_when_alpha, args=list(alpha=0.5), n=1001) +
	geom_point(aes(hot_stove_susceptibility, 1- mean_proportionSafeChoice, colour=as.factor(variation_level)),alpha=3/3)+
	geom_line(aes(hot_stove_susceptibility, 1-convert_alphaRaw_to_alpha(median_proportionSafeChoice_raw), group=variation_level, colour=as.factor(variation_level)), alpha=3/3)+
	geom_point(data=heterogeneous_groups_denrell2007task_beta_summary_global%>%dplyr::filter(indivOrGroup=="Group"), aes(mean_hot_stove_susceptibility, 1-mean_safeChoiceProb, colour=as.factor(variation_level)), size=4, shape=18, alpha = 2/3)+
	# individual condition
	geom_point(data = heterogeneous_groups_denrell2007task_beta_summary%>%filter(variation_level==0), aes(hot_stove_susceptibility, 1- mean_proportionSafeChoice), colour='grey30',alpha=3/3)+
	geom_line(data = heterogeneous_groups_denrell2007task_beta_summary%>%filter(variation_level==0), aes(hot_stove_susceptibility, 1-mean_proportionSafeChoice, group=variation_level), colour='grey30', alpha=3/3)+
	# geom_point(data=heterogeneous_groups_denrell2007task_beta_summary_global%>%filter(variation_level==0), aes(mean_hot_stove_susceptibility, 1-mean_safeChoiceProb), colour='grey30', size=4, shape=18, alpha = 2/3)+
	labs(x=expression(paste("Hot stove susceptibility ",alpha[i],"", (beta[i]+1) ,sep="")),
		y="Proportion of choosing\nthe risky option",
		colour = "Individual\nheterogeneity",
		title=expression(paste("Heterogeneous ", beta[i], sep=""))
		)+
	ylim(c(0,1))+
	xlim(c(0,8))+
	myTheme_Helvetica()+
	scale_colour_viridis_d(direction=-1)+
	theme(legend.position = 'NaN')+
	theme(plot.title = element_text(vjust = - 10, hjust = 0.7))+
	NULL -> heterogeneous_groups_denrell2007task_betaEffect_plot)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_betaEffect_plot.png", plot = heterogeneous_groups_denrell2007task_betaEffect_plot, dpi = 300, width = 4.5, height = 3)

write.csv(heterogeneous_groups_denrell2007task_beta_summary,
			"~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_beta_summary.csv",
			row.names=FALSE)








# Individual variations in the susceptibility to the hot stove effect -- sigma
# where group members are all susceptible to the hot stove effect
repetition = 20000
conditions = c("Individual","Group")
## -- Global parameter setting --
# sigmaList = 0.3
# sigmaRawList = list(
# 	'1'=rep(-0.8472979, groupSize),
# 	'2'=c(-1.371583, -1.10944, -0.8472979, -0.5851553, -0.3230127),
# 	'3'=c(-1.895868, -1.371583, -0.8472979, -0.3230127, 0.2012725),
# 	'4'=c(-2.944439, -1.895868, -0.8472979, 0.2012725, 1.249843),
# 	'5'=c(-5.8472979, -2.944439, -0.8472979, 1.249843, 4.152702)
# 	)
sigmaList = list(
	'1'=rep(0.3, groupSize),
	# '2'=c(0.225, 0.225, 0.3, 0.45, 0.45),
	'3'=c(0.15, 0.225, 0.3, 0.375, 0.45),
	# '4'=c(0.05, 0.15, 0.3, 0.45, 0.55),
	'5'=c(0.05, 0.1, 0.3, 0.5, 0.55),
	# '6'=c(0.01, 0.1, 0.3, 0.5, 0.59),
	'7'=c(0.01, 0.04, 0.1, 0.55, 0.8)
	)
# thetaList = c(1, 2, 4)
thetaList = 2
# alphaList = list(
# 	'1'=rep(0.55, groupSize),
# 	'2'=c(0.45, 0.50, 0.55, 0.60, 0.65),
# 	'3'=c(0.35, 0.45, 0.55, 0.65, 0.75),
# 	'4'=c(0.2, 0.32, 0.55, 0.78, 0.9)
# 	)
alphaList = 0.5
invTemperatureList = 7
# invTemperatureList = list(
# 	'1'=rep(7, groupSize),
# 	'2'=c(5.75, 6.375, 7, 7.625, 8.25),
# 	'3'=c(4.5, 5.75, 7, 8.25, 9.5),
# 	'4'=c(2, 4.5, 7, 9.5, 12)
# 	)
#hot_stove_susceptibility = (invTemperatureList+1) * alphaList
variation_level_list = 1:length(sigmaList)

## -- simulation parameters
groupSizeList = c(5)
horizon = 150 # = number of trials
numOptions = 2

## -- transformed parameters
#alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList)
# sigmaRawList = mapply(convert_alpha_to_alphaRaw, sigmaList)

## -- Individual variation
variationAlphaRaw = 0.01 #1.5
variationBeta = 0.01#0.7
variatonAnnealing = 0.01#1.5
variationSigmaRaw = 0.01
variationTheta = 0.01

# -- Setting a two-armed bandit task --
riskPremium = 1.5
mu_sure = 1
mu_risky = mu_sure * riskPremium
sd_sure = 0.05
sd_risky = 1
initialExpextation = 0

# -- Initial Settings --
heterogeneous_groups_denrell2007task_sigma = list()
groupSize = groupSizeList[1]

sigmaFlag = 0
# -- simulation --
s_time = Sys.time()
for (indivOrGroup in conditions) {
	if (sigmaFlag==1 & indivOrGroup=="Individual") next
	for (variation in variation_level_list) {
		heterogeneous_groups_denrell2007task_sigma[[paste("indivOrGroup=",indivOrGroup)]][[paste("variation_level=",variation)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
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
			alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList[1][[1]])
			thisAlpha = ( alphaRawList + 0 * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
			thisBeta = invTemperatureList[1] + 0  * rt(groupSize, df = 14, ncp = 0)
			thisBeta[which(thisBeta<0)] <- 0
			thisHotStoveSusceptibility <- thisAlpha * (thisBeta+1)
			#thisAnnealing = annealing + variatonAnnealing * rt(groupSize, df = 14, ncp = 0)
			if(indivOrGroup == 'Individual'){
				thisSigma <- rep(0, groupSize)
				thisTheta <- rep(0, groupSize)
				theta <- 0
			} else {
				sigmaRawList = mapply(convert_alpha_to_alphaRaw, sigmaList[variation][[1]]) %>% sample()
				thisSigma = ( sigmaRawList + 0 * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
				thisTheta = thetaList + 0 * rt(groupSize, df = 14, ncp = 0)
			}
			thisAnnealing = 0

			for(t in 1:horizon){
				# each individual chooses one option based on his/her choice probability
				choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
				# each subject earns some money (if lucky)
				payoffs[,t] = payoffGenerateGaussian(groupSize, choices[,t], mu_sure, mu_risky, sd_sure, sd_risky)
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
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# In the original article, March only considered a proportional choice
					# If you want to implement softmax rule, you should modify this
					#proportionalChoiceMatrix = Q[,t+1,] %>% apply(1, divideVector, denominator = apply(Q[,t+1,],2,sum)) %>% t()

					###############
					## Softmax choice base solely on Q values
					###############
					#Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta = (thisBeta + thisAnnealing * (t+1)/horizon) ) %>% t() %>% apply(2,expCeiling)
					Q_exp = ( Q[,t+1,] * rep((thisBeta + thisAnnealing * (t+1)/horizon), each = numOptions) ) %>% apply(2,expCeiling)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
					## The followings update the choice probability matrix
					###############
					## Softmax -- END
					###############

					#soc = 1/(1+exp(-(soc_raw)))
					#soc = 0 # soc = 0 indicates asocial learning (i.e. no social info use)
					##netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-epsilon*numOptions)) %>% t() + epsilon
					#netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
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

			safeChoiceProbMean = safeChoiceProb[,round(horizon/2+1):horizon] %>% apply(MARGIN=1, mean)
			# safeChoiceProbMean = safeChoiceProb[,(horizon/2+1):horizon] %>% mean()

			# Submitting this repetition's result
			print(
				data.frame(
					indivOrGroup = rep(indivOrGroup, groupSize),
					variation_level = variation,
					sub = 1:groupSize,
					alpha = thisAlpha,
					beta = thisBeta,
					hot_stove_susceptibility = thisHotStoveSusceptibility,
					sigma = thisSigma,
					theta = thisTheta,
					#mu_theta = rep(theta, groupSize),
					safeChoiceProbMean = safeChoiceProbMean,
					rep = rep
					)
			)
		}
		if(indivOrGroup=="Individual") sigmaFlag <- 1
		gc();gc() # rubbish collection
	}
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
heterogeneous_groups_denrell2007task_sigma_data = heterogeneous_groups_denrell2007task_sigma[[paste("indivOrGroup=",conditions[1])]][[paste("variation_level=",variation_level_list[1])]] %>% data.frame()

for(variation in variation_level_list) {
	#for(sigma in sigmaList) {
		#for (groupSize in groupSizeList) {
			#for (alpha in alphaList) {
				#for (invTemperature in invTemperatureList) {
					#if(variation!=variation_level_list[1]) {
						heterogeneous_groups_denrell2007task_sigma_data = heterogeneous_groups_denrell2007task_sigma_data %>% rbind(data.frame(heterogeneous_groups_denrell2007task_sigma[[paste("indivOrGroup=",'Group')]][[paste("variation_level=",variation)]]))
					#}
				#}
			#}
		#}
	#}
}

#names(heterogeneous_groups_denrell2007task_sigma_data) = c('groupSize', 'learningRate', 'invTemperature', 'copyRate', 'conformityExp', 'proportionSafeChoice')

heterogeneous_groups_denrell2007task_sigma_data$safeChoiceProbMean_raw <- heterogeneous_groups_denrell2007task_sigma_data$safeChoiceProbMean %>% convert_alpha_to_alphaRaw()

heterogeneous_groups_denrell2007task_sigma_summary <- heterogeneous_groups_denrell2007task_sigma_data %>%
	group_by(indivOrGroup, variation_level, hot_stove_susceptibility, sigma, theta) %>%
	summarise(
		#mean_proportionSafeChoice = mean(safeChoiceProbMean, na.rm = TRUE),
	    #sd_proportionSafeChoice = sd(safeChoiceProbMean, na.rm = TRUE),
	    mean_proportionSafeChoice = median(safeChoiceProbMean_raw, na.rm = TRUE) %>% convert_alphaRaw_to_alpha(),
	    median_proportionSafeChoice_raw = median(safeChoiceProbMean_raw, na.rm = TRUE),
	    sd_proportionSafeChoice_raw = sd(safeChoiceProbMean_raw, na.rm = TRUE)
		)

heterogeneous_groups_denrell2007task_sigma_summary$hot_stove_susceptibility_slided <- heterogeneous_groups_denrell2007task_sigma_summary$hot_stove_susceptibility + rnorm(nrow(heterogeneous_groups_denrell2007task_sigma_summary), 0, 0.01)

heterogeneous_groups_denrell2007task_sigma_summary_global <-
	heterogeneous_groups_denrell2007task_sigma_summary %>%
	group_by(indivOrGroup, variation_level) %>%
	summarise(
		mean_sigma = mean(sigma, na.rm=TRUE),
		mean_safeChoiceProb = mean(median_proportionSafeChoice_raw, na.rm=TRUE) %>% convert_alphaRaw_to_alpha())

Pr_when_beta = function (X, beta) {
	Z = -beta/(2*(X/(beta+1))-2) * (X - 2)
	return_vector <- 1 / (1 + exp(Z))
	return_vector[which(X > beta)] <- NA
	return(return_vector)
}

(heterogeneous_groups_denrell2007task_sigma_summary %>%
	dplyr::filter(indivOrGroup=="Group") %>%
	ggplot()+
	geom_hline(yintercept=0.5, linetype='dashed', colour='grey60')+
	#stat_function(data=data.frame(X=c(2,9), invTemperature_factor='β = 7', indivOrGroup='Individual'), fun=Pr_when_beta, args=list(beta=7)) +
	geom_point(aes(sigma, 1- mean_proportionSafeChoice, colour=as.factor(variation_level)),alpha=3/3)+
	geom_line(aes(sigma, 1- mean_proportionSafeChoice, group=variation_level, colour=as.factor(variation_level)), alpha=3/3)+
	geom_point(data=heterogeneous_groups_denrell2007task_sigma_summary_global%>%dplyr::filter(indivOrGroup=='Group'), aes(mean_sigma, 1-mean_safeChoiceProb, colour=as.factor(variation_level)), size=4, shape=18, alpha = 2/3)+
	geom_point(data=heterogeneous_groups_denrell2007task_sigma_summary_global%>%dplyr::filter(indivOrGroup=='Individual'), aes(mean_sigma, 1-mean_safeChoiceProb), colour='black', size=4, shape=18, alpha = 2/3)+
	labs(x=expression(paste("Social learning weight ", sigma[i], sep="")),
		y="Proportion of choosing\nthe risky option",
		colour = "Individual\nheterogeneity",
		title=expression(paste("Heterogeneous ", sigma[i], sep=""))
		)+
	ylim(c(0,1))+
	xlim(c(0,1))+
	myTheme_Helvetica()+
	scale_colour_viridis_d(direction=-1)+
	theme(legend.position = 'NaN')+
	theme(plot.title = element_text(vjust = - 10))+
	NULL -> heterogeneous_groups_denrell2007task_sigmaEffect_plot)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_sigmaEffect_plot.png", plot = heterogeneous_groups_denrell2007task_sigmaEffect_plot, dpi = 300, width = 4.5, height = 3)

write.csv(heterogeneous_groups_denrell2007task_sigma_summary,
			"~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_sigma_summary.csv",
			row.names=FALSE)







# Individual variations in the susceptibility to the hot stove effect -- theta
# where group members are all susceptible to the hot stove effect
repetition = 20000
conditions = c("Individual","Group")
## -- Global parameter setting --
sigmaList = 0.3
# sigmaRawList = list(
# 	'1'=rep(-0.8472979, groupSize),
# 	'2'=c(-1.371583, -1.10944, -0.8472979, -0.5851553, -0.3230127),
# 	'3'=c(-1.895868, -1.371583, -0.8472979, -0.3230127, 0.2012725),
# 	'4'=c(-2.944439, -1.895868, -0.8472979, 0.2012725, 1.249843),
# 	'5'=c(-5.8472979, -2.944439, -0.8472979, 1.249843, 4.152702)
# 	)
# sigmaList = list(
# 	'1'=rep(0.3, groupSize),
# 	'2'=c(0.225, 0.225, 0.3, 0.45, 0.45),
# 	'3'=c(0.15, 0.225, 0.3, 0.375, 0.45),
# 	'4'=c(0.05, 0.15, 0.3, 0.45, 0.55),
# 	'5'=c(0.05, 0.1, 0.3, 0.5, 0.55),
# 	'6'=c(0.01, 0.1, 0.3, 0.5, 0.59),
# 	'7'=c(0.01, 0.04, 0.1, 0.55, 0.8)
# 	)
# thetaList = c(1, 2, 4)
# thetaList = 2
thetaList = list(
	'1'=rep(2, groupSize),
	# '2'=c(1.75, 1.825, 2, 2.125, 2.25),
	'2'=c(1.5, 1.75, 2, 2.25, 2.5),
	# '3'=c(1, 1.5, 2, 2.5, 3),
	# '3'=c(0.5, 1, 2, 3, 3.5),
	#'3'=c(0, 0.5, 2, 3.5, 4),
	'4'=c(0, 1.8, 2.2, 2.8, 3.2), # skewed to conformist
	'5'=c(0, 0.5, 0.8, 1.2, 7.5), # skewed to non-conformist
	'6'=c(-1, 0.5, 2, 3.5, 5) # skewed to non-conformist
	)
# alphaList = list(
# 	'1'=rep(0.55, groupSize),
# 	'2'=c(0.45, 0.50, 0.55, 0.60, 0.65),
# 	'3'=c(0.35, 0.45, 0.55, 0.65, 0.75),
# 	'4'=c(0.2, 0.32, 0.55, 0.78, 0.9)
# 	)
alphaList = 0.5
invTemperatureList = 7
# invTemperatureList = list(
# 	'1'=rep(7, groupSize),
# 	'2'=c(5.75, 6.375, 7, 7.625, 8.25),
# 	'3'=c(4.5, 5.75, 7, 8.25, 9.5),
# 	'4'=c(2, 4.5, 7, 9.5, 12)
# 	)
#hot_stove_susceptibility = (invTemperatureList+1) * alphaList
variation_level_list = 1:length(thetaList)

## -- simulation parameters
groupSizeList = c(5)
horizon = 150 # = number of trials
numOptions = 2

## -- transformed parameters
#alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList)
# sigmaRawList = mapply(convert_alpha_to_alphaRaw, sigmaList)

## -- Individual variation
variationAlphaRaw = 0.01 #1.5
variationBeta = 0.01#0.7
variatonAnnealing = 0.01#1.5
variationSigmaRaw = 0.01
variationTheta = 0.01

# -- Setting a two-armed bandit task --
riskPremium = 1.5
mu_sure = 1
mu_risky = mu_sure * riskPremium
sd_sure = 0.05
sd_risky = 1
initialExpextation = 0

# -- Initial Settings --
heterogeneous_groups_denrell2007task_theta = list()
groupSize = groupSizeList[1]

sigmaFlag = 0
# -- simulation --
s_time = Sys.time()
for (indivOrGroup in conditions) {
	if (sigmaFlag==1 & indivOrGroup=="Individual") next
	for (variation in variation_level_list) {
		heterogeneous_groups_denrell2007task_theta[[paste("indivOrGroup=",indivOrGroup)]][[paste("variation_level=",variation)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
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
			alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList[1][[1]])
			thisAlpha = ( alphaRawList + 0 * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
			thisBeta = invTemperatureList[1] + 0  * rt(groupSize, df = 14, ncp = 0)
			thisBeta[which(thisBeta<0)] <- 0
			thisHotStoveSusceptibility <- thisAlpha * (thisBeta+1)
			#thisAnnealing = annealing + variatonAnnealing * rt(groupSize, df = 14, ncp = 0)
			if(indivOrGroup == 'Individual'){
				thisSigma <- rep(0, groupSize)
				thisTheta <- rep(0, groupSize)
				theta <- 0
			} else {
				sigmaRawList = mapply(convert_alpha_to_alphaRaw, sigmaList[1][[1]])
				thisSigma = ( sigmaRawList + 0 * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
				thisTheta = (thetaList[variation][[1]] + 0 * rt(groupSize, df = 14, ncp = 0)) %>% sample()
			}
			thisAnnealing = 0

			for(t in 1:horizon){
				# each individual chooses one option based on his/her choice probability
				choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
				# each subject earns some money (if lucky)
				payoffs[,t] = payoffGenerateGaussian(groupSize, choices[,t], mu_sure, mu_risky, sd_sure, sd_risky)
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
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# In the original article, March only considered a proportional choice
					# If you want to implement softmax rule, you should modify this
					#proportionalChoiceMatrix = Q[,t+1,] %>% apply(1, divideVector, denominator = apply(Q[,t+1,],2,sum)) %>% t()

					###############
					## Softmax choice base solely on Q values
					###############
					#Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta = (thisBeta + thisAnnealing * (t+1)/horizon) ) %>% t() %>% apply(2,expCeiling)
					Q_exp = ( Q[,t+1,] * rep((thisBeta + thisAnnealing * (t+1)/horizon), each = numOptions) ) %>% apply(2,expCeiling)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
					## The followings update the choice probability matrix
					###############
					## Softmax -- END
					###############

					#soc = 1/(1+exp(-(soc_raw)))
					#soc = 0 # soc = 0 indicates asocial learning (i.e. no social info use)
					##netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-epsilon*numOptions)) %>% t() + epsilon
					#netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
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

			safeChoiceProbMean = safeChoiceProb[,round(horizon/2+1):horizon] %>% apply(MARGIN=1, mean)
			# safeChoiceProbMean = safeChoiceProb[,(horizon/2+1):horizon] %>% mean()

			# Submitting this repetition's result
			print(
				data.frame(
					indivOrGroup = rep(indivOrGroup, groupSize),
					variation_level = variation,
					sub = 1:groupSize,
					alpha = thisAlpha,
					beta = thisBeta,
					hot_stove_susceptibility = thisHotStoveSusceptibility,
					sigma = thisSigma,
					theta = thisTheta,
					#mu_theta = rep(theta, groupSize),
					safeChoiceProbMean = safeChoiceProbMean,
					rep = rep
					)
			)
		}
		if(indivOrGroup=="Individual") sigmaFlag <- 1
		gc();gc() # rubbish collection
	}
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
heterogeneous_groups_denrell2007task_theta_data = heterogeneous_groups_denrell2007task_theta[[paste("indivOrGroup=",conditions[1])]][[paste("variation_level=",variation_level_list[1])]] %>% data.frame()

for(variation in variation_level_list) {
	#for(sigma in sigmaList) {
		#for (groupSize in groupSizeList) {
			#for (alpha in alphaList) {
				#for (invTemperature in invTemperatureList) {
					#if(variation!=variation_level_list[1]) {
						heterogeneous_groups_denrell2007task_theta_data = heterogeneous_groups_denrell2007task_theta_data %>% rbind(data.frame(heterogeneous_groups_denrell2007task_theta[[paste("indivOrGroup=",'Group')]][[paste("variation_level=",variation)]]))
					#}
				#}
			#}
		#}
	#}
}

#names(heterogeneous_groups_denrell2007task_theta_data) = c('groupSize', 'learningRate', 'invTemperature', 'copyRate', 'conformityExp', 'proportionSafeChoice')

heterogeneous_groups_denrell2007task_theta_data$safeChoiceProbMean_raw <- heterogeneous_groups_denrell2007task_theta_data$safeChoiceProbMean %>% convert_alpha_to_alphaRaw()

heterogeneous_groups_denrell2007task_theta_summary <- heterogeneous_groups_denrell2007task_theta_data %>%
	group_by(indivOrGroup, variation_level, hot_stove_susceptibility, sigma, theta) %>%
	summarise(
		#mean_proportionSafeChoice = mean(safeChoiceProbMean, na.rm = TRUE),
	    #sd_proportionSafeChoice = sd(safeChoiceProbMean, na.rm = TRUE),
	    mean_proportionSafeChoice = median(safeChoiceProbMean_raw, na.rm = TRUE) %>% convert_alphaRaw_to_alpha(),
	    median_proportionSafeChoice_raw = median(safeChoiceProbMean_raw, na.rm = TRUE),
	    sd_proportionSafeChoice_raw = sd(safeChoiceProbMean_raw, na.rm = TRUE)
		)

heterogeneous_groups_denrell2007task_theta_summary$hot_stove_susceptibility_slided <- heterogeneous_groups_denrell2007task_theta_summary$hot_stove_susceptibility + rnorm(nrow(heterogeneous_groups_denrell2007task_theta_summary), 0, 0.01)

heterogeneous_groups_denrell2007task_theta_summary_global <-
	heterogeneous_groups_denrell2007task_theta_summary %>%
	group_by(indivOrGroup, variation_level) %>%
	summarise(
		mean_theta = mean(theta, na.rm=TRUE),
		mean_safeChoiceProb = mean(median_proportionSafeChoice_raw, na.rm=TRUE) %>% convert_alphaRaw_to_alpha())

Pr_when_beta = function (X, beta) {
	Z = -beta/(2*(X/(beta+1))-2) * (X - 2)
	return_vector <- 1 / (1 + exp(Z))
	return_vector[which(X > beta)] <- NA
	return(return_vector)
}

(heterogeneous_groups_denrell2007task_theta_summary %>%
	dplyr::filter(indivOrGroup=="Group") %>%
	ggplot()+
	geom_hline(yintercept=0.5, linetype='dashed', colour='grey60')+
	#stat_function(data=data.frame(X=c(2,9), invTemperature_factor='β = 7', indivOrGroup='Individual'), fun=Pr_when_beta, args=list(beta=7)) +
	geom_point(aes(theta, 1- mean_proportionSafeChoice, colour=as.factor(variation_level)),alpha=3/3)+
	geom_line(aes(theta, 1- mean_proportionSafeChoice, group=variation_level, colour=as.factor(variation_level)), alpha=3/3)+
	geom_point(data=heterogeneous_groups_denrell2007task_theta_summary_global%>%dplyr::filter(indivOrGroup=='Group'), aes(mean_theta, 1-mean_safeChoiceProb, colour=as.factor(variation_level)), size=4, shape=18, alpha = 2/3)+
	geom_point(data=heterogeneous_groups_denrell2007task_theta_summary_global%>%dplyr::filter(indivOrGroup=='Individual'), aes(mean_theta, 1-mean_safeChoiceProb), colour='black', size=4, shape=18, alpha = 2/3)+
	labs(x=expression(paste("Conformity exponent ", theta[i], sep="")),
		y="Proportion of choosing\nthe risky option",
		colour = "Individual\nheterogeneity",
		title=expression(paste("Heterogeneous ", theta[i], sep=""))
		)+
	ylim(c(0,1))+
	xlim(c(-1,8))+
	myTheme_Helvetica()+
	scale_colour_viridis_d(direction=-1)+
	theme(legend.position = 'NaN')+
	theme(plot.title = element_text(vjust = - 10))+
	NULL -> heterogeneous_groups_denrell2007task_thetaEffect_plot)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_thetaE ffect_plot.png", plot = heterogeneous_groups_denrell2007task_thetaEffect_plot, dpi = 300, width = 4.5, height = 3)

write.csv(heterogeneous_groups_denrell2007task_theta_summary,
			"~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_theta_summary.csv",
			row.names=FALSE)

heteroeneity_plot <- plot_grid(
	heterogeneous_groups_denrell2007task_alphaEffect_plot,
	heterogeneous_groups_denrell2007task_betaEffect_plot,
	heterogeneous_groups_denrell2007task_sigmaEffect_plot,
	heterogeneous_groups_denrell2007task_thetaEffect_plot,
	 labels = c('','','',''), ncol = 2, align = 'v')
ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heteroeneity_plot.pdf", plot = heteroeneity_plot, dpi = 600, width = 9, height = 6)








































## -- Global parameter setting --
conditions = c("Individual", "Group")
sigmaList = seq(0.05, 0.4, 0.1)
thetaList = c(1, 2, 4)
alphaList = c(0.1, 0.2, 0.4, 0.6, 0.8)
invTemperatureList = c(7)
hot_stove_susceptibility = (invTemperatureList+1) * alphaList
annealing = 0

## -- simulation parameters
groupSizeList = c(5)
repetition = 20  #1000
horizon = 150 # = number of trials
numOptions = 2

## -- transformed parameters
alphaRawList = mapply(convert_alpha_to_alphaRaw, alphaList)
sigmaRawList = mapply(convert_alpha_to_alphaRaw, sigmaList)

## -- Individual variation
variationAlphaRaw = 0 #1.5
variationBeta = 0#0.7
variatonAnnealing = 0#1.5
variationSigmaRaw = 0
variationTheta = 0

# -- Setting a two-armed bandit task --
riskPremium = 1.5
mu_sure = 1
mu_risky = mu_sure * riskPremium
sd_sure = 0.05
sd_risky = 1
initialExpextation = 0

# -- Initial Settings --
heterogeneous_groups_denrell2007task_vol1 = list()
groupSize = groupSizeList[1]

sigmaFlag = 0
# -- simulation --
s_time = Sys.time()
for (indivOrGroup in conditions) {
	if (sigmaFlag==1 & indivOrGroup=="Individual") next
	for (theta in thetaList) {
		heterogeneous_groups_denrell2007task_vol1[[paste("indivOrGroup=",indivOrGroup)]][[paste("theta=", theta)]] <- foreach(rep = 1:repetition, .combine=rbind) %dopar% {
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
			thisAlpha = ( alphaRawList + 0 * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
			thisBeta = invTemperatureList + 0  * rt(groupSize, df = 14, ncp = 0)
			thisBeta[which(thisBeta<0)] <- 0
			thisHotStoveSusceptibility <- thisAlpha * (thisBeta+1)
			#thisAnnealing = annealing + variatonAnnealing * rt(groupSize, df = 14, ncp = 0)
			if(indivOrGroup == 'Individual'){
				thisSigma <- rep(0, groupSize)
				thisTheta <- rep(0, groupSize)
				theta <- 0
			} else {
				thisSigma = ( runif(groupSize, min=min(sigmaRawList), max=max(sigmaRawList)) + variationSigmaRaw * rt(groupSize, df = 14, ncp = 0) ) %>% mapply(convert_alphaRaw_to_alpha, .)
				thisTheta = theta + variationTheta * rt(groupSize, df = 14, ncp = 0)
			}
			thisAnnealing = 0
			for(t in 1:horizon){
				# each individual chooses one option based on his/her choice probability
				choices[,t] = mapply(function(p1,p2){ sample(1:numOptions, 1, prob=c(p1,p2), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,] )
				# each subject earns some money (if lucky)
				payoffs[,t] = payoffGenerateGaussian(groupSize, choices[,t], mu_sure, mu_risky, sd_sure, sd_risky)
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
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# In the original article, March only considered a proportional choice
					# If you want to implement softmax rule, you should modify this
					#proportionalChoiceMatrix = Q[,t+1,] %>% apply(1, divideVector, denominator = apply(Q[,t+1,],2,sum)) %>% t()

					###############
					## Softmax choice base solely on Q values
					###############
					#Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta = (thisBeta + thisAnnealing * (t+1)/horizon) ) %>% t() %>% apply(2,expCeiling)
					Q_exp = ( Q[,t+1,] * rep((thisBeta + thisAnnealing * (t+1)/horizon), each = numOptions) ) %>% apply(2,expCeiling)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], thisTheta, numOptions)
					## The followings update the choice probability matrix
					###############
					## Softmax -- END
					###############

					#soc = 1/(1+exp(-(soc_raw)))
					#soc = 0 # soc = 0 indicates asocial learning (i.e. no social info use)
					##netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-epsilon*numOptions)) %>% t() + epsilon
					#netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
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

			safeChoiceProbMean = safeChoiceProb[,round(2*horizon/3+1):horizon] %>% apply(MARGIN=1, mean)
			# safeChoiceProbMean = safeChoiceProb[,(horizon/2+1):horizon] %>% mean()

			# Submitting this repetition's result
			print(
				data.frame(
					indivOrGroup = rep(indivOrGroup, groupSize),
					sub = 1:groupSize,
					alpha = thisAlpha,
					beta = thisBeta,
					hot_stove_susceptibility = thisHotStoveSusceptibility,
					sigma = thisSigma,
					theta = thisTheta,
					mu_theta = rep(theta, groupSize),
					safeChoiceProbMean = safeChoiceProbMean,
					rep = rep
					)
			)
		}
		if(indivOrGroup=="Individual") sigmaFlag <- 1
		gc();gc() # rubbish collection
	}
}
e_time = Sys.time()
e_time - s_time
# -- simulation END --

# -- saving the data --
heterogeneous_groups_denrell2007task_vol1_data = heterogeneous_groups_denrell2007task_vol1[[paste("indivOrGroup=",'Individual')]][[paste("theta=",thetaList[1])]] %>% data.frame()

for(theta in thetaList) {
	#for(sigma in sigmaList) {
		#for (groupSize in groupSizeList) {
			#for (alpha in alphaList) {
				#for (invTemperature in invTemperatureList) {
					#if(theta!=thetaList[1]) {
						heterogeneous_groups_denrell2007task_vol1_data = heterogeneous_groups_denrell2007task_vol1_data %>% rbind(data.frame(heterogeneous_groups_denrell2007task_vol1[[paste("indivOrGroup=",'Group')]][[paste("theta=",theta)]]))
					#}
				#}
			#}
		#}
	#}
}

#names(heterogeneous_groups_denrell2007task_vol1_data) = c('groupSize', 'learningRate', 'invTemperature', 'copyRate', 'conformityExp', 'proportionSafeChoice')

heterogeneous_groups_denrell2007task_vol1_data$safeChoiceProbMean_raw <- heterogeneous_groups_denrell2007task_vol1_data$safeChoiceProbMean %>% convert_alpha_to_alphaRaw()

heterogeneous_groups_denrell2007task_vol1_summary <- heterogeneous_groups_denrell2007task_vol1_data %>%
	group_by(indivOrGroup, sub, hot_stove_susceptibility, mu_theta) %>%
	summarise(
		mean_proportionSafeChoice = mean(safeChoiceProbMean, na.rm = TRUE),
	    sd_proportionSafeChoice = sd(safeChoiceProbMean, na.rm = TRUE),
	    median_proportionSafeChoice_raw = median(safeChoiceProbMean_raw, na.rm = TRUE),
	    sd_proportionSafeChoice_raw = sd(safeChoiceProbMean_raw, na.rm = TRUE)
		)

heterogeneous_groups_denrell2007task_vol1_summary$hot_stove_susceptibility_slided <- heterogeneous_groups_denrell2007task_vol1_summary$hot_stove_susceptibility + rnorm(nrow(heterogeneous_groups_denrell2007task_vol1_summary), 0, 0.1)

Pr_when_beta = function (X, beta) {
	Z = -beta/(2*(X/(beta+1))-2) * (X - 2)
	return_vector <- 1 / (1 + exp(Z))
	return_vector[which(X > beta)] <- NA
	return(return_vector)
}

(heterogeneous_groups_denrell2007task_vol1_summary %>%
	ggplot()+
	stat_function(data=data.frame(X=c(2,8), invTemperature_factor='β = 7', mu_theta=0), fun=Pr_when_beta, args=list(beta=7)) +
	geom_point(aes(hot_stove_susceptibility_slided, 1- convert_alphaRaw_to_alpha(median_proportionSafeChoice_raw), colour=as.factor(mu_theta)))+
	geom_line(aes(hot_stove_susceptibility_slided, 1-convert_alphaRaw_to_alpha(median_proportionSafeChoice_raw), group=mu_theta, colour=as.factor(mu_theta)))+
	geom_errorbar(aes(hot_stove_susceptibility_slided, ymin=1- convert_alphaRaw_to_alpha(median_proportionSafeChoice_raw+sd_proportionSafeChoice_raw), ymax=1-convert_alphaRaw_to_alpha(median_proportionSafeChoice_raw-sd_proportionSafeChoice_raw), colour=as.factor(mu_theta)), width=0)+
	geom_hline(yintercept=0.5, linetype='dashed')+
	labs(x=expression(paste("Individual susceptibility",alpha[i],"*", (beta[i]+1) ,sep="")),
		y="Proportion of choosing\nthe risky option",
		colour = "Conformity\nexponent"#expression(theta)
		)+
	ylim(c(0,1))+
	myTheme_Helvetica()+
	scale_colour_viridis_d()+
	NULL -> heterogeneous_groups_denrell2007task_vol1_plot)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_vol1_plot.png", plot = heterogeneous_groups_denrell2007task_vol1_plot, dpi = 300, width = 4.5, height = 3)


write.csv(heterogeneous_groups_denrell2007task_vol1_summary,
			"~/Dropbox/wataru/papers/RiskySocialLearning/riskPrem1.5_Gaussian/results/heterogeneous_groups_denrell2007task_vol1_summary.csv",
			row.names=FALSE)







