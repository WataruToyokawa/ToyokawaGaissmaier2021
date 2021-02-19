###############################################################################################
##
## Figures of the collective rescue paper (ToyokawaGaissmaier2021)
## Wataru Toyokawa
## 07 November. 2020
##
###############################################################################################
library(tidyverse)
library(cowplot)
library(metR)

## Load Functions
# setwd("~/analysis_repo") <------- Set this folder as a working directory
source("agentBasedSim/functions.R")


## ========================================================
#
# figure 1: the collective rescue in the stable environment
#
# =========================================================

# schematic figure
stable_2AB_scheme <- cowplot::ggdraw() + cowplot::draw_image("agentBasedSim/stable_2AB_scheme.png", scale = 1)

# figures of the simulation
socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_data <- read_csv("agentBasedSim/socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_data.csv")

socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_data %>%
	group_by(groupSize, learningRate, invTemperature, copyRate, conformityExp) %>%
	summarise(
		mean_proportionSafeChoice = mean(proportionSafeChoice, na.rm = TRUE),
	    median_proportionSafeChoice = median(proportionSafeChoice, na.rm = TRUE),
	    sd_proportionSafeChoice = sd(proportionSafeChoice, na.rm = TRUE)
		)


socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$copyRate_factor = paste(rep('\U03C3 = ', nrow(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary)), socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$copyRate, sep ='')

socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$conformityExp_factor = paste(rep('\U03b8 = ', nrow(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary)), socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$conformityExp, sep ='')

socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$invTemperature_factor = paste(rep('\U03b2 = ', nrow(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary)), socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$invTemperature, sep ='')
#socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$conformityExp_factor = factor(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$conformityExp_factor, levels = c('θ = 1','θ = 2','θ = 4','θ = 8'))

## Adding Denrell (2007)'s analytical solution
Denrell2007 = function (alpha, beta, mu, sd) {
	1 / ( 1 + exp( (alpha*(beta^2)*(sd^2))/(2*(2-alpha)) - beta*mu ) )
}

## when mu = surePayoff + 0.5 and sd = 1, the above function can be reduced as follows:
Denrell2007Solution = function (alpha) {
	(2 - alpha)/alpha
}

Denrell2007RiskyChoice = c()
alphaArray = c()
betaArray = c()
for (alpha in seq(0,1,0.1)) {
	for (beta in seq(0,10,1)) {
		alphaArray <- append(alphaArray, alpha)
		betaArray <- append(betaArray, beta)
		Denrell2007RiskyChoice <- append(Denrell2007RiskyChoice, Denrell2007(alpha, beta, mu=0.5, sd=1))
	}
}

Denrell2007Simulation = data.frame(alpha = alphaArray, beta = betaArray, riskyChoiceProb = Denrell2007RiskyChoice)

# (Denrell2007Simulation %>%
# 	ggplot(aes(alpha, beta))+
# 	geom_raster(aes(fill = riskyChoiceProb), stat = 'identity')+
# 	stat_function(fun=Denrell2007Solution, color='black', linetype='dashed', size=2/3)+
# 	scale_fill_gradient2(midpoint = 0.5, low = "blue", mid = "grey90", high = "red")+
# 	labs(x=expression(paste('Learning rate ',alpha,sep="")), y=expression(paste('Inverse temperature ',beta,sep="")), title='', fill="Proportion of \nchoosing \nthe risky option")+
# 	myTheme_Helvetica()+
# 	theme(axis.text.x = element_text(angle = 90))+
# 	theme(strip.text.y = element_text(angle = 0))+
# 	theme(legend.text = element_text(angle = 0))+
# 	#theme(legend.position = 'top')+
# 	ylim(c(0,10))+
# 	NULL -> Denrell2007_figure_analytical
# )

(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>%
	dplyr::filter(copyRate ==0) %>%
	ggplot() +
	geom_raster(mapping = aes(learningRate, invTemperature, fill = 1-mean_proportionSafeChoice), stat = 'identity') +
	labs(x=expression(paste('Learning rate ',alpha,sep="")),
		y=expression(paste('Inverse temperature ',beta,sep="")),
		#title='Gaussian noise\n mu=0.5; sigma=1',
		fill = "Proportion of \nchoosing \nthe risky option"
		#title='Gaussian noise\n mu=0.5; sigma=1', fill = "Proportion of\nsafe choice"
		)+
	#scale_fill_viridis(limits = c(0.45, 1), option="magma")+
	stat_function(fun=Denrell2007Solution, color='black', linetype='dashed', size=2/3)+
	scale_fill_gradient2(midpoint = 0.5, high = "red", mid = "grey90", low = "blue")+
	ylim(c(0,10))+
	#myTheme_gillsansMT()+
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle = 90))+
	theme(strip.text.y = element_text(angle = 0))+
	theme(legend.text = element_text(angle = 0))+
	# theme(legend.position = 'top')+
	#facet_grid(copyRate_factor ~ conformityExp_factor)+
	NULL -> Denrell2007_figure
)

# plot with Denrell (2007) curve
(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>%
	dplyr::filter(copyRate %in% c(0.25, 0.5)) %>%
	dplyr::filter(conformityExp %in% c(1, 4)) %>%
	ggplot() +
	geom_raster(mapping = aes(learningRate, invTemperature, fill = 1-mean_proportionSafeChoice), stat = 'identity') +
	labs(x=expression(paste('Learning rate ',alpha,sep="")),
		y=expression(paste('Inverse temperature ',beta,sep="")),
		#title='Gaussian noise\n mu=0.5; sigma=1',
		fill = "Proportion of \nchoosing \nthe risky option"
		#title='Gaussian noise\n mu=0.5; sigma=1', fill = "Proportion of\nsafe choice"
		)+
	#scale_fill_viridis(limits = c(0.45, 1), option="magma")+
	stat_function(fun=Denrell2007Solution, color='black', linetype='dashed', size=2/3)+
	#geom_contour(mapping = aes(learningRate, invTemperature, z = mean_proportionSafeChoice), breaks = c(0.5), colour = 'black')+
	scale_fill_gradient2(midpoint = 0.5, high = "red", mid = "grey90", low = "blue", breaks=c(0.1,0.5,0.9), labels=c(0.1,0.5,0.9) )+
	ylim(c(0,10))+
	#myTheme_gillsansMT()+
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle = 90))+
	theme(strip.text = element_text(size=12))+
	theme(legend.text = element_text(angle = 0))+
	theme(legend.position = 'top')+
	# theme(legend.position = NaN)+
	facet_grid(copyRate_factor ~ conformityExp_factor)+
	NULL -> Denrell2007_figure_social_learning
)

socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>%
	dplyr::filter(copyRate !=0) %>%
	# dplyr::filter(conformityExp %in% c(1, 4)) %>%
	ggplot() +
	geom_raster(mapping = aes(learningRate, invTemperature, fill = 1-mean_proportionSafeChoice), stat = 'identity') +
	labs(x=expression(paste('Learning rate ',alpha,sep="")),
		y=expression(paste('Inverse temperature ',beta,sep="")),
		#title='Gaussian noise\n mu=0.5; sigma=1',
		fill = "Proportion of \nchoosing \nthe risky option"
		#title='Gaussian noise\n mu=0.5; sigma=1', fill = "Proportion of\nsafe choice"
		)+
	#scale_fill_viridis(limits = c(0.45, 1), option="magma")+
	stat_function(fun=Denrell2007Solution, color='black', linetype='dashed', size=2/3)+
	#geom_contour(mapping = aes(learningRate, invTemperature, z = mean_proportionSafeChoice), breaks = c(0.5), colour = 'black')+
	scale_fill_gradient2(midpoint = 0.5, high = "red", mid = "grey90", low = "blue")+
	scale_y_continuous(limits = c(0, 10), breaks = c(1,5,10))+
	#myTheme_gillsansMT()+
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle = 90))+
	theme(strip.text.y = element_text(angle = 0))+
	theme(legend.text = element_text(angle = 90))+
	theme(legend.position = 'top')+
	# theme(legend.position = NaN)+
	facet_grid(copyRate_factor ~ conformityExp_factor)+
	NULL -> Denrell2007_figure_social_learning_full


## difference between social learners' performance and sole reinforcement learners
sigmaList = c(0, 0.25, 0.5, 0.75, 0.9)
thetaList = c(1, 2, 4, 8)
socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary_baseline <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>% dplyr::filter(copyRate==0)

socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>% arrange(copyRate, conformityExp) # re-ordering the data frame
socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$mean_proportionSafeChoice_diff_from_baseline <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$mean_proportionSafeChoice - rep(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary_baseline$mean_proportionSafeChoice, (length(sigmaList)-1) * length(thetaList) + 1)
# plot with Denrell (2007) curve
(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>%
	dplyr::filter(copyRate %in% c(0.25, 0.5)) %>%
	dplyr::filter(conformityExp %in% c(1, 4)) %>%
	ggplot() +
	geom_raster(mapping = aes(learningRate, invTemperature, fill = -mean_proportionSafeChoice_diff_from_baseline), stat = 'identity') +
	labs(x=expression(paste('Learning rate ',alpha,sep="")),
		y=expression(paste('Inverse temperature ',beta,sep="")),
		#title='Gaussian noise\n mu=0.5; sigma=1',
		fill = "Increases in \nthe risky choice by \nsocial learning"
		#title='Gaussian noise\n mu=0.5; sigma=1', fill = "Proportion of\nsafe choice"
		)+
	#scale_fill_viridis(limits = c(0.45, 1), option="magma")+
	stat_function(fun=Denrell2007Solution, color='black', linetype='dashed', size=2/3)+
	#geom_contour(mapping = aes(learningRate, invTemperature, z = mean_proportionSafeChoice), breaks = c(0.5), colour = 'black')+
	scale_fill_gradient2(midpoint = 0, high = "darkorange", mid = "grey90", low = "darkorchid3", breaks=c(-0.2,0,0.4),labels=c(-0.2,0,0.4))+
	ylim(c(0,10))+
	#myTheme_gillsansMT()+
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle = 90))+
	theme(strip.text = element_text(size=12))+
	theme(legend.text = element_text(angle = 0))+
	theme(legend.position = 'top')+
	facet_grid(copyRate_factor ~ conformityExp_factor)+
	NULL -> Denrell2007_figure_diff_from_baseline
)

socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>%
	dplyr::filter(copyRate != 0) %>%
	#dplyr::filter(conformityExp %in% c(1, 4)) %>%
	ggplot() +
	geom_raster(mapping = aes(learningRate, invTemperature, fill = -mean_proportionSafeChoice_diff_from_baseline), stat = 'identity') +
	labs(x=expression(paste('Learning rate ',alpha,sep="")),
		y=expression(paste('Inverse temperature ',beta,sep="")),
		#title='Gaussian noise\n mu=0.5; sigma=1',
		fill = "Increases in \nthe risky choice by \nsocial learning"
		#title='Gaussian noise\n mu=0.5; sigma=1', fill = "Proportion of\nsafe choice"
		)+
	#scale_fill_viridis(limits = c(0.45, 1), option="magma")+
	stat_function(fun=Denrell2007Solution, color='black', linetype='dashed', size=2/3)+
	#geom_contour(mapping = aes(learningRate, invTemperature, z = mean_proportionSafeChoice), breaks = c(0.5), colour = 'black')+
	scale_fill_gradient2(midpoint = 0, high = "darkorange", mid = "grey90", low = "darkorchid3")+
	scale_y_continuous(limits = c(0, 10), breaks = c(1,5,10))+
	#myTheme_gillsansMT()+
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle = 90))+
	# theme(axis.text.y = element_text(size=9))+
	theme(strip.text.y = element_text(angle = 0))+
	theme(legend.text = element_text(angle = 90))+
	theme(legend.position = 'top')+
	facet_grid(copyRate_factor ~ conformityExp_factor)+
	NULL -> Denrell2007_figure_diff_from_baseline_full

## Figure 1
figure1_left <- plot_grid(stable_2AB_scheme, Denrell2007_figure, labels = c('',''), ncol = 1, align = 'v')

(figure1 <- plot_grid(figure1_left, Denrell2007_figure_social_learning, Denrell2007_figure_diff_from_baseline, labels = c('','',''), ncol = 3, align = 'v')
)

## ========================================================
#
# figure 2: the collective rescue in the stable environment
#
# =========================================================
# Plot the individual learning from a different angle
# According to Denrell 2007, the asymptotic equilibrium of risky choice rate is
# Pr^* = 1 / (1 + exp(Z)), where Z = -beta/(2*(alpha-2)) * (alpha*(beta+1) - 2)
#
socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$hot_stove_suceptibility <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$learningRate * (socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary$invTemperature + 1)

Pr_when_beta = function (X, beta) {
	Z = -beta/(2*(X/(beta+1))-2) * (X - 2)
	return_vector <- 1 / (1 + exp(Z))
	return_vector[which(X >= beta+1)] <- NA
	return_vector[which(X == 0)] <- 1/2
	return(return_vector)
	# if(X <= beta) {
	# 	1 / (1 + exp(Z))
	# }else{
	# 	return(0)
	# }
}

ggplot(mapping = aes(x=X)) +
	mapply(
		function(b, co){ stat_function(data=data.frame(X=c(0,8)), fun=Pr_when_beta, args=list(beta=b), aes_q(color=co)) },
		seq(2,10,1),
		seq(2,10,1)
	)+
	geom_vline(xintercept=2, linetype='dashed')+
	geom_hline(yintercept=0.5, linetype='dashed')+
	myTheme_Helvetica()+
	scale_colour_viridis_c(expression(beta))+
	labs(
		y='Probability of choosing\nthe risky alternative',
		x=expression(alpha * (beta + 1))
	)+
	NULL

socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary_added <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>% rbind(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>% dplyr::filter(copyRate==0))

added_length <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary %>% dplyr::filter(copyRate==0) %>% nrow()
all_length <- socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary_added %>% nrow()
socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary_added$conformityExp[(all_length-added_length+1):all_length] <- 4
socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary_added$conformityExp_factor[(all_length-added_length+1):all_length] <- '\U03b8 = 4'

# Figure 2
(socialLearningParamSearch_riskPrem150_Gaussian_indivDiff_summary_added %>%
	dplyr::filter(copyRate %in% c(0, 0.25, 0.5)) %>%
	dplyr::filter(conformityExp %in% c(1, 4)) %>%
	dplyr::filter(invTemperature %in% c(3,5,7)) %>%
	ggplot() +
	geom_point(aes(hot_stove_suceptibility, 1-mean_proportionSafeChoice, colour=copyRate_factor)) +
	stat_function(data=data.frame(X=c(0,8),invTemperature_factor='β = 3'), fun=Pr_when_beta, n = 1001, args=list(beta=3)) +
	stat_function(data=data.frame(X=c(0,8),invTemperature_factor='β = 5'), fun=Pr_when_beta, n = 1001, args=list(beta=5)) +
	stat_function(data=data.frame(X=c(0,8),invTemperature_factor='β = 7'), fun=Pr_when_beta, n = 1001, args=list(beta=7)) +
	geom_vline(xintercept=2, linetype='dashed')+
	geom_hline(yintercept=0.5, linetype='dashed')+
	facet_grid(invTemperature_factor ~ conformityExp_factor) +
	scale_colour_viridis_d(expression(sigma))+
	labs(
		y='Probability of choosing\nthe risky alternative',
		x=expression(alpha * (beta + 1))
	)+
	xlim(c(0,8))+
	myTheme_Helvetica()+
	NULL -> collective_rescue_simplified_view
)


## ========================================================
#
# figure 2: Individual heterogeneity
#
# =========================================================
heterogeneous_groups_denrell2007task_alpha_summary <- read.csv('agentBasedSim/heterogeneous_groups_denrell2007task_alpha_summary.csv')
heterogeneous_groups_denrell2007task_beta_summary <- read.csv('agentBasedSim/heterogeneous_groups_denrell2007task_beta_summary.csv')
heterogeneous_groups_denrell2007task_sigma_summary <- read.csv('agentBasedSim/heterogeneous_groups_denrell2007task_sigma_summary.csv')
heterogeneous_groups_denrell2007task_theta_summary <- read.csv('agentBasedSim/heterogeneous_groups_denrell2007task_theta_summary.csv')

# alpha
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

# beta

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

# sigma
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

# theta
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

(
	heteroeneity_plot <- plot_grid(
	heterogeneous_groups_denrell2007task_alphaEffect_plot,
	heterogeneous_groups_denrell2007task_betaEffect_plot,
	heterogeneous_groups_denrell2007task_sigmaEffect_plot,
	heterogeneous_groups_denrell2007task_thetaEffect_plot,
	 labels = c('','','',''), ncol = 2, align = 'v')
)


## ========================================================
#
# figure 3: The dynamical model of collective behaviour
#
# =========================================================
## Analytical solution of asocial behavioural dynamics
noSocialCurveSimplest= function (n, e, pH, pL) {
    return ( -(n*(pH - pL)*((-1 + e)*pH + e*pL))/((pH + pL)*((-1 + e)*pH - e*pL)) )
}
zeroIsoclineSimplest = function (pH, pL) {
    return ( pH/(pH + pL) )
}

# Fig. 3a
# schematic figure
schematic_simplest <- cowplot::ggdraw() + cowplot::draw_image("dynamicsModel/schematic_simplest.001.Times.png", scale = 1)

### Fig. 3b
## function plot
diagonalLine = function(x){return(x)}
pLs <- c(0.1, 0.2, 0.4, 0.6)
(noSocialCurveSimplest_plot <- ggplot(data.frame(X=c(0,1)), aes(x=X)) +
    #stat_function(fun = diagonalLine, linetype = 'dashed', color = 'grey')+
    geom_vline(xintercept = 0.5, linetype='dashed', color='grey60')+
    geom_hline(yintercept = 0, linetype='dashed', color='grey60')+
    mapply(
        function(explorationrate) {
          stat_function(fun = noSocialCurveSimplest, args=list(n = 20, pH=0.7, pL=explorationrate), aes_q(color=explorationrate), size = 1)
        },
        pLs
    ) +
    annotate(geom="text", x=0.1, y=-13.5, label=expression(paste(italic(pl),' = 0.1',sep="")) ) +
    annotate(geom="text", x=0.1, y=-9.5, label=expression(paste(italic(pl),' = 0.2',sep=""))) +
    annotate(geom="text", x=0.1, y=-3.5, label=expression(paste(italic(pl),' = 0.4',sep=""))  ) +
    annotate(geom="text", x=0.1, y=0.5, label=expression(paste(italic(pl),' = 0.6',sep=""))) +
    scale_colour_viridis_c(option="cividis", direction=-1)+
    labs(color="pL",
       # y = expression(paste('Risk-seeking bias: ', N[R]^'*' - N[S]^'*',sep="")),
       y = expression(atop('Risk-seeking bias: ', paste(N[R]^'*' - N[S]^'*'))),
       x = expression(italic(e))
       # x = expression(paste('The rate of getting enchanted with R: ', italic(e),sep=""))
       # x = expression(atop('The rate of', paste('getting enchanted with R: ', italic(e), sep="")))
       )+
    myTheme_Helvetica()+
    theme(legend.position = 'none')+
    NULL)

# Fig. 3c

FigSocialLearningSimplest <- read_csv("dynamicsModel/FigSocialLearningSimplest.csv", col_names = FALSE)
names(FigSocialLearningSimplest) = c('e','c','f','pl','maxS','minS','diffS','maxR','minR','diffR','diffRS')
FigSocialLearningSimplest_sample = dplyr::sample_frac(tbl = FigSocialLearningSimplest, size = 0.00001)

FigSocialLearningSimplest$f_category = paste('f = ', FigSocialLearningSimplest$f, sep='')
FigSocialLearningSimplest$c_category = paste('c = ', FigSocialLearningSimplest$c, sep='')
FigSocialLearningSimplest$pl_category = paste('pl = ', FigSocialLearningSimplest$pl, sep='')
FigSocialLearningSimplest$f_category = factor(FigSocialLearningSimplest$f_category, levels = c('f = 1','f = 2','f = 10'))
(FigSocialLearningSimplest %>%
	dplyr::filter(c == 0) %>%
    ggplot(aes(x=pl))+
    #geom_line(mapping = aes(y = diffRS, group = pl, colour = pl)) +
    geom_raster(mapping = aes(y=e, fill = diffRS), stat = 'identity') +
    stat_function(fun=zeroIsoclineSimplest, args=list(pH=0.7), color='black', linetype='dashed', size=1)+
    labs(
        fill = expression(paste(N[R]^'*' - N[S]^'*')),
        # fill = 'Differences in\nrisk seeking and\nrisk aversion',
        y = expression(italic(e)),
        x = expression(italic(p[l]))
        # y = expression(atop('The rate of getting', paste('enchanted with R: ', italic(e), sep=""))),
        # x = expression(paste('The rate of exploration: ', italic(p[l]), sep=""))
        )+
    scale_fill_gradient2(midpoint = 0, low = "blue", mid = "grey90", high = "red", breaks=c(-10,0,10))+
    myTheme_Helvetica()+
    # theme(axis.text.x = element_text(angle = 90), legend.position='top')+
    theme(strip.text.y = element_text(angle = 0))+
    # facet_grid(c_category ~ f_category)+
    #geom_hline(yintercept = 0, linetype = 'dotted')+
    NULL -> FigSocialLearningSimplest_individual)

# Fig 3d
(FigSocialLearningSimplest %>%
	dplyr::filter((c == 0.25 | c == 0.5) & (f == 1 | f == 2)) %>%
    ggplot(aes(x=pl))+
    #geom_line(mapping = aes(y = diffRS, group = pl, colour = pl)) +
    geom_raster(mapping = aes(y=e, fill = diffRS), stat = 'identity') +
    stat_function(fun=zeroIsoclineSimplest, args=list(pH=0.7), color='black', linetype='dashed', size=1)+
    labs(
        fill = expression(paste(N[R]^'*' - N[S]^'*')),
        # fill = 'Differences in\nrisk seeking and\nrisk aversion',
        y = expression(italic(e)),
        x = expression(italic(p[l]))
        # y = expression(atop('The rate of', paste('getting enchanted with R: ', italic(e), sep=""))),
        # x = expression(paste('The rate of exploration: ', italic(p[l]), sep=""))
        )+
    scale_fill_gradient2(midpoint = 0, low = "blue", mid = "grey90", high = "red", breaks=c(-10,0,10))+
    myTheme_Helvetica()+
    theme(axis.text.x = element_text(angle = 90), legend.position='top')+
    theme(strip.text.y = element_text(angle = 270))+
    facet_grid(c_category ~ f_category)+
    #geom_hline(yintercept = 0, linetype = 'dotted')+
    NULL -> FigSocialLearningSimplest_social)

## Figure 3
(
	figure3_centre <- plot_grid(noSocialCurveSimplest_plot, FigSocialLearningSimplest_individual, labels = c('',''), ncol = 1, align = 'v')
)

# Figure 3 - additional (the bifurcation analysis)
sleqtableSimplest <- read_csv("dynamicsModel/sleqtableSimplest.csv", col_names = FALSE)
names(sleqtableSimplest) = c('f','S_0','c','e','R_eq')

sleqtableSimplest$RPreferingInitial = 0
sleqtableSimplest$direction = NA
sleqtableSimplest$R_0 = 20 - sleqtableSimplest$S_0
sleqtableSimplest$conformityExponent = paste(rep('θ = ', nrow(sleqtableSimplest)), sleqtableSimplest$f, sep ='')
sleqtableSimplest$conformityExponent = factor(sleqtableSimplest$conformityExponent, levels = c('θ = 0','θ = 1','θ = 2','θ = 10'))
sleqtableSimplest$e_factor = paste(rep('e = ', nrow(sleqtableSimplest)), sleqtableSimplest$e, sep ='')

for (i in 1:nrow(sleqtableSimplest)) {
    if (sleqtableSimplest$R_eq[i]>10) {
        sleqtableSimplest$RPreferingInitial[i] = 1
    }
    if (sleqtableSimplest$R_eq[i] - sleqtableSimplest$R_0[i] > 0) {
        sleqtableSimplest$direction[i] = 'upward'
    } else {
        sleqtableSimplest$direction[i] = 'downward'
    }
}

sleqtableSimplest <- sleqtableSimplest %>% dplyr::filter(c != 1)

(sleqtableSimplest %>%
    ggplot(aes(x=c))+
    geom_point(aes(y=R_0, colour=as.factor(RPreferingInitial), shape=direction), alpha=1/2)+
    geom_point(aes(y=R_eq))+
    labs(
        title = '',
        x=expression(paste('Social influence ', sigma, sep="")),
        y=expression(paste('Equilibrium density of ',N[R]^'*',sep=""))#,
        #title='Social influence\n (pH = 0.5; pL = 0.1; d=0.5; l=0.25)'
        )+
    scale_shape_manual(values=c('upward'=2, 'downward'=6), name='Stream\'s direction')+
    scale_color_manual(values=c('0'='#56B4E9','1'='#D55E00'), name='Risky choice regime')+
    myTheme_Helvetica()+
    theme(axis.text.x = element_text(angle = 90), legend.position='top')+
    facet_grid(e_factor ~ conformityExponent)+
    xlim(c(0,1))+
    #ylim(c(0,20))+
    geom_hline(yintercept=10, linetype='dashed')+
    theme(legend.position = 'none')+
    theme(strip.text.y = element_text(angle = 270))+
    NULL -> sleqtableSimplest_plot)



## ========================================================
#
# Figure for the experimental results
#
# =========================================================
library(ggpubr)
social_learning_model_validation_0820_data <- read.csv('experimentalAnalysis/social_learning_model_validation_0820_data.csv')
social_learning_model_validation_1022_riskID11_data <- read.csv('experimentalAnalysis/social_learning_model_validation_1022_riskID11_data.csv')
social_learning_model_validation_1022_riskID12_data <- read.csv('experimentalAnalysis/social_learning_model_validation_1022_riskID12_data.csv')

# ======================================
# 1-risky 1-safe (2-armed) task
# ======================================
fit_AL00_multiVar_LKJ_indiv_0820_globalparameters <- read.csv('experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_0820_globalparameters.csv')
fit_AL00_multiVar_LKJ_indiv_0820_parameters <- read.csv('experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_0820_parameters.csv')
fit_SL00_multiVar_LKJ_0820_globalparameters <- read.csv('experimentalAnalysis/fit_SL00_multiVar_LKJ_0820_globalparameters.csv')
fit_SL00_multiVar_LKJ_0820_parameters <- read.csv('experimentalAnalysis/fit_SL00_multiVar_LKJ_0820_parameters.csv')

behaviour_main_0820 <- read.csv("experimentalAnalysis/behaviour_main_0820.csv")
behaviour_indiv_0820 <-read.csv("experimentalAnalysis/behaviour_indiv_0820.csv")
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
	dplyr::filter(soc_mean > 3/10 & soc_mean < 6/10 & hot_stove_susceptibility_rounded < 6) %>%
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

# ======================================
# 1-risky 3-safe (4-armed) task
# ======================================

# fit result -- global parameters
fit_SL00_multiVar_LKJ_1022_globalparameters <- read.csv('experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_globalparameters.csv')
fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters <- read.csv('experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID11_indiv_riskID11Condition_globalparameters.csv')

## behavioural data summary
allBehaviour1022_group <- read.csv("experimentalAnalysis/allBehaviour1022_group.csv")
allBehaviour1022_group_riskID11 <- allBehaviour1022_group %>% dplyr::filter(riskDistributionId == 11)

fit_SL00_multiVar_LKJ_1022_parameters <- read.csv("experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_parameters.csv")

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


allBehaviour1022_indiv <- read.csv("experimentalAnalysis/allBehaviour1022_indiv.csv")
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
fit_AL00_multiVar_LKJ_indiv_riskID11_parameters <- read.csv('experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID11_parameters.csv')

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
	dplyr::filter(soc_mean > 3/10 & soc_mean < 6/10 & hot_stove_susceptibility_rounded < 6) %>%
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

# ======================================
# 2-risky 2-safe (4-armed) task
# ======================================

# fit result -- global parameters
fit_SL00_multiVar_LKJ_1022_globalparameters <- read.csv('experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_globalparameters.csv')
fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters <- read.csv('experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID12_indiv_riskID12Condition_globalparameters.csv')

## behavioural data summary
allBehaviour1022_group <- read.csv("experimentalAnalysis/allBehaviour1022_group.csv")
allBehaviour1022_group_riskID12 <- allBehaviour1022_group %>% dplyr::filter(riskDistributionId == 12)

fit_SL00_multiVar_LKJ_1022_parameters <- read.csv("experimentalAnalysis/fit_SL00_multiVar_LKJ_1022_parameters.csv")

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


allBehaviour1022_indiv <- read.csv("experimentalAnalysis/allBehaviour1022_indiv.csv")
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
fit_AL00_multiVar_LKJ_indiv_riskID12_parameters <- read.csv('experimentalAnalysis/fit_AL00_multiVar_LKJ_indiv_riskID12_parameters.csv')

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
	dplyr::filter(soc_mean > 3/10 & soc_mean < 6/10 & hot_stove_susceptibility_rounded < 6) %>%
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

ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.5,yend=0.5),colour="grey30", size=0.5) +
	geom_ribbon(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='grey20', alpha=1/2)+
	geom_ribbon(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==1&soc_mean_category=='mild'), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='orange', alpha=1/2)+
	geom_point(data = parameterfit_indiv_AL00_0820, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
	geom_point(data = fit_parameters_group_SL00_mcmc, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean, colour=soc_mean), shape = 20) +
	geom_line(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid), size=1.0)+
	geom_line(data=social_learning_model_validation_0820_summary%>%dplyr::filter(condition_dummy==1), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), size=1.0)+
	geom_line(data=social_learning_model_validation_0820_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), linetype = 'dashed', size=1.0)+
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Helvetica()+
	xlim(c(0,6.5))+
	labs(
		x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		y = 'Mean proportion of choosing\nthe optimal risky option',
		title = 'The 1-risky-1-safe task \n(N = 168)') +
	#theme(legend.position = c(0.85, 0.5))+
	theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12))+
	theme(legend.text = element_text(size=11))+
	NULL -> fig6_a

ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
	geom_ribbon(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='grey20', alpha=1/2)+
	geom_ribbon(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==1&soc_mean_category=='mild'), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='orange', alpha=1/2)+
	geom_point(data = fit_AL_indiv_riskID11_parameters, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
	geom_point(data = fit_SL00_riskID11_parameters, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean, colour=soc_mean_SL00_multiVar_LKJ), shape = 20) +
	geom_line(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid), size=1.0)+
	geom_line(data=social_learning_model_validation_1022_riskID11_summary%>%dplyr::filter(condition_dummy==1), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), size=1.0)+
	geom_line(data=social_learning_model_validation_1022_riskID11_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), linetype = 'dashed', size=1.0)+
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Helvetica()+
	xlim(c(0,6.5))+
	labs(
		x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		y = 'Mean proportion of choosing\nthe optimal risky option',
		title = 'The 1-risky-3-safe task \n(N = 148)') +
	theme(legend.position = NaN)+
	theme(legend.title = element_text(size=12))+
	theme(legend.text = element_text(size=11))+
	NULL -> fig6_b

ggplot() +
	geom_segment(aes(x=0,xend=6,y=0.25,yend=0.25),colour="grey30", size=0.5) +
  geom_ribbon(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='grey20', alpha=1/2)+
	geom_ribbon(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==1&soc_mean_category=='mild'), mapping=aes(hot_stove_susceptibility_rounded, ymin=proportionRiskyChoice_b2_lower, ymax=proportionRiskyChoice_b2_upper), fill='orange', alpha=1/2)+
	geom_point(data = fit_AL_indiv_riskID12_parameters, mapping=aes(hot_stove_susceptibility_trancated, risky_choice_mean), colour='grey20', shape = 17)+ # shape=5: diamond
  geom_point(data = fit_SL00_riskID12_parameters, mapping=aes(hot_stove_susceptibility_trancated,risky_choice_mean, colour=soc_mean_SL00_multiVar_LKJ), shape = 20) +
	geom_line(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==0), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid), size=1.0)+
	geom_line(data=social_learning_model_validation_1022_riskID12_summary%>%dplyr::filter(condition_dummy==1), mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), size=1.0)+
	geom_line(data=social_learning_model_validation_1022_riskID12_summary_reallyHighSigma, mapping=aes(hot_stove_susceptibility_rounded, proportionRiskyChoice_b2_mid, group=soc_mean_category, colour=mean(soc_mean)), linetype = 'dashed', size=1.0)+
	scale_colour_viridis_c(expression('Copying weight \U03C3'[i]), begin = 0.2, end = 0.9, option='plasma', direction=-1)+
	myTheme_Helvetica()+
	xlim(c(0,6.5))+
	labs(
		x = expression(atop('Susceptibility to the hot stove effect', paste(alpha[i], '(', beta[i], '+1)'))),
		y = 'Mean proportion of choosing\nthe optimal risky option',
		title = 'The 2-risky-2-safe task \n(N = 151)') +
	theme(legend.position = c(0.75, 0.7))+
	theme(legend.title = element_text(size=12))+
	theme(legend.text = element_text(size=11))+
	NULL -> fig6_c

(
  figure_exp_model_pred <- ggarrange(fig6_a, fig6_b, fig6_c
	# ,  common.legend = TRUE
	# ,  legend = 'right'
	, labels = c('','',''), ncol = 3, align = 'v'
	)
)
