## Plotting the results of parameter recovery test

library(tidyverse)
library(cowplot)

dropbox_path <- "~/Dropbox/wataru/papers/RiskySocialLearning/experiment/parameter_recovery/SL00/"
source(paste0(dropbox_path , 'functions.R'))

#######################################################
# recovery_param_set1
#######################################################
true_global_parameters_4ab_riskID11 <- read.csv(paste0(dropbox_path , "recovery_param_set1/true_global_parameters_4ab_riskID11.csv"))
true_global_parameters_4ab_riskID12 <- read.csv(paste0(dropbox_path , "recovery_param_set1/true_global_parameters_4ab_riskID12.csv"))
true_global_parameters_set1 <- true_global_parameters_4ab_riskID11 %>% rbind(true_global_parameters_4ab_riskID12)
true_global_parameters_set1$task <- rep(1:2, each= nrow(true_global_parameters_4ab_riskID11))

true_individual_parameters_4ab_riskID11 <- read.csv(paste0(dropbox_path , "recovery_param_set1/true_individual_parameters_4ab_riskID11.csv"))
true_individual_parameters_4ab_riskID12 <- read.csv(paste0(dropbox_path , "recovery_param_set1/true_individual_parameters_4ab_riskID12.csv"))
true_individual_set1 <- rbind(true_individual_parameters_4ab_riskID11, true_individual_parameters_4ab_riskID12)
true_individual_set1$task <- rep(1:2, each= nrow(true_individual_parameters_4ab_riskID11))
true_individual_set1$indiv <- 1:nrow(true_individual_set1)

glogal_set1 <- read.csv(paste0(dropbox_path , "recovery_param_set1/fit_SL00_multiVar_LKJ_param_recov_4ab_globalparameters.csv"))
individual_set1 <- read.csv(paste0(dropbox_path , "recovery_param_set1/fit_SL00_multiVar_LKJ_param_recov_4ab_parameters.csv"))


## =======================================================================
## Global parameters
#  =======================================================================
true_global_parameters_set1_long <- true_global_parameters_set1 %>% pivot_longer(-task, names_to = "variable", values_to = "true")
true_global_parameters_set1_long$task_factor <- true_global_parameters_set1_long$task
true_global_parameters_set1_long$task_factor[which(true_global_parameters_set1_long$task_factor==1)] <- '1-Risky 3-Safe Task'
true_global_parameters_set1_long$task_factor[which(true_global_parameters_set1_long$task_factor==2)] <- '2-Risky 2-Safe Task'

fit_global_parameters_set1_long <-
  glogal_set1[1:(2 * 4 * 2),] %>%
  separate(variable, into=c('variable', 'task]'), sep='\\[') %>%
  separate('task]', into=c('task', ']'), sep='\\]')

fit_global_parameters_set1_long$task_factor <- fit_global_parameters_set1_long$task
fit_global_parameters_set1_long$task_factor[which(fit_global_parameters_set1_long$task_factor==1)] <- '1-Risky 3-Safe Task'
fit_global_parameters_set1_long$task_factor[which(fit_global_parameters_set1_long$task_factor==2)] <- '2-Risky 2-Safe Task'

ggplot(data = fit_global_parameters_set1_long) +
	geom_point(data=true_global_parameters_set1_long, aes(variable, true), size=2, colour='red')+
	geom_errorbar(aes(variable, ymin=q5, ymax=q95)) +
	geom_point(aes(variable, mean)) +
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle=90))+
	labs(title='Global parameters' , x = '' , y = 'Parameter recovery 1')+
	facet_grid(task_factor ~ .) +
	NULL -> global_recovery_set1


## =======================================================================
## Individual parameters (showing correlations)
#  =======================================================================

# correlations between indiv parameters
individual_parameters_set1 <- true_individual_set1 %>% left_join(individual_set1, by = c("indiv" = "sub"))

## Correlation coefficients
# cor.test(individual_parameters_set1$alpha, individual_parameters_set1$alpha_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_set1$log_beta, individual_parameters_set1$log_beta_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_set1$soc0, individual_parameters_set1$logit_soc0_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_set1$theta, individual_parameters_set1$theta_mean_SL00_multiVar_LKJ)

cors_alpha = individual_parameters_set1 %>%
  summarize(cor = round(cor.test(alpha, alpha_mean_SL00_multiVar_LKJ)$estimate,3)
  	, p = round(cor.test(alpha, alpha_mean_SL00_multiVar_LKJ)$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 0.1, fit = 0.75)

cors_log_beta = individual_parameters_set1 %>%
  summarize(cor = round(cor.test(log_beta, log_beta_mean_SL00_multiVar_LKJ)$estimate,3)
  	, p = round(cor.test(log_beta, log_beta_mean_SL00_multiVar_LKJ)$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 2, fit = 0)

cors_soc0 = individual_parameters_set1 %>%
  summarize(cor = round(cor.test(convert_alphaRaw_to_alpha(soc0), convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ))$estimate,3)
  	, p = round(cor.test(convert_alphaRaw_to_alpha(soc0), convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ))$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 0.1, fit = 0.75)

cors_theta = individual_parameters_set1 %>%
  summarize(cor = round(cor.test(theta, theta_mean_SL00_multiVar_LKJ)$estimate,3)
  	, p = round(cor.test(theta, theta_mean_SL00_multiVar_LKJ)$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 2, fit = 0)

alpha_plot <- ggplot(individual_parameters_set1) +
	geom_point(aes(alpha, alpha_mean_SL00_multiVar_LKJ, colour=sqrt((alpha - alpha_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('Learning rate ', alpha[i])) , x = 'True' , y = 'Fit' ) +
	geom_text(data=cors_alpha, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	xlim(c(0,1)) +
	ylim(c(0,1)) +
	myTheme_Helvetica() +
	NULL

log_beta_plot <- ggplot(individual_parameters_set1) +
	geom_point(aes(log_beta, log_beta_mean_SL00_multiVar_LKJ, colour=sqrt((log_beta - log_beta_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('log inv. temp. ', beta[i])) , x = 'True' , y = 'Fit' ) +
	geom_text(data=cors_log_beta, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	myTheme_Helvetica() +
	NULL

soc0_plot <- ggplot(individual_parameters_set1) +
	geom_point(aes(convert_alphaRaw_to_alpha(soc0), convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ), colour=sqrt((convert_alphaRaw_to_alpha(soc0) - convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ))^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('Copying weight ', sigma[i])) , x = 'True' , y = 'Fit' ) +
	xlim(c(0,1)) +
	ylim(c(0,1)) +
	geom_text(data=cors_soc0, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	myTheme_Helvetica() +
	NULL

theta_plot <- ggplot(individual_parameters_set1) +
	geom_point(aes(theta, theta_mean_SL00_multiVar_LKJ, colour=sqrt((theta - theta_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('Conformity exponent ', theta[i])) , x = 'True' , y = 'Fit' ) +
	geom_text(data=cors_theta, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	myTheme_Helvetica() +
	NULL


# merging figures
indiv_parameters_true_vs_fit_set1_l <- plot_grid(
	alpha_plot
	# , log_beta_plot
	, soc0_plot
	# , theta_plot
	, labels = c('','','',''), ncol = 1, align = 'v')

indiv_parameters_true_vs_fit_set1_r <- plot_grid(
	# alpha_plot
	log_beta_plot
	# , soc0_plot
	, theta_plot
	, labels = c('','','',''), ncol = 1, align = 'v')

param_recov_res_set1 <- plot_grid(global_recovery_set1, indiv_parameters_true_vs_fit_set1_l, indiv_parameters_true_vs_fit_set1_r
	, labels = c('','')
	# , ncol=1, align = 'v', axis = 'lr'
	, ncol=3, align = 'v', axis = 'bt'
	)

ggsave(file = paste0(dropbox_path , 'param_recov_res_set1.pdf')
	, plot = param_recov_res_set1
	, width=12, height=6)






#######################################################
# recovery_param_set3
#######################################################
true_global_parameters_4ab_riskID11 <- read.csv(paste0(dropbox_path , "recovery_param_set3/true_global_parameters_4ab_riskID11.csv"))
true_global_parameters_4ab_riskID12 <- read.csv(paste0(dropbox_path , "recovery_param_set3/true_global_parameters_4ab_riskID12.csv"))
true_global_parameters_set3 <- true_global_parameters_4ab_riskID11 %>% rbind(true_global_parameters_4ab_riskID12)
true_global_parameters_set3$task <- rep(1:2, each= nrow(true_global_parameters_4ab_riskID11))

true_individual_parameters_4ab_riskID11 <- read.csv(paste0(dropbox_path , "recovery_param_set3/true_individual_parameters_4ab_riskID11.csv"))
true_individual_parameters_4ab_riskID12 <- read.csv(paste0(dropbox_path , "recovery_param_set3/true_individual_parameters_4ab_riskID12.csv"))
true_individual_set3 <- rbind(true_individual_parameters_4ab_riskID11, true_individual_parameters_4ab_riskID12)
true_individual_set3$task <- rep(1:2, each= nrow(true_individual_parameters_4ab_riskID11))
true_individual_set3$indiv <- 1:nrow(true_individual_set3)

glogal_set3 <- read.csv(paste0(dropbox_path , "recovery_param_set3/fit_SL00_multiVar_LKJ_param_recov_4ab_globalparameters.csv"))
individual_set3 <- read.csv(paste0(dropbox_path , "recovery_param_set3/fit_SL00_multiVar_LKJ_param_recov_4ab_parameters.csv"))


## =======================================================================
## Global parameters
#  =======================================================================
true_global_parameters_set3_long <- true_global_parameters_set3 %>% pivot_longer(-task, names_to = "variable", values_to = "true")
true_global_parameters_set3_long$task_factor <- true_global_parameters_set3_long$task
true_global_parameters_set3_long$task_factor[which(true_global_parameters_set3_long$task_factor==1)] <- '1-Risky 3-Safe Task'
true_global_parameters_set3_long$task_factor[which(true_global_parameters_set3_long$task_factor==2)] <- '2-Risky 2-Safe Task'

fit_global_parameters_set3_long <-
  glogal_set3[1:(2 * 4 * 2),] %>%
  separate(variable, into=c('variable', 'task]'), sep='\\[') %>%
  separate('task]', into=c('task', ']'), sep='\\]')

fit_global_parameters_set3_long$task_factor <- fit_global_parameters_set3_long$task
fit_global_parameters_set3_long$task_factor[which(fit_global_parameters_set3_long$task_factor==1)] <- '1-Risky 3-Safe Task'
fit_global_parameters_set3_long$task_factor[which(fit_global_parameters_set3_long$task_factor==2)] <- '2-Risky 2-Safe Task'

ggplot(data = fit_global_parameters_set3_long) +
	geom_point(data=true_global_parameters_set3_long, aes(variable, true), size=2, colour='red')+
	geom_errorbar(aes(variable, ymin=q5, ymax=q95)) +
	geom_point(aes(variable, mean)) +
	myTheme_Helvetica()+
	theme(axis.text.x = element_text(angle=90))+
	labs(title='Global parameters' , x = '' , y = 'Parameter recovery 2')+
	facet_grid(task_factor ~ .) +
	NULL -> global_recovery_set3


## =======================================================================
## Individual parameters (showing correlations)
#  =======================================================================

# correlations between indiv parameters
individual_parameters_set3 <- true_individual_set3 %>% left_join(individual_set3, by = c("indiv" = "sub"))

## Correlation coefficients
# cor.test(individual_parameters_set3$alpha, individual_parameters_set3$alpha_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_set3$log_beta, individual_parameters_set3$log_beta_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_set3$soc0, individual_parameters_set3$logit_soc0_mean_SL00_multiVar_LKJ)
# cor.test(individual_parameters_set3$theta, individual_parameters_set3$theta_mean_SL00_multiVar_LKJ)

cors_alpha = individual_parameters_set3 %>%
  summarize(cor = round(cor.test(alpha, alpha_mean_SL00_multiVar_LKJ)$estimate,3)
  	, p = round(cor.test(alpha, alpha_mean_SL00_multiVar_LKJ)$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 0.1, fit = 0.75)

cors_log_beta = individual_parameters_set3 %>%
  summarize(cor = round(cor.test(log_beta, log_beta_mean_SL00_multiVar_LKJ)$estimate,3)
  	, p = round(cor.test(log_beta, log_beta_mean_SL00_multiVar_LKJ)$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 2, fit = 0)

cors_soc0 = individual_parameters_set3 %>%
  summarize(cor = round(cor.test(convert_alphaRaw_to_alpha(soc0), convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ))$estimate,3)
  	, p = round(cor.test(convert_alphaRaw_to_alpha(soc0), convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ))$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 0.1, fit = 0.75)

cors_theta = individual_parameters_set3 %>%
  summarize(cor = round(cor.test(theta, theta_mean_SL00_multiVar_LKJ)$estimate,3)
  	, p = round(cor.test(theta, theta_mean_SL00_multiVar_LKJ)$p.value,3)
  	) %>%
  as.data.frame() %>%
  mutate(p = ifelse(p < 0.01, "p < 0.01", paste("p =", p )), true = 2, fit = 0)

alpha_plot <- ggplot(individual_parameters_set3) +
	geom_point(aes(alpha, alpha_mean_SL00_multiVar_LKJ, colour=sqrt((alpha - alpha_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('Learning rate ', alpha[i])) , x = 'True' , y = 'Fit' ) +
	geom_text(data=cors_alpha, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	xlim(c(0,1)) +
	ylim(c(0,1)) +
	myTheme_Helvetica() +
	NULL

log_beta_plot <- ggplot(individual_parameters_set3) +
	geom_point(aes(log_beta, log_beta_mean_SL00_multiVar_LKJ, colour=sqrt((log_beta - log_beta_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('log inv. temp. ', beta[i])) , x = 'True' , y = 'Fit' ) +
	geom_text(data=cors_log_beta, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	myTheme_Helvetica() +
	NULL

soc0_plot <- ggplot(individual_parameters_set3) +
	geom_point(aes(convert_alphaRaw_to_alpha(soc0), convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ), colour=sqrt((convert_alphaRaw_to_alpha(soc0) - convert_alphaRaw_to_alpha(logit_soc0_mean_SL00_multiVar_LKJ))^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('Copying weight ', sigma[i])) , x = 'True' , y = 'Fit' ) +
	ylim(c(0,1)) +
	ylim(c(0,1)) +
	geom_text(data=cors_soc0, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	myTheme_Helvetica() +
	NULL

theta_plot <- ggplot(individual_parameters_set3) +
	geom_point(aes(theta, theta_mean_SL00_multiVar_LKJ, colour=sqrt((theta - theta_mean_SL00_multiVar_LKJ)^2)))+
	scale_colour_viridis_c('Difference') +
	labs(title = expression(paste('Conformity exponent ', theta[i])) , x = 'True' , y = 'Fit' ) +
	geom_text(data=cors_theta, mapping = aes(x = true, y = fit, label = paste("r =",cor)), hjust = 0, vjust = 0) +
	myTheme_Helvetica() +
	NULL


# merging figures
indiv_parameters_true_vs_fit_set3_l <- plot_grid(
	alpha_plot
	# , log_beta_plot
	, soc0_plot
	# , theta_plot
	, labels = c('','','',''), ncol = 1, align = 'v')

indiv_parameters_true_vs_fit_set3_r <- plot_grid(
	# alpha_plot
	log_beta_plot
	# , soc0_plot
	, theta_plot
	, labels = c('','','',''), ncol = 1, align = 'v')

param_recov_res_set3 <- plot_grid(global_recovery_set3, indiv_parameters_true_vs_fit_set3_l, indiv_parameters_true_vs_fit_set3_r
	, labels = c('','')
	# , ncol=1, align = 'v', axis = 'lr'
	, ncol=3, align = 'v', axis = 'bt'
	)

ggsave(file = paste0(dropbox_path , 'param_recov_res_set3.pdf')
	, plot = param_recov_res_set3
	, width=12, height=6)



###################################
## Merging set 1 and 2
##################################

param_recov_res_both <- plot_grid(
	param_recov_res_set1
	, param_recov_res_set3
	, labels = c('','')
	, ncol=1, align = 'v', axis = 'lr'
	)

ggsave(file = paste0(dropbox_path , 'param_recov_res_both.pdf')
	, plot = param_recov_res_both
	, width=12, height=12)
