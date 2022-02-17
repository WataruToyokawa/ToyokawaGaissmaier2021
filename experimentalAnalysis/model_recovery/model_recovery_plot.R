##############################################################################################
##
## Model fitting and Bayesian model comparison for simulated data. 
## I considered the following models
## (1) Asocial RL
## (2) Decision biasing social influence 
## (3) Value shaping social influence
##
## Author: 豊川 航 / TOYOKAWA, Wataru (11 January 2022, University of Konstanz, Germany)
###############################################################################################

if(FALSE) rm(list=ls(all=TRUE)) # cleaning the workspace 
library(tidyverse)
library(LaplacesDemon)
library(cowplot)
#rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores()) # regacy from rstan?

## Load Functions
source('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/data_4ab_1022/analysis/functions.R')
source('~/Dropbox/wataru/papers/RiskySocialLearning/experiment/parameter_recovery/SL00/functions.R')

# path -- drop box
dropbox_path <- "~/Dropbox/wataru/papers/RiskySocialLearning/experiment/parameter_recovery/SL00/recovery_param_expfit/model_recovery/"

###############################################################
# ---- Bayesian Model Selection -----
# Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ. 
#       Bayesian model selection for group studies. NeuroImage. 2009;46(4):1004–17.
###############################################################

# ---- AL was the ground truth --------
# Loading data (if you start directly from here)
fit_SL_to_AL_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_SL_to_AL_pseudodata_parameters.csv"))
fit_VS_to_AL_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_VS_to_AL_pseudodata_parameters.csv"))
fit_AL_to_AL_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_AL_to_AL_pseudodata_parameters.csv"))

## STEP 1: Using WAIC as an approximation of model log-evidence
#-----------------------------------
WAICtable_AL_is_true <- data.frame(
	subject=as.character(fit_SL_to_AL_pseudodata_parameters$sub),
	WAIC_AL00=fit_AL_to_AL_pseudodata_parameters$waic_AL00_multiVar_LKJ,
	WAIC_VS00=fit_VS_to_AL_pseudodata_parameters$waic_VS00_multiVar_LKJ,
	WAIC_SL00=fit_SL_to_AL_pseudodata_parameters$waic_SL00_multiVar_LKJ
	)
# STEP 2. Running the algorithm which updates estimates of the Dirichlet density parameter, α
#         iteretively until convergence
#------------------------------------------------------------------------------------------
num_model = 3
D_alpha0 = rep(1, num_model)
D_alpha = D_alpha0
for(i in 1:500){
	u_WAIC_1 = exp(-WAICtable_AL_is_true$WAIC_AL00 + digamma(D_alpha[1]) - digamma(sum(D_alpha)))
	u_WAIC_2 = exp(-WAICtable_AL_is_true$WAIC_VS00 + digamma(D_alpha[2]) - digamma(sum(D_alpha)))
	u_WAIC_3 = exp(-WAICtable_AL_is_true$WAIC_SL00 + digamma(D_alpha[3]) - digamma(sum(D_alpha)))
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
rWAIC_AL_is_true_for_fig = data.frame(
	models = c('Asocial','Value shaping','Decision biasing'),
	r_WAIC = r_WAIC)
rWAIC_AL_is_true_for_fig$models = factor(c('Asocial','Value shaping','Decision biasing'),levels=c('Asocial','Value shaping','Decision biasing'))

write.csv(rWAIC_AL_is_true_for_fig, paste0(dropbox_path , "rWAIC_AL_is_true_for_fig.csv"), row.names=FALSE)

# STEP 4: Plot (Exceedance Probability)
# (1) sample: (r1, r2, r3) ~ Dirichlet(alpha_1, alpha_2, alpha_3); 
# (2) Count cases where r3 > r1 and r3 > r2
# -------------------------------------------------------------------------------------------
total_sample_N = 10000
dir_sample = rdirichlet(total_sample_N, D_alpha)
xp_AL = (which(apply(dir_sample, MARGIN=1, FUN=which.max) == 1) %>% length()) / total_sample_N
xp_AL


# ---- VS was the ground truth --------
# Loading data (if you start directly from here)
fit_SL_to_VS_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_SL_to_VS_pseudodata_parameters.csv"))
fit_VS_to_VS_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_VS_to_VS_pseudodata_parameters.csv"))
fit_AL_to_VS_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_AL_to_VS_pseudodata_parameters.csv"))

## STEP 1: Using WAIC as an approximation of model log-evidence
#-----------------------------------
WAICtable_VS_is_true <- data.frame(
	subject=as.character(fit_SL_to_VS_pseudodata_parameters$sub),
	WAIC_AL00=fit_AL_to_VS_pseudodata_parameters$waic_AL00_multiVar_LKJ,
	WAIC_VS00=fit_VS_to_VS_pseudodata_parameters$waic_VS00_multiVar_LKJ,
	WAIC_SL00=fit_SL_to_VS_pseudodata_parameters$waic_SL00_multiVar_LKJ
	)
# STEP 2. Running the algorithm which updates estimates of the Dirichlet density parameter, α
#         iteretively until convergence
#------------------------------------------------------------------------------------------
num_model = 3
D_alpha0 = rep(1, num_model)
D_alpha = D_alpha0
for(i in 1:500){
	u_WAIC_1 = exp(-WAICtable_VS_is_true$WAIC_AL00 + digamma(D_alpha[1]) - digamma(sum(D_alpha)))
	u_WAIC_2 = exp(-WAICtable_VS_is_true$WAIC_VS00 + digamma(D_alpha[2]) - digamma(sum(D_alpha)))
	u_WAIC_3 = exp(-WAICtable_VS_is_true$WAIC_SL00 + digamma(D_alpha[3]) - digamma(sum(D_alpha)))
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
rWAIC_VS_is_true_for_fig = data.frame(
	models = c('Asocial','Value shaping','Decision biasing'),
	r_WAIC = r_WAIC)
rWAIC_VS_is_true_for_fig$models = factor(c('Asocial','Value shaping','Decision biasing'),levels=c('Asocial','Value shaping','Decision biasing'))

# STEP 4: Plot (Exceedance Probability)
# (1) sample: (r1, r2, r3) ~ Dirichlet(alpha_1, alpha_2, alpha_3); 
# (2) Count cases where r3 > r1 and r3 > r2
# -------------------------------------------------------------------------------------------
total_sample_N = 10000
dir_sample = rdirichlet(total_sample_N, D_alpha)
xp_VS = (which(apply(dir_sample, MARGIN=1, FUN=which.max) == 2) %>% length()) / total_sample_N
xp_VS



# ---- SL was the ground truth --------
# Loading data (if you start directly from here)
fit_SL_to_SL_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_SL_to_SL_pseudodata_parameters.csv"))
fit_VS_to_SL_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_VS_to_SL_pseudodata_parameters.csv"))
fit_AL_to_SL_pseudodata_parameters = read_csv(paste0(dropbox_path , "fit_AL_to_SL_pseudodata_parameters.csv"))

## STEP 1: Using WAIC as an approximation of model log-evidence
#-----------------------------------
WAICtable_SL_is_true <- data.frame(
	subject=as.character(fit_SL_to_SL_pseudodata_parameters$sub),
	WAIC_AL00=fit_AL_to_SL_pseudodata_parameters$waic_AL00_multiVar_LKJ,
	WAIC_VS00=fit_VS_to_SL_pseudodata_parameters$waic_VS00_multiVar_LKJ,
	WAIC_SL00=fit_SL_to_SL_pseudodata_parameters$waic_SL00_multiVar_LKJ
	)
# STEP 2. Running the algorithm which updates estimates of the Dirichlet density parameter, α
#         iteretively until convergence
#------------------------------------------------------------------------------------------
num_model = 3
D_alpha0 = rep(1, num_model)
D_alpha = D_alpha0
for(i in 1:500){
	u_WAIC_1 = exp(-WAICtable_SL_is_true$WAIC_AL00 + digamma(D_alpha[1]) - digamma(sum(D_alpha)))
	u_WAIC_2 = exp(-WAICtable_SL_is_true$WAIC_VS00 + digamma(D_alpha[2]) - digamma(sum(D_alpha)))
	u_WAIC_3 = exp(-WAICtable_SL_is_true$WAIC_SL00 + digamma(D_alpha[3]) - digamma(sum(D_alpha)))
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
rWAIC_SL_is_true_for_fig = data.frame(
	models = c('Asocial','Value shaping','Decision biasing'),
	r_WAIC = r_WAIC)
rWAIC_SL_is_true_for_fig$models = factor(c('Asocial','Value shaping','Decision biasing'),levels=c('Asocial','Value shaping','Decision biasing'))

write.csv(rWAIC_SL_is_true_for_fig, paste0(dropbox_path , "rWAIC_SL_is_true_for_fig.csv"), row.names=FALSE)

# STEP 4: Plot (Exceedance Probability)
# (1) sample: (r1, r2, r3) ~ Dirichlet(alpha_1, alpha_2, alpha_3); 
# (2) Count cases where r3 > r1 and r3 > r2
# -------------------------------------------------------------------------------------------
total_sample_N = 10000
dir_sample = rdirichlet(total_sample_N, D_alpha)
xp_SL = (which(apply(dir_sample, MARGIN=1, FUN=which.max) == 3) %>% length()) / total_sample_N
xp_SL





# --- Plot ----
rWAIC_merged = rbind(rWAIC_AL_is_true_for_fig, rWAIC_VS_is_true_for_fig, rWAIC_SL_is_true_for_fig)
rWAIC_merged$truth = factor(rep(c('Asocial','Value shaping','Decision biasing'), each = 3), levels=c('Asocial','Value shaping','Decision biasing'))
rWAIC_merged$models_ab = factor(rep(c('AL','VS','DB'), 3), levels=c('AL','VS','DB'))

(rWAIC_merged %>% ggplot(aes(models_ab, truth)) + 
    geom_raster(aes(fill = r_WAIC)) +
    scale_fill_gradient2(low = "white", high = "black", name = 'Model frequency')+
    annotate(geom="text", x=1, y=1, label=paste('XP = \n', round(10*xp_AL)/10), size=6, colour='white') +
    annotate(geom="text", x=2, y=2, label=paste('XP = \n', xp_VS), size=6, colour='white') +
    annotate(geom="text", x=3, y=3, label=paste('XP = \n', xp_SL), size=6, colour='white') +
    labs(x='Recovered by fit', y='Simulated model', title='Model recovery') +
    myTheme_Arial()+
    theme(
		plot.title = element_text(size=16, family='Arial', face='bold'),
	)+
    NULL -> model_recovery_plot
)
ggsave(file = paste0(dropbox_path, "model_recovery_plot.png"), plot = model_recovery_plot, width=9, height=6)




# --- Merging model recovery result with model selection result ---
# Note: First run the file "Bayesian_model_selection.R" in the eLife/Revision2/model_fitting/fitting_results/ folder

# Merging all three plots
(model_comparison_result =  plot_grid(
    model_recovery_plot
    , r_plot_0820
    , r_plot_task12
    , r_plot_revision_exp
    , axis = 'bt'
    , labels = c('a','b','c','d'), ncol = 2, align = 'v'
    )
)

ggsave(file = "~/Dropbox/wataru/papers/RiskySocialLearning/draft/submissions/eLife/Revision2/model_fitting/fitting_results/model_comparison_result.png"
, plot = model_comparison_result, width=15, height=9)
