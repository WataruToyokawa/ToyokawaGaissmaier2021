# ToyokawaGaissmaier2021
This is a repository on the paper entitled "Conformist social learning leads to self-organised prevention against adverse bias in risky decision making" written by Wataru Toyokawa co-authored with Wolfgang Gaissmaier. The figures shown in the main text can be reproduced by running `drawing_figures.R`. 

## Individual-based simulation (R)
* R code found in the folder agentBasedSim will generate data used for Figure 1, 2, and 3. 

## Differential equation model
* Model is run by Mathematica
* A numerical simulation generates data for Figures 4 & 5

## Analaysis code for the experiment (Fig. 6 and Table 1)
* Code for the data analyses (i.e. model fitting) are available in `experimental_analyses.Rmd`. 
* The post-hoc siumulation is done by `posthoc_simulation.R`
* For the model analysis, we used `cmdstanr` package
