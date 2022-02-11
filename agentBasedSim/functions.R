###############################################################################################
##
## Searching for conditions under which teaching strategy change collective behaviour
## Wataru Toyokawa
## 18 Jan. 2020
##
##
###############################################################################################

## Functions
payoffGenerateGaussian = function(groupSize, m, mu_sure, mu_risky, sigma_sure, sigma_risky) {
    payoff = numeric(groupSize)
    payoff[which(m==1)] = rnorm(length(which(m==1)), mean=mu_sure, sd=sigma_sure)
    payoff[which(m==2)] = rnorm(length(which(m==2)), mean=mu_risky, sd=sigma_risky)
    return (payoff)
}
payoffGenerate = function(groupSize, m, r1, r2, payoff_H1, payoff_L1, payoff_H2, payoff_L2) {
	payoff = numeric(groupSize)
	payoff[which(m==1|m==2)] = sample(c(payoff_L1,payoff_H1), length(which(m==1|m==2)), prob=c(1-r1, r1), replace=TRUE)
	payoff[which(m==3|m==4)] = sample(c(payoff_L2,payoff_H2), length(which(m==3|m==4)), prob=c(1-r2, r2), replace=TRUE)
	# Add a small noise to the payoff
	payoff = payoff + rnorm(groupSize,0,0.2) #-0.66 ~ +0.66
	return (payoff)
}
payoffGenerate_twoOptionVer = function(groupSize, m, r1, r2, payoff_L1, payoff_H1, payoff_L2, payoff_H2) {
	payoff = numeric(groupSize)
	payoff[which(m==1)] = sample(c(payoff_L1,payoff_H1), length(which(m==1)), prob=c(1-r1, r1), replace=TRUE)
	payoff[which(m==2)] = sample(c(payoff_L2,payoff_H2), length(which(m==2)), prob=c(1-r2, r2), replace=TRUE)
	# Add a small noise to the payoff
	payoff = payoff + rnorm(groupSize,0,0.02) #-0.066 ~ +0.066
	return (payoff)
}
payoffGenerate_org2 = function(groupSize, m, r1, r2, payoff_H1, payoff_L1, payoff_H2, payoff_L2) {
    payoff = numeric(groupSize)
    payoff[which(m==1)] = sample(c(payoff_L1,payoff_H1), length(which(m==1)), prob=c(1-r1, r1), replace=TRUE)
    payoff[which(m==2)] = sample(c(payoff_L2,payoff_H2), length(which(m==2)), prob=c(1-r2, r2), replace=TRUE)
    # Add a small noise to the payoff
    payoff = payoff + rnorm(groupSize,0,0.2) #-0.6~0.6
    return (payoff)
}
payoffGenerate_org = function(groupSize, m, r1, r2, k1, k2) {
    payoff = numeric(groupSize)
    payoff[which(m==1)] = sample(c(0,k1), length(which(m==1)), prob=c(1-r1, r1), replace=TRUE)
    payoff[which(m==2)] = sample(c(0,k2), length(which(m==2)), prob=c(1-r2, r2), replace=TRUE)
    return (payoff)
}
divideVector = function(numerator, denominator) {
    return( numerator/denominator )
}
multiplyBeta = function(Q, beta) {
    return( beta*Q )
}
expCeiling = function(x) {
    if(length(which(is.infinite(exp(x))))==0) {
        return( exp(x) )
    }else{
        result = exp(x)
        result[which(is.infinite(exp(x)))] = 8.988466e+307
        return( result ) # 2^1023
    }
}
powerTheta = function(f, theta) {
    return( f^theta )
}
frequencyDependentCopy_org = function(F, lastChoice, theta, numOptions) {
    totalN=length(lastChoice)
    f = rep(F, totalN)
    f[lastChoice+((1:totalN)-1)*numOptions] = f[lastChoice+((1:totalN)-1)*numOptions] - 1
    f_matrix = matrix(f, ncol=totalN)
    ftheta = apply(f_matrix, 1, powerTheta, theta=theta) %>% t()
    denom = apply(ftheta,2,sum)
    if(denom == 0) denom = 1
    return ( apply(ftheta, 1, divideVector, denominator=denom) %>% t() )
}
frequencyDependentCopy = function(F, lastChoice, theta, numOptions) {
    totalN=length(lastChoice)
    f = rep(F, totalN)
    f[lastChoice+((1:totalN)-1)*numOptions] = f[lastChoice+((1:totalN)-1)*numOptions] - 1
    f_matrix = matrix(f, ncol=totalN)
    ftheta = f_matrix ^ rep(theta, each= numOptions)
    denom = apply(ftheta,2,sum)
    if(length(which(denom == 0)>0)) denom[which(denom==0)] = 1
    return ( apply(ftheta, 1, divideVector, denominator=denom) %>% t() )
}
myTheme_legend = function() {
  theme(
    #legend.position = 'right',
    legend.title = element_text(size=16, family="Times" ,colour='black'),
    legend.text = element_text(size=15, family="Times" ,colour='black'),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    panel.background = element_rect(fill = "white", colour = NA),
    #panel.grid = element_blank(),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=15, family="Times" ),
    plot.background = element_rect(colour = "white")
  )
}
myTheme_gillsansMT = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Gill Sans MT" ), #"Gill Sans MT"
    axis.line = element_line(colour = "#000000", size = 0.5, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=15, family="Gill Sans MT" ,colour='black'),
    legend.title = element_text(size=16, family="Gill Sans MT" ,colour='black'),
    axis.text.x = element_text(size=15, family="Gill Sans MT" ,colour='black'),
    axis.text.y = element_text(size=15, family="Gill Sans MT" ,colour='black'),
    axis.title.x=element_text(size=16, family="Gill Sans MT" ),
    axis.title.y=element_text(size=16, family="Gill Sans MT" ),
    plot.title = element_text(size=15, family="Gill Sans MT", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}
myTheme_gillsans = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Gill Sans" ), #"Gill Sans"
    axis.line = element_line(colour = "#000000", size = 0.5, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=15, family="Gill Sans" ,colour='black'),
    legend.title = element_text(size=16, family="Gill Sans" ,colour='black'),
    axis.text.x = element_text(size=15, family="Gill Sans" ,colour='black'),
    axis.text.y = element_text(size=15, family="Gill Sans" ,colour='black'),
    axis.title.x=element_text(size=16, family="Gill Sans" ),
    axis.title.y=element_text(size=16, family="Gill Sans" ),
    plot.title = element_text(size=15, family="Gill Sans", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}

myTheme_Times = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    axis.line = element_line(colour = "#000000", size = 0.5, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=15, family="Times" ,colour='black'),
    legend.title = element_text(size=16, family="Times" ,colour='black'),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    plot.title = element_text(size=15, family="Times", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}


myTheme_Helvetica = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=13, family="Helvetica" ), #"Helvetica"
    axis.line = element_line(colour = "#000000", size = 0.2, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=13, family="Helvetica" ,colour='black'),
    legend.title = element_text(size=13, family="Helvetica" ,colour='black'),
    axis.text.x = element_text(size=13, family="Helvetica" ,colour='black'),
    axis.text.y = element_text(size=13, family="Helvetica" ,colour='black'),
    axis.title.x=element_text(size=13, family="Helvetica" ),
    axis.title.y=element_text(size=13, family="Helvetica" ),
    plot.title = element_text(size=13, family="Helvetica", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}

myTheme_Arial = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=13, family="Arial" ), #"Arial"
    axis.line = element_line(colour = "#000000", size = 0.2, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=13, family="Arial" ,colour='black'),
    legend.title = element_text(size=13, family="Arial" ,colour='black'),
    axis.text.x = element_text(size=13, family="Arial" ,colour='black'),
    axis.text.y = element_text(size=13, family="Arial" ,colour='black'),
    axis.title.x=element_text(size=13, family="Arial" ),
    axis.title.y=element_text(size=13, family="Arial" ),
    plot.title = element_text(size=13, family="Arial", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}


## Converting function
convert_alpha_to_alphaRaw = function (alpha) {
    alpha[which(alpha < 1e-5)] <- 1e-5
    alpha[which(alpha > 1 - 1e-5)] <- 1 - 1e-5
    return (log(-alpha/(alpha-1)))
    # if(alpha != 0 & alpha != 1) {
    #     return (log(-alpha/(alpha-1)))
    # } else if (alpha == 0) {
    #     return (-10)
    # } else {
    #     return (10)
    # }
}

convert_alphaRaw_to_alpha = function (alphaRaw) {
    return( 1 / (1 + exp(-alphaRaw)))
}
