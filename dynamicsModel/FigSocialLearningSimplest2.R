###############################################################################################
##
## Population dynamics model
## Figure
## Wataru Toyokawa
## 5 May 2020
##
###############################################################################################

rm(list=ls(all=TRUE)) # cleaning the workspace

# Loading
library(readr)
library(tidyverse)

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

## Analytical solution of asocial behavioural dynamics
noSocialCurveSimplest= function (n, e, pH, pL) {
    return ( -(n*(pH - pL)*((-1 + e)*pH + e*pL))/((pH + pL)*((-1 + e)*pH - e*pL)) )
}

zeroIsoclineSimplest = function (pH, pL) {
    return ( pH/(pH + pL) )
}

wd <- "~/Dropbox/wataru/papers/RiskySocialLearning/populationDynamics"

####
## Equilibrium

sleqtableSimplest <- read_csv(paste0(wd, "/equilibrium/sleqtableSimplest.csv"), col_names = FALSE)
names(sleqtableSimplest) = c('f','S_0','c','e','R_eq')

sleqtableSimplest$RPreferingInitial = 0
sleqtableSimplest$direction = NA
sleqtableSimplest$R_0 = 20 - sleqtableSimplest$S_0
sleqtableSimplest$conformityExponent = paste(rep('f = ', nrow(sleqtableSimplest)), sleqtableSimplest$f, sep ='')
sleqtableSimplest$conformityExponent = factor(sleqtableSimplest$conformityExponent, levels = c('f = 0','f = 1','f = 2','f = 10'))
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
        x=expression(paste('Social influence ', italic(c), sep="")),
        y=expression(paste('Equilibrium density of ',N[R]^'*',sep=""))#,
        #title='Social influence\n (pH = 0.5; pL = 0.1; d=0.5; l=0.25)'
        )+
    scale_shape_manual(values=c('upward'=2, 'downward'=6), name='Stream\'s direction')+
    scale_color_manual(values=c('0'='#56B4E9','1'='#D55E00'), name='Risky choice regime')+
    myTheme_Times()+
    theme(axis.text.x = element_text(angle = 90), legend.position='top')+
    facet_grid(e_factor ~ conformityExponent)+
    xlim(c(0,1))+
    #ylim(c(0,20))+
    geom_hline(yintercept=10, linetype='dashed')+
    theme(legend.position = 'none')+
    theme(strip.text.y = element_text(angle = 0))+
    NULL -> sleqtableSimplest_plot)


ggsave(file = paste0(wd,"/equilibrium/sleqtableSimplest_plot.png"), plot = sleqtableSimplest_plot, dpi = 300, width = 9, height = 6)


####
## Equilibrium - l = 0.1; d = 0.2

sleqtableSimplest2 <- read_csv(paste0(wd,"/equilibrium/sleqtableSimplest2.csv"), col_names = FALSE)
names(sleqtableSimplest2) = c('f','S_0','c','e','R_eq')

sleqtableSimplest2$RPreferingInitial = 0
sleqtableSimplest2$direction = NA
sleqtableSimplest2$R_0 = 20 - sleqtableSimplest2$S_0
sleqtableSimplest2$conformityExponent = paste(rep('f = ', nrow(sleqtableSimplest2)), sleqtableSimplest2$f, sep ='')
sleqtableSimplest2$conformityExponent = factor(sleqtableSimplest2$conformityExponent, levels = c('f = 0','f = 1','f = 2','f = 10'))
sleqtableSimplest2$e_factor = paste(rep('e = ', nrow(sleqtableSimplest2)), sleqtableSimplest2$e, sep ='')

for (i in 1:nrow(sleqtableSimplest2)) {
    if (sleqtableSimplest2$R_eq[i]>10) {
        sleqtableSimplest2$RPreferingInitial[i] = 1
    }
    if (sleqtableSimplest2$R_eq[i] - sleqtableSimplest2$R_0[i] > 0) {
        sleqtableSimplest2$direction[i] = 'upward'
    } else {
        sleqtableSimplest2$direction[i] = 'downward'
    }
}

sleqtableSimplest2 <- sleqtableSimplest2 %>% dplyr::filter(c != 1)

(sleqtableSimplest2 %>%
    ggplot(aes(x=c))+
    geom_point(aes(y=R_0, colour=as.factor(RPreferingInitial), shape=direction), alpha=1/2)+
    geom_point(aes(y=R_eq))+
    labs(
        title = '',
        x=expression(paste('Social influence ', italic(c), sep="")),
        y=expression(paste('Equilibrium density of ',N[R]^'*',sep=""))#,
        #title='Social influence\n (pH = 0.5; pL = 0.1; d=0.5; l=0.25)'
        )+
    scale_shape_manual(values=c('upward'=2, 'downward'=6), name='Stream\'s direction')+
    scale_color_manual(values=c('0'='#56B4E9','1'='#D55E00'), name='Risky choice regime')+
    myTheme_Times()+
    theme(axis.text.x = element_text(angle = 90), legend.position='top')+
    facet_grid(e_factor ~ conformityExponent)+
    xlim(c(0,1))+
    #ylim(c(0,20))+
    geom_hline(yintercept=10, linetype='dashed')+
    theme(legend.position = 'none')+
    theme(strip.text.y = element_text(angle = 0))+
    NULL -> sleqtableSimplest2_plot)


ggsave(file = paste0(wd,"/equilibrium/sleqtableSimplest2_plot.png"), plot = sleqtableSimplest2_plot, dpi = 300, width = 9, height = 6)



######################################################
##
## Phase plot : e v.s. d (the heat-map against analytical curve)
##
######################################################

FigSocialLearningSimplest <- read_csv(paste0(wd,"/FigSocialLearningSimplest.csv"), col_names = FALSE)

names(FigSocialLearningSimplest) = c('e','c','f','pl','maxS','minS','diffS','maxR','minR','diffR','diffRS')
FigSocialLearningSimplest_sample = dplyr::sample_frac(tbl = FigSocialLearningSimplest, size = 0.00001)

FigSocialLearningSimplest$f_category = paste('f = ', FigSocialLearningSimplest$f, sep='')
FigSocialLearningSimplest$c_category = paste('c = ', FigSocialLearningSimplest$c, sep='')
FigSocialLearningSimplest$pl_category = paste('pl = ', FigSocialLearningSimplest$pl, sep='')
FigSocialLearningSimplest$f_category = factor(FigSocialLearningSimplest$f_category, levels = c('f = 1','f = 2','f = 10'))

(FigSocialLearningSimplest %>%
    ggplot(aes(x=pl))+
    #geom_line(mapping = aes(y = diffRS, group = pl, colour = pl)) +
    geom_raster(mapping = aes(y=e, fill = diffRS), stat = 'identity') +
    stat_function(fun=zeroIsoclineSimplest, args=list(pH=0.7), color='black', linetype='dashed', size=1)+
    labs(
        fill = expression(paste('Degree of risk seeking: ', N[R]^'*' - N[S]^'*',sep="")),
        y = expression(paste(italic(e), ': the rate of getting enchanted with R',sep="")),
        x = expression(paste('The rate of exploration: ', italic(p[l]), sep=""))
        )+
    scale_fill_gradient2(midpoint = 0, low = "blue", mid = "grey90", high = "red")+
    myTheme_Times()+
    theme(axis.text.x = element_text(angle = 90), legend.position='top')+
    theme(strip.text.y = element_text(angle = 0))+
    facet_grid(c_category ~ f_category)+
    #geom_hline(yintercept = 0, linetype = 'dotted')+
    NULL -> FigSocialLearningSimplest_plot)

#FigSocialLearningSimplest2 %>%
#    ggplot(aes(x=e))+
#    geom_line(mapping = aes(y = diffRS, group = pl, colour = pl)) +
#    #geom_tile(mapping = aes(d, l, fill = diffRS), stat = 'identity') +
#    stat_function(fun=noSocialCurveSimplest, args=list(n = 20, pH=0.7, pL=0.1), color='black', linetype='dashed', size=1)+
#    labs(
#        y = expression(paste('Risk seeking bias: ', N[R]^'*' - N[S]^'*',sep="")),
#        x = expression(paste(italic(e), ': the rate of getting enchanted with R',sep=""))
#        #title='Social influence\n (pH = 0.5; pL = 0.1; S0 = 10; R0 = 10)',
#        )+
#    scale_colour_viridis_c()+
#    myTheme_Times()+
#    theme(axis.text.x = element_text(angle = 90),legend.position='top')+
#    theme(strip.text.y = element_text(angle = 0))+
#    theme(legend.position = 'right')+
#    facet_grid(c_category ~ f_category)+
#    geom_hline(yintercept = 0, linetype = 'dotted')+
#    NULL -> FigSocialLearningSimplest_plot


ggsave(file = paste0(wd,"/FigSocialLearningSimplest_plot.png"), plot = FigSocialLearningSimplest_plot, dpi = 300, width = 6.5, height = 6.5)

#write.csv(FigSocialLearningSimplest_sample2, paste0(wd,"/FigSocialLearningSimplest_sample2.csv"), row.names=FALSE)


FigSocialLearningSimplestBiased <- read_csv(paste0(wd,"/FigSocialLearningSimplestBiased.csv"), col_names = FALSE)

names(FigSocialLearningSimplestBiased) = c('e','c','f','pl','maxS','minS','diffS','maxR','minR','diffR','diffRS')
FigSocialLearningSimplestBiased_sample = dplyr::sample_frac(tbl = FigSocialLearningSimplestBiased, size = 0.00001)

FigSocialLearningSimplestBiased$f_category = paste('f = ', FigSocialLearningSimplestBiased$f, sep='')
FigSocialLearningSimplestBiased$c_category = paste('c = ', FigSocialLearningSimplestBiased$c, sep='')
FigSocialLearningSimplestBiased$pl_category = paste('pl = ', FigSocialLearningSimplestBiased$pl, sep='')
FigSocialLearningSimplestBiased$f_category = factor(FigSocialLearningSimplestBiased$f_category, levels = c('f = 1','f = 2','f = 10'))

(FigSocialLearningSimplestBiased %>%
    ggplot(aes(x=pl))+
    #geom_line(mapping = aes(y = diffRS, group = pl, colour = pl)) +
    geom_raster(mapping = aes(y=e, fill = diffRS), stat = 'identity') +
    stat_function(fun=zeroIsoclineSimplest, args=list(pH=0.7), color='black', linetype='dashed', size=1)+
    labs(
        fill = expression(paste('Degree of risk seeking: ', N[R]^'*' - N[S]^'*',sep="")),
        y = expression(paste(italic(e), ': the rate of getting enchanted with R',sep="")),
        x = expression(paste('The rate of exploration: ', italic(p[l]), sep=""))
        )+
    scale_fill_gradient2(midpoint = 0, low = "blue", mid = "grey90", high = "red")+
    myTheme_Times()+
    theme(axis.text.x = element_text(angle = 90), legend.position='top')+
    theme(strip.text.y = element_text(angle = 0))+
    facet_grid(c_category ~ f_category)+
    #geom_hline(yintercept = 0, linetype = 'dotted')+
    NULL -> FigSocialLearningSimplestBiased_plot)


ggsave(file = paste0(wd,"/FigSocialLearningSimplestBiased_plot.png"), plot = FigSocialLearningSimplestBiased_plot, dpi = 300, width = 6.5, height = 6.5)

#write.csv(FigSocialLearningSimplestBiased_sample2, paste0(wd,"/FigSocialLearningSimplestBiased_sample2.csv"), row.names=FALSE)



###
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
       y = expression(paste('Risk-seeking bias: ', N[R]^'*' - N[S]^'*',sep="")),
       x = expression(paste(italic(e), ': the frequency of ', N[R]^'+',sep="")))+
    myTheme_Times()+
    theme(legend.position = 'none')+
    NULL)



ggsave(file = paste0(wd,"/noSocialCurveSimplest_plot.png"), plot = noSocialCurveSimplest_plot, dpi = 300, width = 3.5, height = 3.5)



