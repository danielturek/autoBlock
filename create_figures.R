



## Figure: 'samplingEfficiency'
## Section on 'Efficiency loss from not blocking sampling'
## shows lines, attenuation of sampling efficiency, function of correlation, model size
library(ggplot2); library(grid); library(gridExtra)
load('~/GitHub/autoBlock/results/results_samplingEfficiency.RData')
df <- dfSamplingEfficiency[dfSamplingEfficiency$expDecay==FALSE, ]
dfExpDecay <- dfSamplingEfficiency[dfSamplingEfficiency$expDecay==TRUE, ]
df$d <- as.factor(df$N)
dfExpDecay$d <- as.factor(dfExpDecay$N)
p1 <- qplot(data=df, x=-log(1-rho), y=essPerN, color=d, geom='line', ylab=expression(S(Psi)), xlab=expression(-log(1-rho))) + theme(legend.position=c(.8, .8)) + scale_y_log10()
p2 <- qplot(data=dfExpDecay, x=-log(1-rho), y=essPerN, color=d, geom='line', ylab=expression(S(Psi)), xlab=expression(-log(1-rho))) + theme(legend.position=c(.8, .8)) + scale_y_log10()
dev.new(width=6, height=5)
grid.arrange(p1, p2, ncol = 2)
dev.copy2pdf(file='~/GitHub/autoBlock/figures/samplingEfficiency.pdf')
system('cp ~/GitHub/autoBlock/figures/samplingEfficiency.pdf ~/GitHub/nimblePapers/autoBlock/')





## Figure: 'computationalRequirement'
## Section on Numerical Timing Results of Computational Efficiency
## shows timing of algorithms, for different model structings, and all/non blocked
library(ggplot2); library(grid); library(gridExtra)
load('~/GitHub/autoBlock/results/results_computationalRequirement.RData')
df <- dfComputationalRequirement
df <- df[df$blocking != 'blockNoAdapt', ]  ## remove non-adaptive blocking
df$dist <- factor(df$dist, levels=c('uni','gamma','multi'), labels=c('Normal','Gamma','MV Normal'))  ## rename dist factor
df$blocking <- factor(df$blocking, levels=c('scalar','blockAdapt'), labels=c('All Scalar', 'All Blocked'))  ## rename blocking factor
dev.new(width=4, height=4.5)
qplot(data=df, x=N, y=timePer10kN, geom='line', linetype=dist, color=blocking, ylim=c(0,16), xlim=c(0,500)) + theme(legend.position=c(.47, .7)) + labs(x='Model dimension (d)', y='Runtime (seconds per 10,000 MCMC samples)', color='MCMC\nAlgorithm', linetype='Model\nStructure')
dev.copy2pdf(file='~/GitHub/autoBlock/figures/computationalRequirement.pdf')
system('cp ~/GitHub/autoBlock/figures/computationalRequirement.pdf ~/GitHub/nimblePapers/autoBlock/')



## Figure 'contrivedModels'
## Double bar-chart of Overall Efficiency for contrived model structures,
## combining two Simulated Data examples:
## left pane: 'partitions'
## right pane: 'mixedRhos'
library(ggplot2); library(grid); library(gridExtra)
red<-'#D55E00';   green<-'#009E73';   pink='#CC79A7';   blue<-'blue';   lightblue<-'#56B4E9'
load('~/GitHub/autoBlock/results/results_varyingBlksFixedCorr.RData')
dfVary <- dfVaryingBlksFixedCorr_summary
dfVary <- dfVary[dfVary$mcmc %in% c('all', 'auto0', 'autoMax'), ]
dfVary$mcmc <- factor(dfVary$mcmc, levels=c('all','auto0','autoMax'), labels=c('All Blocked','All Scalar','Auto Blocking'))
p1 <- ggplot(dfVary, aes(x=as.factor(rho),y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.69, .81)) + labs(y='Efficiency (effective samples / time)', x='Correlation', color='MCMC\nAlgorithm') + scale_color_manual(values=c(red,green,blue,lightblue))
load('~/GitHub/autoBlock/results/results_fixedBlksVaryingCorr.RData')
dfFixed <- dfFixedBlksVaryingCorr_summary
dfFixed <- dfFixed[dfFixed$N %in% c(20, 50, 100), ]  # remove everything BUT N = 20,50,100
dfFixed$mcmc <- factor(dfFixed$mcmc, levels=c('all','auto0','autoMax'), labels=c('All Blocked','All Scalar','Auto Blocking'))
p2 <- ggplot(dfFixed, aes(x=as.factor(N),y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.68, .81)) + labs(y='Efficiency (effective samples / time)', x='Model size (N)', color='MCMC\nAlgorithm') + scale_color_manual(values=c(red,green,blue,lightblue))
dev.new(width=6, height=4)
##multiplot(p1, p2, cols=2)
grid.arrange(p1, p2, ncol = 2)
dev.copy2pdf(file='~/GitHub/autoBlock/figures/contrivedModels.pdf')
system('cp ~/GitHub/autoBlock/figures/contrivedModels.pdf ~/GitHub/nimblePapers/autoBlock/')




## Figure 'exampleModels'
## Line chart of Efficiency results for: SSM (both), Litters, and Spatial
library(ggplot2); library(grid); library(gridExtra)
red<-'#D55E00';   green<-'#009E73';   black<-'black';   blue<-'blue';   lightblue<-'#56B4E9'
load('~/GitHub/autoBlock/results/results_litters.RData')
load('~/GitHub/autoBlock/results/results_ice.RData')
load('~/GitHub/autoBlock/results/results_SSMindependent.RData')
load('~/GitHub/autoBlock/results/results_SSMcorrelated.RData')
load('~/GitHub/autoBlock/results/results_spatial.RData')
load('~/GitHub/autoBlock/results/results_mhp.RData')
dfFig <- rbind(dflitters_summary, dfice_summary, dfSSMindependent_summary, dfSSMcorrelated_summary, dfspatial_summary, dfmhp_summary)
dfFig[grepl('^block.*', dfFig$mcmc), ]$mcmc <- 'informed'  ## informed blockings to 'informed'
dfFig <- dfFig[dfFig$mcmc %in% c('all','auto0','default','autoMax'), ]  ## removed 'informed'
dfFig$model <- factor(dfFig$model, levels=c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp'), labels=c('Random\nEffects', 'Auto\nRegressive', 'St. Space\nIndep.', 'St. Space\nCorr.', 'Spatial', 'GLMM'))
dfFig$mcmc <- factor(dfFig$mcmc, levels=c('all', 'default', 'auto0', 'autoMax'), labels=c('All Blocked', 'Default', 'All Scalar', 'Auto Blocking'))
dev.new(width=5.0, height=4)
ggplot(dfFig, aes(x=model,y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.83, .76)) + labs(y='Efficiency (effective samples / time)', x='', color='MCMC Algorithm') + scale_color_manual(values=c(red,black,green,blue))
dev.copy2pdf(file='~/GitHub/autoBlock/figures/exampleModels.pdf')
system('cp ~/GitHub/autoBlock/figures/exampleModels.pdf ~/GitHub/nimblePapers/autoBlock/')



## efficiency improvement factors
df <- dfFig
##df <- df[df$model != 'independent',]
for(mod in unique(df$model)) df[df$model==mod,]$Efficiency <- df[df$model==mod & df$mcmc=='Auto Blocking',]$Efficiency / df[df$model==mod,]$Efficiency
df[, c('model', 'mcmc', 'Efficiency')]














