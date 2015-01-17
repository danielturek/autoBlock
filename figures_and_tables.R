

path <- '~/GitHub/autoBlock';     setwd(path)
source('autoBlock_utils.R')


## Figure: 'dfsampEff'
## Section on 'Efficiency loss from not blocking sampling'
## shows lines, attenuation of sampling efficiency, function of correlation, model size
path <- '~/GitHub/autoBlock';     setwd(path)
library(ggplot2); library(grid); library(gridExtra)
load('dfsampEff.RData')
df <- dfsampEff
#load('dfsampEffMostlyBlocked.RData')
load('dfsampEffExpDecay.RData')
dfExpDecay <- dfsampEff
df$d <- as.factor(df$N)
dfExpDecay$d <- as.factor(dfExpDecay$N)
p1 <- qplot(data=df, x=-log(1-rho), y=essPerN, color=d, geom='line', ylab=expression(S(Psi)), xlab=expression(-log(1-rho))) + theme(legend.position=c(.8, .8)) + scale_y_log10()
p2 <- qplot(data=dfExpDecay, x=-log(1-rho), y=essPerN, color=d, geom='line', ylab=expression(S(Psi)), xlab=expression(-log(1-rho))) + theme(legend.position=c(.8, .8)) + scale_y_log10()
dev.new(width=6, height=5)
grid.arrange(p1, p2, ncol = 2)
dev.copy2pdf(file='dfsampEff.pdf')
system('cp dfsampEff.pdf ~/GitHub/nimblePapers/autoBlock/')

## fitting regressions to these sampling effieincy lines
path <- '~/GitHub/autoBlock';     setwd(path)
load('dfsampEff.RData')
load('dfsampEffExpDecay.RData')
df <- dfsampEff[dfsampEff$k!=0,] # remove k=0, rho=0, where lines converge
qplot(data=df, x=log(1-rho), y=log(essPerN), color=factor(N), geom='line')
logESS <- log(df$essPerN)
logN <- log(df$N)
loglogN <- log(logN)
logRho <- log(1-df$rho)
## from dfsampEff
m <- lm(logESS ~ logN + logRho)
summary(m)
m$coeff
## eff = exp(intercept) * N^Bn * (1-rho)^Br
##                    (Bn)        (Br)
## (Intercept)        logN      logRho     ## from dfsampEff, with expDecay=FALSE,
##  -0.4346747  -1.1768437   1.0921665      ## alpha 0.2, rhos 0.8 0.96 ...
##
cbind(exp(-0.4346747) * dfsampEff$N^-1.1768437 * (1-dfsampEff$rho)^1.0921665, dfsampEff$essPerN)
##
## from dfsampEffExpDecay
m <- lm(logESS ~ loglogN + logRho)
summary(m)
m$coeff
## eff = exp(intercept) * (logN)^Bn * (1-rho)^Br
##                    (Bn)        (Br)
## (Intercept)     loglogN      logRho     ## from dfsampEffExpDecay, with expDecay=TRUE,
##   -1.261946   -1.222333    1.211014      ## alpha 0.2, rhos 0.8 0.96 ...



## Figure: 'blockTiming'
## Section on Numerical Timing Results of Computational Efficiency
## shows timing of algorithms, for different model structings, and all/non blocked
path <- '~/GitHub/autoBlock';     setwd(path)
library(ggplot2); library(grid); library(gridExtra)
load('dfblockTesting.RData')
df <- dfblockTesting[dfblockTesting$blocking != 'blockNoAdapt', ]  ## remove non-adaptive blocking
df$dist <- factor(df$dist, levels=c('uni','gamma','multi'), labels=c('Normal','Gamma','MV Normal'))  ## rename dist factor
df$blocking <- factor(df$blocking, levels=c('scalar','blockAdapt'), labels=c('All Scalar', 'All Blocked'))  ## rename blocking factor
dev.new(width=4, height=4.5)
qplot(data=df, x=N, y=timePer10kN, geom='line', linetype=dist, color=blocking, ylim=c(0,16), xlim=c(0,500)) + theme(legend.position=c(.47, .7)) + labs(x='Model dimension (d)', y='Runtime (seconds per 10,000 MCMC samples)', color='MCMC\nAlgorithm', linetype='Model\nStructure')
dev.copy2pdf(file='blockTiming.pdf')
system('cp blockTiming.pdf ~/GitHub/nimblePapers/autoBlock/')



## Figure 'contrivedMCMCefficiencyBars'
## Double bar-chart of Overall Efficiency for contrived model structures,
## combining two Simulated Data examples:
## left pane: 'partitions'
## right pane: 'mixedRhos'
path <- '~/GitHub/autoBlock';     setwd(path)
library(ggplot2); library(grid); library(gridExtra)
red<-'#D55E00';   green<-'#009E73';   pink='#CC79A7';   blue<-'blue';   lightblue<-'#56B4E9'
load('dfPartitionsN64.RData')
dfN64 <- dfN64[dfN64$mcmc != 'autoMax', ]  # don't include autoMax
dfN64$mcmc <- factor(dfN64$mcmc, levels=c('all','auto0','auto1','auto2'), labels=c('All Blocked','All Scalar','Auto 1','Auto 2'))
p1 <- ggplot(dfN64, aes(x=as.factor(rho),y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.71, .77)) + labs(y='Efficiency (effective samples / time)', x='Correlation', color='MCMC\nAlgorithm') + scale_color_manual(values=c(red,green,blue,lightblue))
load('dfmixedRhos.RData')
dfMix <- dfMix[dfMix$mcmc != 'autoMax', ]  # don't include autoMax
dfMix <- dfMix[dfMix$N %in% c(20, 50, 100), ]  # remove everything BUT N = 20,50,100
dfMix$mcmc <- factor(dfMix$mcmc, levels=c('all','auto0','auto1','auto2'), labels=c('All Blocked','All Scalar','Auto 1','Auto 2'))
p2 <- ggplot(dfMix, aes(x=as.factor(N),y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.70, .77)) + labs(y='Efficiency (effective samples / time)', x='Model Size (N)', color='MCMC\nAlgorithm') + scale_color_manual(values=c(red,green,blue,lightblue))
dev.new(width=6, height=4)
##multiplot(p1, p2, cols=2)
grid.arrange(p1, p2, ncol = 2)
dev.copy2pdf(file='contrivedMCMCefficiencyBars.pdf')
system('cp contrivedMCMCefficiencyBars.pdf ~/GitHub/nimblePapers/autoBlock/')




## Figure 'SSMLittersSpatialEfficiency'
## Line chart of Efficiency results for: SSM (both), Litters, and Spatial
path <- '~/GitHub/autoBlock';     setwd(path)
library(ggplot2); library(grid); library(gridExtra)
red<-'#D55E00';   green<-'#009E73';   black<-'black';   blue<-'blue';   lightblue<-'#56B4E9'
load('dfSSM.RData')
dfS <- dfS[c(1, 3:14), ]  ## remove the very-poor 'independent' model 'blockMUB' informed blocking
load('dflittersGAMMA-UNIFprior.RData')
load('dfspatial.RData')
dfFig <- rbind(dfS, dfLit, dfSpat)
dfFig[grepl('^block.*', dfFig$mcmc), ]$mcmc <- 'informed'  ## informed blockings to 'informed'
dfFig <- dfFig[dfFig$mcmc %in% c('all','auto0','default','autoMax'), ]  ## removed 'informed'
dfFig$model <- factor(dfFig$model, levels=c('independent', 'correlated', 'litters', 'spatial'), labels=c('State Space\nIndependent', 'State Space\nCorrelated', 'Random\nEffects', 'Spatial'))
dfFig$mcmc <- factor(dfFig$mcmc, levels=c('all', 'default', 'auto0', 'autoMax'), labels=c('All Blocked', 'Default', 'All Scalar', 'Auto Block'))
## normalize all Efficiencies by that of 'auto0'
## for(mod in unique(dfFig$model)) {
##     norm <- dfFig[dfFig$model==mod & dfFig$mcmc=='all', ]$Efficiency
##     dfFig[dfFig$model==mod, ]$Efficiency <- dfFig[dfFig$model==mod, ]$Efficiency / norm
## }
dev.new(width=4.5, height=4)
ggplot(dfFig, aes(x=model,y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.76, .73)) + labs(y='Efficiency (effective samples / time)', x='', color='MCMC Algorithm') + scale_color_manual(values=c(red,black,green,blue))
dev.copy2pdf(file='SSMLittersSpatialEfficiency.pdf')
system('cp SSMLittersSpatialEfficiency.pdf ~/GitHub/nimblePapers/autoBlock/')

## efficiency improvement factors
df <- dfFig
df <- df[df$model != 'independent',]
df <- df[df$mcmc %in% c('all', 'default', 'auto0', 'autoMax'), ]
for(mod in unique(df$model)) df[df$model==mod,]$Efficiency <- df[df$model==mod & df$mcmc=='autoMax',]$Efficiency / df[df$model==mod,]$Efficiency
df[, c('model', 'mcmc', 'Efficiency')]




## Figure 'spatialBoxplots'
## boxplots for some parameter distributions under the spatial model, various MCMCs
path <- '~/GitHub/autoBlock';     setwd(path)
library(ggplot2); library(grid); library(gridExtra)
load('dfspatialWithSamples.RData')   ## CURRENTLY, I DON'T HAVE THIS SAMPLES FILE
params <- c('mu', 'sigma', 'rho', 'g[66]')
df <- dfSpatSamples[dfSpatSamples$param %in% params, ]
figs <- list()
for(i in seq_along(params)) {
    param <- params[i]
    dfTemp <- df[df$param == param, ]
    figs[[i]] <- ggplot(dfTemp, aes(x=as.factor(param), y=samp, color=mcmc, ylab='', xlab=param)) + geom_boxplot() + theme(legend.position=c(.2, .8))
}
#dev.new(width=6, height=4)
figs$ncol = length(params)
do.call(grid.arrange, figs)

dev.new(width=6, height=4)
ggplot(df, aes(x=as.factor(param), y=samp, color=mcmc)) + geom_boxplot() + theme(legend.position=c(.2, .8))

dev.copy2pdf(file=NO_SAVE_YET)
system('cp ????? ~/GitHub/nimblePapers/autoBlock/')



## Figure for grant proposal
## 'mixedRhos'
## line chart, model-size on x-axis,
## efficiency of all, auto, none, on y-axis
## plot of efficiencies
path <- '~/GitHub/autoBlock';     setwd(path)
library(ggplot2)
load('dfmixedRhos.RData')
dfMix <- dfMix[dfMix$mcmc %in% c('all','auto0','autoMax'), ]  # only keep all blocked, no blocks (auto0), and best automatic performance
qplot(data=dfMix, x=N, y=Efficiency, color=mcmc, geom='line')
## this will give us the factor of improvement resulting from automatic blocking,
## for each value of N.
uniqueN <- unique(dfMix$N)
improvementFactor <- dfMix[dfMix$mcmc=='autoMax',]$Efficiency / unlist(lapply(uniqueN, function(N) max(dfMix[dfMix$N==N & dfMix$mcmc %in% c('all','auto0'),]$Efficiency)))
improvementFactor




## Table data for 'Varying Size Blocks of Fixed Correlation'
## This comes from 'partitions'
## path <- '~/GitHub/autoBlock';     setwd(path)
## load('dfpartitionsN64.RData')
## dfN64$essPer10k <- round(dfN64$essPer10k, 2)
## dfN64$essPer10k <- round(dfN64$essPer10k, 0)
## dfN64$Efficiency <- round(dfN64$Efficiency, 1)
## dfN64

## Table data for 'Fixed Size Blocks of Mixed Correlations'
## This comes from 'mixedRhos'
## path <- '~/GitHub/autoBlock';     setwd(path)
## load('dfmixedRhos.RData')
## dfMix$essPer10k <- round(dfMix$essPer10k, 2)
## dfMix$essPer10k <- round(dfMix$essPer10k, 0)
## dfMix$Efficiency <- round(dfMix$Efficiency, 1)
## dfMix

## Table data for State Space Models
path <- '~/GitHub/autoBlock';     setwd(path)
load('dfSSM.RData')
dfS$essPer10k <- round(dfS$essPer10k, 1)
dfS$timePer10k <- round(dfS$timePer10k, 2)
dfS$Efficiency <- round(dfS$Efficiency, 1)
dfS[, c(1, 2, 9, 8,10)]


## Table data for Litters
path <- '~/GitHub/autoBlock';     setwd(path)
load('dflittersGAMMA-UNIFprior.RData')
dfLit$essPer10k <- round(dfLit$essPer10k, 1)
dfLit$timePer10k <- round(dfLit$timePer10k, 2)
dfLit$Efficiency <- round(dfLit$Efficiency, 1)
dfLit[, c(1, 2, 9, 8, 10)]
## Table data for Spatial model
path <- '~/GitHub/autoBlock';     setwd(path)
load('dfspatial.RData')
dfSpat$essPer10k <- round(dfSpat$essPer10k, 1)
dfSpat$timePer10k <- round(dfSpat$timePer10k, 2)
dfSpat$Efficiency <- round(dfSpat$Efficiency, 1)
dfSpat[, c(1, 2, 9, 8, 10)]







