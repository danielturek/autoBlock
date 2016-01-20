

## before making exampleModels figure, we need to process and create the
## *_summary data frames......
## these should only need to be created *once*, after simulations are run
rm(list=ls())
## FAILED: hclust method 'median'
exampleModelNames <- c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp')
for(thisResultDir in c('results', 'results_hclust_single', 'results_hclust_average')) {
######for(thisResultDir in c('results', 'results_hclust_single', 'results_hclust_average', 'results_hclust_wardd')) {
    for(exModelName in exampleModelNames) {
        dataFileName <- paste0('results_', exModelName, '.RData')
        loadDir <- file.path('~/GitHub/legacy/autoBlock', thisResultDir, dataFileName)
        load(loadDir)
        if(thisResultDir == 'results') {
            thisDF <- eval(parse(text = paste0('df', exModelName, '_summary')))
            } else thisDF <- eval(parse(text = paste0('df', exModelName)))
        thisDF$model <- factor(exModelName)
        thisDF[thisDF$mcmc == 'All Blocked', 'mcmc'] <- 'all'
        thisDF[thisDF$mcmc == 'All Scalar', 'mcmc'] <- 'auto0'
        thisDF[thisDF$mcmc == 'Default', 'mcmc'] <- 'default'
        thisDF[thisDF$mcmc == 'Auto-Blocking', 'mcmc'] <- 'autoMax'
        eval(parse(text = paste0('df', exModelName, '_summary <- thisDF')))
        eval(parse(text = paste0('save(\'df', exModelName, '_summary\', file = \'', thisResultDir, '/results_', exModelName, '_summary.RData\')')))
    }
}


rm(list=ls())

## initial (baseline) runs, using hclust method 'complete'
resultsDirName <- 'results'
figFileNameTag <- ''

## next set of runs, for BA reviewers, using hclust method 'single'
resultsDirName <- 'results_hclust_single'
figFileNameTag <- '_hclust_single'



## Figure: 'algorithmicEfficiency'
## Section on 'Efficiency loss from not blocking sampling'
## shows lines, attenuation of algorithmic efficiency, function of correlation, model size
library(ggplot2); library(grid); library(gridExtra)
loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, 'results_samplingEfficiency.RData')
load(loadDir)
df <- dfsamplingEfficiency[dfsamplingEfficiency$expDecay==FALSE, ]
dfExpDecay <- dfsamplingEfficiency[dfsamplingEfficiency$expDecay==TRUE, ]
df$d <- as.factor(df$N)
dfExpDecay$d <- as.factor(dfExpDecay$N)
p1 <- qplot(data=df, x=-log(1-rho), y=essPerN, color=d, geom='line', ylab='Algorithmic Efficiency', xlab='Correlation Strength') + theme(legend.position=c(.14, .27)) + scale_y_log10()
p2 <- qplot(data=dfExpDecay, x=-log(1-rho), y=essPerN, color=d, geom='line', ylab='Algorithmic Efficiency', xlab='Correlation Strength') + theme(legend.position=c(.14, .27)) + scale_y_log10()
dev.new(width=6, height=3.5)
grid.arrange(p1, p2, ncol = 2)
thisFigureName <- 'algorithmicEfficiency'
figureFile <- paste0('~/GitHub/legacy/autoBlock/figures/', thisFigureName, figFileNameTag, '.pdf')
dev.copy2pdf(file=figureFile)
systemCopyCall <- paste0('cp ', figureFile, ' ~/GitHub/nimble/nimblePapers/autoBlock/')
system(systemCopyCall)





## Figure: 'computationalRequirement'
## Section on Numerical Timing Results of Computational Efficiency
## shows timing of algorithms, for different model structings, and all/non blocked
library(ggplot2); library(grid); library(gridExtra)
loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, 'results_computationalRequirement.RData')
load(loadDir)
df <- dfcomputationalRequirement
df <- df[df$blocking != 'blockNoAdapt', ]  ## remove non-adaptive blocking
df$dist <- factor(df$dist, levels=c('uni','gamma','multi'), labels=c('Normal','Gamma','MV Normal'))  ## rename dist factor
df$blocking <- factor(df$blocking, levels=c('scalar','blockAdapt'), labels=c('All Scalar', 'All Blocked'))  ## rename blocking factor
dev.new(width=4, height=4)
qplot(data=df, x=N, y=timePer10kN, geom='line', linetype=dist, color=blocking, ylim=c(0,16), xlim=c(0,500)) + theme(legend.position=c(.43, .65)) + labs(x='Model dimension (d)', y='Runtime (seconds per 10,000 MCMC samples)', color='MCMC\nAlgorithm', linetype='Model\nStructure')
thisFigureName <- 'computationalRequirement'
figureFile <- paste0('~/GitHub/legacy/autoBlock/figures/', thisFigureName, figFileNameTag, '.pdf')
dev.copy2pdf(file=figureFile)
systemCopyCall <- paste0('cp ', figureFile, ' ~/GitHub/nimble/nimblePapers/autoBlock/')
system(systemCopyCall)






## Figure 'contrivedModels'
## Double bar-chart of Overall Efficiency for contrived model structures,
## combining two Simulated Data examples:
## left pane: 'partitions'
## right pane: 'mixedRhos'
library(ggplot2); library(grid); library(gridExtra)
red<-'#D55E00';   green<-'#009E73';   pink='#CC79A7';   blue<-'blue';   lightblue<-'#56B4E9'
loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, 'results_varyingBlksFixedCorr.RData')
load(loadDir)
dfVary <- dfVaryingBlksFixedCorr_summary
dfVary <- dfVary[dfVary$mcmc %in% c('all', 'auto0', 'autoMax'), ]
dfVary$mcmc <- factor(dfVary$mcmc, levels=c('all','auto0','autoMax'), labels=c('All Blocked','All Scalar','Auto Blocking'))
p1 <- ggplot(dfVary, aes(x=as.factor(rho),y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.69, .81)) + labs(y='Efficiency (effective samples / time)', x='Correlation', color='MCMC\nAlgorithm') + scale_color_manual(values=c(red,green,blue,lightblue))
loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, 'results_fixedBlksVaryingCorr.RData')
load(loadDir)
dfFixed <- dfFixedBlksVaryingCorr_summary
dfFixed <- dfFixed[dfFixed$N %in% c(20, 50, 100), ]  # remove everything BUT N = 20,50,100
dfFixed$mcmc <- factor(dfFixed$mcmc, levels=c('all','auto0','autoMax'), labels=c('All Blocked','All Scalar','Auto Blocking'))
p2 <- ggplot(dfFixed, aes(x=as.factor(N),y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.68, .81)) + labs(y='Efficiency (effective samples / time)', x='Model size (N)', color='MCMC\nAlgorithm') + scale_color_manual(values=c(red,green,blue,lightblue))
dev.new(width=6, height=4)
##multiplot(p1, p2, cols=2)
grid.arrange(p1, p2, ncol = 2)
thisFigureName <- 'contrivedModels'
figureFile <- paste0('~/GitHub/legacy/autoBlock/figures/', thisFigureName, figFileNameTag, '.pdf')
dev.copy2pdf(file=figureFile)
systemCopyCall <- paste0('cp ', figureFile, ' ~/GitHub/nimble/nimblePapers/autoBlock/')
system(systemCopyCall)







## Figure 'exampleModels'
## Line chart of Efficiency results for: SSM (both), Litters, and Spatial
library(ggplot2); library(grid); library(gridExtra)
red<-'#D55E00';   green<-'#009E73';   black<-'black';   blue<-'blue';   lightblue<-'#56B4E9'
exampleModelNames <- c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp')
for(exName in exampleModelNames) {
    dataFileName <- paste0('results_', exName, '.RData')
    loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, dataFileName)
    load(loadDir)
}
dfFig <- rbind(dflitters_summary, dfice_summary, dfSSMindependent_summary, dfSSMcorrelated_summary, dfspatial_summary, dfmhp_summary)
dfFig[grepl('^block.*', dfFig$mcmc), ]$mcmc <- 'informed'  ## informed blockings to 'informed'
dfFig <- dfFig[dfFig$mcmc %in% c('all','auto0','default','autoMax'), ]  ## removed 'informed'
dfFig$model <- factor(dfFig$model, levels=c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp'), labels=c('Random\nEffects', 'Auto\nRegressive', 'St. Space\nIndep.', 'St. Space\nCorr.', 'Spatial', 'GLMM'))
dfFig$mcmc <- factor(dfFig$mcmc, levels=c('all', 'default', 'auto0', 'autoMax'), labels=c('All Blocked', 'Default', 'All Scalar', 'Auto Blocking'))
dev.new(width=4.5, height=3.5)
ggplot(dfFig, aes(x=model,y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position=c(.82, .74)) + labs(y='Efficiency (effective samples / time)', x='', color='MCMC Algorithm') + scale_color_manual(values=c(red,black,green,blue))
thisFigureName <- 'exampleModels'
figureFile <- paste0('~/GitHub/legacy/autoBlock/figures/', thisFigureName, figFileNameTag, '.pdf')
dev.copy2pdf(file=figureFile)
systemCopyCall <- paste0('cp ', figureFile, ' ~/GitHub/nimble/nimblePapers/autoBlock/')
system(systemCopyCall)


## make one plot with the performance of the 3 hclust methods:
## complete, single, average
rm(list=ls())
red<-'#D55E00';   green<-'#009E73';   black<-'black';   blue<-'blue';   lightblue<-'#56B4E9'
exampleModelNames <- c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp')
resultsDirName <- 'results'
for(exName in exampleModelNames) {
    dataFileName <- paste0('results_', exName, '_summary.RData')
    loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, dataFileName)
    load(loadDir)
}
dfFig <- rbind(dflitters_summary, dfice_summary, dfSSMindependent_summary, dfSSMcorrelated_summary, dfspatial_summary, dfmhp_summary)
dfFig$method <- 'complete'
dfFig <- dfFig[, c('model', 'mcmc', 'method', 'Efficiency')]
dfAllMethods <- dfFig
addOneHclustMethod <- function(name) {
    resultsDirName <- paste0('results_hclust_', name)
    for(exName in exampleModelNames) {
        dataFileName <- paste0('results_', exName, '_summary.RData')
        loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, dataFileName)
        load(loadDir)
    }
    dfFig <- rbind(dflitters_summary, dfice_summary, dfSSMindependent_summary, dfSSMcorrelated_summary, dfspatial_summary, dfmhp_summary)
    dfFig$method <- name
    dfFig <- dfFig[, c('model', 'mcmc', 'method', 'Efficiency')]
    dfAllMethods <<- rbind(dfAllMethods, dfFig)
}
## FAILED: hclust method 'median'
otherMethodsToAdd <- c('single', 'average')
##otherMethodsToAdd <- c('single', 'average', 'wardd')
for(method in otherMethodsToAdd)
    addOneHclustMethod(method)
dfAllMethods$method <- factor(dfAllMethods$method)
dfAllMethods$model <- factor(dfAllMethods$model, levels=c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp'), labels=c('Random\nEffects', 'Auto\nRegressive', 'St. Space\nIndep.', 'St. Space\nCorr.', 'Spatial', 'GLMM'))
dfAllMethods <- dfAllMethods[dfAllMethods$mcmc %in% c('all','auto0','default','autoMax'), ]
dfAllMethods$mcmc <- factor(dfAllMethods$mcmc, levels=c('all', 'default', 'auto0', 'autoMax'), labels=c('All Blocked', 'Default', 'All Scalar', 'Auto Blocking'))
dfAllMethods <- dfAllMethods[, c('model', 'mcmc', 'method', 'Efficiency')]
library(ggplot2)
ggplot(dfAllMethods[dfAllMethods$mcmc=='Auto Blocking',], aes(model, Efficiency, group=method, color=method)) + geom_line()
figureFile <- '~/GitHub/legacy/autoBlock/figures/hclustMethodComparison.pdf'
dev.copy2pdf(file=figureFile)



## efficiency reduction factors
## for the 'single' method
dfAllMethods[dfAllMethods$mcmc=='Auto Blocking' & dfAllMethods$method=='single', ]$Efficiency / dfAllMethods[dfAllMethods$mcmc=='Auto Blocking' & dfAllMethods$method=='complete', ]$Efficiency
##[1] 1.08940080 0.50979445 1.12148930 0.78565123 0.08720652 0.51997152
## for the 'average' method
dfAllMethods[dfAllMethods$mcmc=='Auto Blocking' & dfAllMethods$method=='average', ]$Efficiency / dfAllMethods[dfAllMethods$mcmc=='Auto Blocking' & dfAllMethods$method=='complete', ]$Efficiency
##[1] 1.05321486 0.49618926 1.23319142 0.75028737 0.08655572 0.49997262




## efficiency improvement factors
##df <- df[df$model != 'independent',]
for(mod in unique(df$model)) df[df$model==mod,]$Efficiency <- df[df$model==mod & df$mcmc=='Auto Blocking',]$Efficiency / df[df$model==mod,]$Efficiency
df[, c('model', 'mcmc', 'Efficiency')]








## NEW figures for presentations on autoBlock
## only showing All Blocked and Scalar, then also with AutoBlock
library(ggplot2); library(grid); library(gridExtra)
red<-'#D55E00';   green<-'#009E73';   black<-'black';   blue<-'blue';   lightblue<-'#56B4E9'
exampleModelNames <- c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp')
for(exName in exampleModelNames) {
    dataFileName <- paste0('results_', exName, '.RData')
    loadDir <- file.path('~/GitHub/legacy/autoBlock', resultsDirName, dataFileName)
    load(loadDir)
}
dfFig <- rbind(dflitters_summary, dfice_summary, dfSSMindependent_summary, dfSSMcorrelated_summary, dfspatial_summary, dfmhp_summary)
dfFig[grepl('^block.*', dfFig$mcmc), ]$mcmc <- 'informed'  ## informed blockings to 'informed'
dfFig$model <- factor(dfFig$model, levels=c('litters', 'ice', 'SSMindependent', 'SSMcorrelated', 'spatial', 'mhp'), labels=c('Random\nEffects', 'Auto\nRegressive', 'St. Space\nIndep.', 'St. Space\nCorr.', 'Spatial', 'GLMM'))
dfFig$mcmc <- factor(dfFig$mcmc, levels=c('auto0','all','autoMax'), labels=c('Univariate','All Blocked','Auto Blocking'))
hei <- 4
wid <- 6
dfFig3 <- dfFig[dfFig$mcmc %in% c('Univariate','All Blocked','Auto Blocking'), ]  ## all, none, auto
dfFig2 <- dfFig[dfFig$mcmc %in% c('Univariate','All Blocked'), ]  ## all, none
dev.new(width=wid, height=hei)
ggplot(dfFig3, aes(x=model,y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position='right') + labs(y='Efficiency (effective samples / time)', x='', color='MCMC Algorithm') + scale_color_manual(values=c(green,red,blue)) + ylim(0,40)
dev.copy2pdf(file='~/Downloads/autoBlock_exampleModels3.pdf')
dev.new(width=wid, height=hei)
ggplot(dfFig2, aes(x=model,y=Efficiency,group=mcmc,color=mcmc)) + geom_line() + geom_point(size=2.5) + theme(legend.position='right') + labs(y='Efficiency (effective samples / time)', x='', color='MCMC Algorithm') + scale_color_manual(values=c(green,red)) + ylim(0,40)
dev.copy2pdf(file='~/Downloads/autoBlock_exampleModels2.pdf')








