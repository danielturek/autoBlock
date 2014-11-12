preCode <- Quote({
    path <- '~/GitHub/autoBlock'
    control <- list(setSeed0 = TRUE)
    source(file.path(path, 'autoBlock_utils.R'))
})
eval(preCode)
preCode[[length(preCode)+1]] <- quote(control$makePlots <- FALSE)



## litters
littersCode <- quote({
    control$niter <- 400000
    runList <- list('all',
                    blockAB = list(c('a[1]','b[1]'), c('a[2]','b[2]')),
                    crossLevel = quote({
                        spec <- MCMCspec(Rmodel, nodes=NULL)
                        spec$addSampler('crossLevel', list(topNodes = c('a[1]', 'b[1]')), print=FALSE)
                        spec$addSampler('crossLevel', list(topNodes = c('a[2]', 'b[2]')), print=FALSE)
                        spec
                    }),
                    'default',
                    'auto')
    ablitters <- autoBlock(code=code_litters, constants=constants_litters, data=data_litters, inits=inits_litters, control=control)
    ablitters$run(runList)
    abList <- list(litters=ablitters)
    dflitters <- createDFfromABlist(abList)
    filename <- file.path(path, 'dflittersGAMMA-UNIFprior.RData')
    save(dflitters, file = filename)
    if(ablitters$makePlots) plotABS(dflitters, xlimToMin=FALSE)
    if(ablitters$makePlots) plotABS(dflitters, xlimToMin=TRUE)
    printMinTimeABS(dflitters)
})
filename <- file.path(path, 'runLitters.R')
cat(codeToText(preCode), file=filename)
cat(codeToText(littersCode), file=filename, append=TRUE)



## state space models
SSMCode <- quote({
    control$niter <- 400000
    runListMUB <- list('all', blockMUB = list(c('mu','b')), 'default', 'auto')
    runListAB  <- list('all', blockAB  = list(c('a', 'b')), 'default', 'auto')
    abSSMmub <- autoBlock(code=code_SSMmub, constants=constants_SSMmub, data=data_SSMmub, inits=inits_SSMmub, control=control)
    abSSMab <- autoBlock(code=code_SSMab, constants=constants_SSMab, data=data_SSMab, inits=inits_SSMab, control=control)
    abSSMmub$run(runListMUB)
    abSSMab$run(runListAB)
    abList <- list(independent=abSSMmub, correlated=abSSMab)
    dfSSM <- createDFfromABlist(abList)
    filename <- file.path(path, 'dfSSM.RData')
    save(dfSSM, file = filename)
    if(abSSMab$makePlots) plotABS(dfSSM, xlimToMin=FALSE)
    if(abSSMab$makePlots) plotABS(dfSSM, xlimToMin=TRUE)
    printMinTimeABS(dfSSM)
})
filename <- file.path(path, 'runSSM.R')
cat(codeToText(preCode), file=filename)
cat(codeToText(SSMCode), file=filename, append=TRUE)



## spatial
spatialCode <- quote({
    control$niter <- 400000
    runList <- list('all', 'default', 'auto')
    abspatial <- autoBlock(code=code_spatial, constants=constants_spatial, data=data_spatial, inits=inits_spatial, control=control)
    abspatial$run(runList)
    abList <- list(spatial=abspatial)
    dfspatial <- createDFfromABlist(abList)
    filename <- file.path(path, 'dfspatial.RData')
    save(dfspatial, file = filename)
    if(abspatial$makePlots) plotABS(dfspatial, xlimToMin=FALSE)
    if(abspatial$makePlots) plotABS(dfspatial, xlimToMin=TRUE)
    printMinTimeABS(dfspatial)
})
filename <- file.path(path, 'runspatial.R')
cat(codeToText(preCode), file=filename)
cat(codeToText(spatialCode), file=filename, append=TRUE)



## partitions of N = 2^k
## constant rho
k <- 4
rhoVector <- c(0.2, 0.5, 0.8)
for(rho in rhoVector) {
    partitionsCode <- substitute({
        rho <- RHO
        k <- KKK
        N <- 2^k
        tag <- paste0('N', N, 'rho', rho)
        control$niter <- 100000
        runList <- list('all', 'auto')
        blockSizes <- 2^(0:k)
        data <- list()
        inits <- list(x=rep(0,N))
        abList <- list()
        for(blockSize in blockSizes) {
            numberOfBlocks <- N / blockSize
            listOfBlockIndexes <- lapply(((1:numberOfBlocks)-1)*blockSize, function(x) x+(1:blockSize))
            codeAndConstants <- createCodeAndConstants(N, listOfBlockIndexes, rep(rho,numberOfBlocks))
            code <- codeAndConstants$code
            constants <- codeAndConstants$constants
            ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
            ab$run(runList)
            abList[[paste0('blockSz', blockSize)]] <- ab
        }
        if(k > 1) {
            blockLengths <- c(1, 2^(0:(k-1)))
            indList <- list(); cur <- 1
            for(len in blockLengths) { indList <- c(indList, list(cur:(cur+len-1))); cur <- cur+len }
            codeAndConstants <- createCodeAndConstants(N, indList, rep(rho,length(indList)))
            code <- codeAndConstants$code
            constants <- codeAndConstants$constants
            ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
            ab$run(runList)
            abList[['blockSzMixed']] <- ab
        }
        dfText <- paste0('df', tag)
        eval(substitute(DF <- createDFfromABlist(abList), list(DF=as.name(dfText))))
        filename <- file.path(path, paste0(dfText, '.RData'))
        eval(substitute(save(DF, file = filename), list(DF=as.name(dfText))))
        if(ab$makePlots) eval(substitute(plotABS(DF), list(DF=as.name(dfText))))
        eval(substitute(printMinTimeABS(DF), list(DF=as.name(dfText))))
    }, list(RHO=rho, KKK=k))
    filename <- file.path(path, paste0('runPartitionsN', 2^k, 'rho', rho, '.R'))
    cat(codeToText(preCode), file=filename)
    cat(codeToText(partitionsCode), file=filename, append=TRUE)
}



## mixed, overlapping, rhos
mixedRhosCode <- substitute({
    control$niter <- 400000
    abList <- list()
    Nvalues <- c(20, 30, 40, 50, 100)   ## multiples of 10
    for(N in Nvalues) {
        tag <- paste0('mixedRhosN', N)
        blockSize <- N/10
        numberOfBlocks <- 9
        indList <- lapply(((1:numberOfBlocks)-1)*blockSize, function(x) x+(1:blockSize))
        rhoVector <- seq(from=0.9, to=0.1, by=-0.1)
        runList <- list('all', 'auto')
        codeAndConstants <- createCodeAndConstants(N, indList, rhoVector)
        code <- codeAndConstants$code
        constants <- codeAndConstants$constants
        data <- list()
        inits <- list(x=rep(0,N))
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[[tag]] <- ab
    }
    dfText <- 'dfmixedRhos'
    eval(substitute(DF <- createDFfromABlist(abList), list(DF=as.name(dfText))))
    filename <- file.path(path, paste0(dfText, '.RData'))
    eval(substitute(save(DF, file = filename), list(DF=as.name(dfText))))
    if(control$makePlots) eval(substitute(plotABS(DF), list(DF=as.name(dfText))))
    eval(substitute(printMinTimeABS(DF), list(DF=as.name(dfText))))
})
filename <- file.path(path, 'runmixedRhos.R')
cat(codeToText(preCode), file=filename)
cat(codeToText(mixedRhosCode), file=filename, append=TRUE)




## scalar or block samplers, for various N
scalarOrBlockCode <- substitute({
    control$niter <- 400000
    abList <- list()
    ## Nvalues <- c(2, 3, 4, 5, 10, 20, 50)
    Nvalues <- c(100, 200, 500)
    for(N in Nvalues) {
        tag <- paste0('scalarOrBlockN', N)
        code <- createCodeAndConstants(N)$code
        constants <- list()
        data <- list()
        inits <- list(x=rep(0,N))
        runList <- list(
            scalar = quote({
                spec <- MCMCspec(Rmodel, nodes=NULL)
                for(node in Rmodel$expandNodeNames('x')) spec$addSampler('RW', list(targetNode=node), print=FALSE)
                spec
            }),
            blockNoAdapt = quote({
                spec <- MCMCspec(Rmodel, nodes=NULL)
                spec$addSampler('RW_block', list(targetNodes=Rmodel$expandNodeNames('x'), adaptScaleOnly=TRUE), print=FALSE)
                spec
            }),
            blockAdaptive = quote({
                spec <- MCMCspec(Rmodel, nodes=NULL)
                spec$addSampler('RW_block', list(targetNodes=Rmodel$expandNodeNames('x')), print=FALSE)
                spec
            }))
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[[tag]] <- ab
    }
    dfText <- 'dfscalarOrBlock'     # optionally: 'dfscalarOrBlockTEMP'
    eval(substitute(DF <- createDFfromABlist(abList), list(DF=as.name(dfText))))
    filename <- file.path(path, paste0(dfText, '.RData'))
    eval(substitute(save(DF, file = filename), list(DF=as.name(dfText))))
    if(control$makePlots) eval(substitute(plotABS(DF), list(DF=as.name(dfText))))
    eval(substitute(printMinTimeABS(DF), list(DF=as.name(dfText))))
    })
filename <- file.path(path, 'runscalarOrBlock.R')
cat(codeToText(preCode), file=filename)
cat(codeToText(scalarOrBlockCode), file=filename, append=TRUE)

## load('dfscalarOrBlock.RData')
## head(dfscalarOrBlock)
## dim(dfscalarOrBlock)
## unique(dfscalarOrBlock$model)
## printMinTimeABS(dfscalarOrBlock)
## load('dfscalarOrBlockTEMP.RData')
## head(dfscalarOrBlockTEMP)
## dim(dfscalarOrBlockTEMP)
## unique(dfscalarOrBlockTEMP$model)
## printMinTimeABS(dfscalarOrBlockTEMP)
## dfscalarOrBlock <- rbind(dfscalarOrBlock, dfscalarOrBlockTEMP)
## head(dfscalarOrBlock)
## dim(dfscalarOrBlock)
## unique(dfscalarOrBlock$model)
## save(dfscalarOrBlock, file = 'dfscalarOrBlock.RData')
## rm(list=ls())
## load('dfscalarOrBlock.RData')
## head(dfscalarOrBlock)
## dim(dfscalarOrBlock)
## unique(dfscalarOrBlock$model)
## printMinTimeABS(dfscalarOrBlock)
## plotABS(dfscalarOrBlock)

testBlockSamplerPerformance <- function(N, niter) {
    code <- createCodeAndConstants(N)$code
    constants <- list()
    data <- list()
    inits <- list(x=rep(0,N))
    Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
    spec <- MCMCspec(Rmodel, NULL)
    spec$addSampler('RW_block', list(targetNodes='x', adaptScaleOnly=TRUE), print=FALSE)
    Rmcmc <- buildMCMC(spec)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    set.seed(0); Cmcmc(niter)
    samples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))
    burnedSamples <- samples[(niter/2+1):niter, , drop = FALSE]
    ess <- apply(burnedSamples, 2, effectiveSize)
    meanESS <- mean(ess)
    eff <- meanESS / (niter/2)
    adaptedScale <- nfVar(Cmcmc, 'samplerFunctions')$contentsList[[1]]$scale
    acceptanceRateHistory <- nfVar(Cmcmc, 'samplerFunctions')$contentsList[[1]]$acceptanceRateHistory
    aRate <- acceptanceRateHistory[length(acceptanceRateHistory)]
    optimalRates <- c(0.44, 0.35, 0.32, 0.25, 0.234)
    optRate <- optimalRates[if(N>5) 5 else N]
    retDF <- data.frame(
        N = N,
        eff = eff,
        adaptedScale = adaptedScale,
        aRate = aRate,
        optRate = optRate
    )
    return(retDF)
}

Nvalues <- 1:10
niter <- 1000000
dfblockSamplerPerformance <- data.frame()
for(i in Nvalues) {
    dfout <- testBlockSamplerPerformance(N = i, niter = niter)
    dfblockSamplerPerformance <- rbind(dfblockSamplerPerformance, dfout)
}
save(dfblockSamplerPerformance, file = 'dfblockSamplerPerformance.RData')
df <- dfblockSamplerPerformance
df
y <- sqrt(df$N) * df$adaptedScale  ###### Interesting ! ! !
y
plot(y, ylim=c(2,3))
y <- df$meanESS * df$N
plot(y)





## testing of litters priors
if(FALSE) {
    
    library(nimble)
    library(coda)

    suite <- MCMCsuite(
        model = code_litters,
        constants = constants_litters,
        data = data_litters,
        inits = inits_litters,
        monitors = c('a', 'b', 'p'),
        niter = 200000,
        thin = 1,
        burnin = 100000,
        MCMCs = c('ABfirst', 'ABlast', 'crossLevel'),
        MCMCdefs = list(
            ABfirst = quote({
                spec <- MCMCspec(Rmodel, useConjugacy = FALSE)
                spec$addSampler('RW_block', list(targetNodes = c('a[1]', 'b[1]')), print=FALSE)
                spec$addSampler('RW_block', list(targetNodes = c('a[2]', 'b[2]')), print=FALSE)
                spec$setSamplers(c(37,38,5:36), print=FALSE)
                spec$getSamplers()
                spec
            }),
            ABlast = quote({
                spec <- MCMCspec(Rmodel, useConjugacy = FALSE)
                spec$addSampler('RW_block', list(targetNodes = c('a[1]', 'b[1]')), print=FALSE)
                spec$addSampler('RW_block', list(targetNodes = c('a[2]', 'b[2]')), print=FALSE)
                spec$setSamplers(c(5:36,37,38), print=FALSE)
                spec$getSamplers()
                spec
            }),
            crossLevel = quote({
                spec <- MCMCspec(Rmodel, nodes=NULL)
                spec$addSampler('crossLevel', list(topNodes = c('a[1]', 'b[1]')), print=FALSE)
                spec$addSampler('crossLevel', list(topNodes = c('a[2]', 'b[2]')), print=FALSE)
                spec$getSamplers()
                spec
            })
            ),
        summaryStats = c('mean', 'sd', 'function(x) effectiveSize(x)'),
        makePlot = FALSE,
        savePlot = FALSE
        )

    suite$output$summary

    samples <- suite$output$samples
    
    littersTraceplots <- function(node, samples) {
        num <- dim(samples)[1]
        par(mfrow = c(num, 1))
        tsplot <- function(x) plot(1:length(x), x, type='l')
        for(i in 1:num) tsplot(samples[i, node, ])
    }

    littersTraceplots('p[1, 15]', samples)
    littersTraceplots('p[1, 8]', samples)
    littersTraceplots('a[1]', samples)   # 80,000
    littersTraceplots('b[1]', samples)   # 10,000
    littersTraceplots('a[2]', samples)   # 100  (200 too high)
    littersTraceplots('b[2]', samples)   # 50 (100 to high)

}




## path <- '~/GitHub/autoBlock'
## control <- list(setSeed0 = TRUE, makePlots = FALSE, niter = 2000)
## source(file.path(path, 'autoBlock_utils.R'))
## abtester <- autoBlock(code=code_tester, constants=constants_tester, data=data_tester, inits=inits_tester, control=control)
## runList <- list('all', 'auto', 'default')
## abtester$run(runList)
## abList <- list(tester = abtester)
## df <- createDFfromABlist(abList)
## df
## printMinTimeABS(df)





## setwd('~/GitHub/autoBlock/')
## source('autoBlock_utils.R')
## library(Imap)
## myscallops <- read.table("http://www.biostat.umn.edu/~brad/data/myscallops.txt", header = TRUE)
## myscallops <- myscallops[1:10,]
## N <- dim(myscallops)[1]
## catch <- myscallops$tcatch
## lat <- myscallops$lat
## long <- myscallops$long
## dist <- array(NA, c(N,N))
## for(i in 1:N) for(j in 1:N) dist[i,j] <- gdist(long[i], lat[i], long[j], lat[j])
## code_spatial<- modelCode({
##     mu ~ dunif(-100, 100)
##     sigma ~ dunif(0, 100)
##     rho ~ dunif(10, 500)
##     muVec[1:N] <- mu * onesVector[1:N]
##     Cov[1:N,1:N] <- sigma^2 * exp(-dist[1:N,1:N] / rho)
##     g[1:N] ~ dmnorm(muVec[1:N], cov = Cov[1:N,1:N])
##     for(i in 1:N) {
##         y[i] ~ dpois(exp(g[i]))
##     }
## })
## constants_spatial <- list(N=N, onesVector=rep(1,N), dist=dist)
## data_spatial <- list(y=catch)
## inits_spatial <- list(mu=2, sigma=5, rho=60, g=rep(0,N))

## suite <- MCMCsuite(
##     model = code_spatial,
##     constants = constants_spatial,
##     data = data_spatial,
##     inits = inits_spatial,
##     niter = 50000,
##     MCMCs = c('custom'),
##     MCMCdefs = list(
##         custom = quote({
##             spec <- MCMCspec(Rmodel, nodes=NULL)
##             spec$addSampler('RW_block', list(targetNodes = c('rho', 'sigma')), print=FALSE)
##             spec$addSampler('RW', list(targetNode = 'mu'), print=FALSE)
##             for(kk in 1:N) { spec$addSampler('RW', list(targetNode = paste0('g[',kk,']')), print=FALSE) }
##             spec$getSamplers()
##             spec
##         })),
##     savePlot = FALSE
## )
## suite$output$summary



