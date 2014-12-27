preCode <- quote({
    path <- '~/GitHub/autoBlock'
    control <- list(setSeed0 = TRUE)
    source(file.path(path, 'autoBlock_utils.R'))
})
eval(preCode)
preCode[[length(preCode)+1]] <- quote(control$makePlots <- FALSE)



## shows how bad sampling efficiency can be, as correlation increases
## output used in Figure 'dfsampEff'
sampEffCode <- quote({
    kValues <- 0:3
    Nvalues <- c(2, 4, 8, 16, 32)
    sampOption <- 1
    expDecay <- TRUE
    niter <- 4000000
    keepInd <- (niter/2+1):niter
    dfsampEff <- data.frame()
    for(k in kValues) {
        for(N in Nvalues) {
            rho <- 1 - (1-0.8)^k
            cat(paste0('\nk = ', k, '\nrho = ', rho, '\nN = ', N, '\n\n'))
            candc <- createCodeAndConstants(N, list(1:N), rho, expDecay=expDecay)
            print(candc)
            code <- candc$code
            constants <- candc$constants
            data <- list()
            inits <- list(x = rep(0,N))
            Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
            nodeNames <- Rmodel$expandNodeNames('x', returnScalarComponents = TRUE)
            spec <- MCMCspec(Rmodel, nodes = NULL)
            if(sampOption==1) { # scalar RW samplers on all nodes
                for(node in nodeNames) spec$addSampler('RW', list(targetNode=node), print=FALSE)
                spec$getSamplers() }
            if(sampOption==2) { # scalar RW sampler on x[1]; block sampler on x[2:N]
                spec$addSampler('RW', list(targetNode='x[1]'), print=FALSE)
                spec$addSampler('RW_block', list(targetNodes=nodeNames[-1]), print=FALSE)
                spec$getSamplers() }
            Rmcmc <- buildMCMC(spec)
            compiledList <- compileNimble(list(Rmodel, Rmcmc))
            Cmodel <- compiledList[[1]]; Cmcmc <- compiledList[[2]]
            Cmodel$setInits(inits)
            set.seed(0)
            timing <- as.numeric(system.time(Cmcmc(niter))[1])
            timePer10kN <- timing / (niter/10000)
            samples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))
            samples <- samples[keepInd, , drop = FALSE]
            ess <- apply(samples, 2, effectiveSize)
            if(sampOption==1) meanESS <- mean(ess)
            if(sampOption==2) meanESS <- as.numeric(ess[1])
            essPerN <- meanESS / length(keepInd)
            samples <- NULL; Cmcmcs <- NA; gc()
            thisDF <- data.frame(k=k, rho=rho, N=N, timePer10kN=timePer10kN, essPerN=essPerN)
            dfsampEff <- rbind(dfsampEff, thisDF)
            ##save(dfsampEff, file = paste0('dfsampEff.RData'))
            save(dfsampEff, file = paste0('dfsampEffExpDecay.RData'))
        }
    }
})
filename <- file.path(path, paste0('runsampEffExpDecay.R'))
cat(codeToText(preCode), file=filename)
cat(codeToText(sampEffCode), file=filename, append=TRUE)




## assesses the adapted scale, acceptance rates, ESS, and timing
## achieved by scalar/block samplers of various sizes, and underlying
## univariate or multivariate distributions
## used in Figure: 'blockTiming'
tagValues <- LETTERS[1:13]
for(tag in tagValues) {
    blockTestingCode <- substitute({
        tag <- TAG
        switch(tag,
               A = { dist <- 'gamma'; Nvalues <- c(2, 3, 4, 5, 10, 20, 30, 40, 50) },
               B = { dist <- 'gamma'; Nvalues <- c(100, 150) },
               C = { dist <- 'gamma'; Nvalues <- c(200, 250) },
               D = { dist <- 'gamma'; Nvalues <- c(300) },
               E = { dist <- 'gamma'; Nvalues <- c(350) },
               F = { dist <- 'gamma'; Nvalues <- c(400) },
               G = { dist <- 'gamma'; Nvalues <- c(450) },
               H = { dist <- 'gamma'; Nvalues <- c(500) },
               I = { dist <- 'gamma'; Nvalues <- c(600) },
               J = { dist <- 'gamma'; Nvalues <- c(700) },
               K = { dist <- 'gamma'; Nvalues <- c(800) },
               L = { dist <- 'gamma'; Nvalues <- c(900) },
               M = { dist <- 'gamma'; Nvalues <- c(1000) }
               )
        niter <- 50000
        keepInd <- (niter/2+1):niter
        dfblockTesting <- data.frame()
        for(N in Nvalues) {
            cat(paste0('\nN = ', N, '\n'))
            cat(paste0('\ndist = ', dist, '\n'))
            if(dist == 'uni')   canda <- createCodeAndConstants(N)
            if(dist == 'multi') candc <- createCodeAndConstants(N, list(1:N), 0)
            if(dist == 'gamma') candc <- createCodeAndConstants(N, gammaScalars = TRUE)
            code <- candc$code
            constants <- candc$constants
            print(code)
            print(constants)
            data <- list()
            inits <- list(x = rep(1,N))
            cat('\ncreating R model.....\n')
            Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
            nodeNames <- Rmodel$expandNodeNames('x', returnScalarComponents = TRUE)
            specList <- list()  # ordering: scalar, blockNoAdapt, blockAdapt
            for(i in 1:3) specList[[i]] <- MCMCspec(Rmodel, nodes = NULL)
            for(node in nodeNames) specList[[1]]$addSampler('RW', list(targetNode=node), print=FALSE)
            specList[[2]]$addSampler('RW_block', list(targetNodes=nodeNames, adaptScaleOnly=TRUE), print=FALSE)
            specList[[3]]$addSampler('RW_block', list(targetNodes=nodeNames), print=FALSE)
            toCompileList <- list(Rmodel)
            for(i in 1:3) toCompileList[[i+1]] <- buildMCMC(specList[[i]])
            cat('\ncompiling.....\n\n')
            compiledList <- compileNimble(toCompileList)
            cat('\ndone compiling!\n\n')
            Cmodel <- compiledList[[1]]
            Cmcmcs <- compiledList[2:4]  # ordering: scalar, blockNoAdapt, blockAdapt
            timePer10kN <- adaptedScale <- adaptedPropSD <- essPerN <- numeric(0)
            for(i in 1:3) {
                cat(paste0('running ', i, '\n'))
                Cmodel$setInits(inits)
                set.seed(0)
                timing <- as.numeric(system.time(Cmcmcs[[i]](niter))[1])
                timePer10kN[i] <- timing / (niter/10000)
                sampler1 <- nfVar(Cmcmcs[[i]], 'samplerFunctions')$contentsList[[1]]
                adaptedScale[i] <- sampler1$scale
                adaptedPropSD[i] <- if(i==1) as.numeric(NA) else sqrt(mean(diag(sampler1$propCov)))
                samples <- as.matrix(nfVar(Cmcmcs[[i]], 'mvSamples'))
                samples <- samples[keepInd, , drop = FALSE]
                ess <- apply(samples, 2, effectiveSize)
                meanESS <- mean(ess)
                essPerN[i] <- meanESS / length(keepInd)
                samples <- NULL
                sampler1 <- NULL
                Cmcmcs[[i]] <- NA
                gc()
            }
            thisDF <- data.frame(
                N = rep(N, 3),
                dist = rep(dist, 3),
                blocking = c('scalar', 'blockNoAdapt', 'blockAdapt'),
                timePer10kN = timePer10kN,
                derivedScale = c(adaptedScale[1], adaptedScale[2:3] * adaptedPropSD[2:3]),
                essPerN = essPerN
            )
            dfblockTesting <- rbind(dfblockTesting, thisDF)
            save(dfblockTesting, file = paste0('dfblockTesting', TAG, '.RData'))
            cat('\n'); print(dfblockTesting)
            gc()
        }
    },list(TAG = tag))
    filename <- file.path(path, paste0('runblockTesting', tag, '.R'))
    cat(codeToText(preCode), file=filename)
    cat(codeToText(blockTestingCode), file=filename, append=TRUE)
}
filename <- 'runblockTesting.sh'
cat('#!/bin/bash\n\n', file=filename)
for(tag in tagValues) {
    cat(paste0('R CMD BATCH --vanilla runblockTesting', tag, '.R\n'), file=filename, append=TRUE)
    cat('git add --all\n', file=filename, append=TRUE)
    cat(paste0('git commit -a -m\'ran runBlockTesting', tag, '.R\'\n'), file=filename, append=TRUE)
    cat('git push\n\n', file=filename, append=TRUE)
}
system(paste0('chmod 777 ', filename))

## combining the A, B, C, ...  dataframes from blockTesting
rm(list=ls())
dfCombined <- data.frame()
tagValues <- LETTERS[1:13]
for(tag in tagValues) {
    load(paste0('dfblockTesting', tag, '.RData'))
    dfCombined <- rbind(dfCombined, dfblockTesting)
}
load('dfblockTesting.RData')
dfblockTestingGamma <- dfCombined
save(dfblockTesting, dfblockTestingUni, dfblockTestingMulti, dfblockTestingGamma,
     file = 'dfblockTesting.RData')



## 'partitions' Simulated Data example
## N = 2^k, constant values of rho
## Used in a Figure of Simulated Data results, and probably a table
partitionsCode <- quote({
    k <- 6
    N <- 2^k
    rhoVector <- c(0.2, 0.5, 0.8)
    control$niter <- 200000
    abList <- list()
    for(rho in rhoVector) {
        tag <- paste0('N', N, 'rho', rho)
        runList <- list('all', 'auto')
        blockLengths <- c(1, 2^(0:(k-1)))
        indList <- list(); cur <- 1
        for(len in blockLengths) { indList <- c(indList, list(cur:(cur+len-1))); cur <- cur+len }
        data <- list()
        inits <- list(x=rep(0,N))
        codeAndConstants <- createCodeAndConstants(N, indList, rep(rho,length(indList)))
        code <- codeAndConstants$code
        constants <- codeAndConstants$constants
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[[paste0('blockSzMixed',tag)]] <- ab
    }
    dfText <- paste0('dfPartitionsN', N)
    eval(substitute(DF <- createDFfromABlist(abList), list(DF=as.name(dfText))))
    filename <- file.path(path, paste0(dfText, '.RData'))
    eval(substitute(save(DF, file = filename), list(DF=as.name(dfText))))
    if(ab$makePlots) eval(substitute(plotABS(DF), list(DF=as.name(dfText))))
    eval(substitute(printMinTimeABS(DF), list(DF=as.name(dfText))))
})
filename <- file.path(path, paste0('runPartitions.R'))
cat(codeToText(preCode), file=filename)
cat(codeToText(partitionsCode), file=filename, append=TRUE)

load('dfpartitionsN64.RData')
niter <- 200000
dfPartitionsN64$timePer10k <- dfPartitionsN64$timing *10000/niter
dfPartitionsN64$essPer10k  <- dfPartitionsN64$ess    *10000/niter * 2
dfPartitionsN64$Efficiency <- dfPartitionsN64$essPer10k / dfPartitionsN64$timePer10k
print(max(abs(dfPartitionsN64$Efficiency - dfPartitionsN64$essPT*2)))
dfPartitionsN64$rho <- as.numeric(gsub('.*rho(.+)', '\\1', dfPartitionsN64$model))
dfPartitionsN64$mcmc <- gsub('-.+', '', dfPartitionsN64$blocking)
dfN64 <- printMinTimeABS(dfPartitionsN64, round=FALSE)
save(dfPartitionsN64, dfN64, file='dfpartitionsN64.RData')



## 'mixedRhos' Simualted Data example
## mixed, overlapping, rhos
## used in a figure of simulated results, and also a table
tagValues <- LETTERS[1:3]
for(tag in tagValues) {
    mixedRhosCode <- substitute({
        tag <- TAG
        switch(tag,   ## multiples of 10
               A = { Nvalues <- c(20, 30, 40, 50, 60, 70, 80) },
               B = { Nvalues <- c(90, 100, 110, 120, 130, 140, 150) },
               C = { Nvalues <- c(160, 170, 180, 190, 200) },
               )
        control$niter <- 100000
        abList <- list()
        for(N in Nvalues) {
            DFtag <- paste0('mixedRhosN', N)
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
            abList[[DFtag]] <- ab
        }
        dfmixedRhos <- createDFfromABlist(abList)
        filename <- file.path(path, paste0('dfmixedRhos', tag, '.RData'))
        save(dfmixedRhos, file = filename)
        ##if(control$makePlots) eval(substitute(plotABS(DF), list(DF=as.name(dfText))))
        ##eval(substitute(printMinTimeABS(DF), list(DF=as.name(dfText))))
    }, list(TAG=tag))
    filename <- file.path(path, paste0('runmixedRhos', tag, '.R'))
    cat(codeToText(preCode), file=filename)
    cat(codeToText(mixedRhosCode), file=filename, append=TRUE)
}
filename <- 'runmixedRhos.sh'
cat('#!/bin/bash\n\n', file=filename)
for(tag in tagValues) {
    cat(paste0('R CMD BATCH --vanilla runmixedRhos', tag, '.R\n'), file=filename, append=TRUE)
    cat('git add --all\n', file=filename, append=TRUE)
    cat(paste0('git commit -a -m\'ran runmixedRhos', tag, '.R\'\n'), file=filename, append=TRUE)
    cat('git push\n\n', file=filename, append=TRUE)
}
system(paste0('chmod 777 ', filename))

## combining the A, B, C, ...  dataframes from mixedRhos
rm(list=ls())
dfCombined <- data.frame()
tagValues <- LETTERS[1:3]
for(tag in tagValues) {
    load(paste0('dfmixedRhos', tag, '.RData'))
    dfCombined <- rbind(dfCombined, dfmixedRhos)
}
dfmixedRhos <- dfCombined
save(dfmixedRhos, file = 'dfmixedRhos.RData')

load('dfmixedRhos.RData')
niter <- 100000
dfmixedRhos$timePer10k <- dfmixedRhos$timing *10000/niter
dfmixedRhos$essPer10k  <- dfmixedRhos$ess    *10000/niter * 2
dfmixedRhos$Efficiency <- dfmixedRhos$essPer10k / dfmixedRhos$timePer10k
print(max(abs(dfmixedRhos$Efficiency - dfmixedRhos$essPT*2)))
dfmixedRhos$N <- as.numeric(gsub('mixedRhosN(.+)', '\\1', dfmixedRhos$model))
dfmixedRhos$mcmc <- gsub('-.+', '', dfmixedRhos$blocking)
dfmixedRhos$model <- gsub('mixedRhos', '', dfmixedRhos$model)
dfmixedRhos$model <- gsub('N([23456789]0)$', 'N0\\1', dfmixedRhos$model)
dfMix <- printMinTimeABS(dfmixedRhos, round=FALSE)
save(dfmixedRhos, dfMix, file='dfmixedRhos.RData')




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

load('dfSSM.RData')
niter <- 400000
dfSSM$timePer10k <- dfSSM$timing *10000/niter
dfSSM$essPer10k  <- dfSSM$ess    *10000/niter * 2
dfSSM$Efficiency <- dfSSM$essPer10k / dfSSM$timePer10k
print(max(abs(dfSSM$Efficiency - dfSSM$essPT*2)))
dfSSM$mcmc <- gsub('-.+', '', dfSSM$blocking)
dfS <- printMinTimeABS(dfSSM, round=FALSE)
save(dfSSM, dfS, file='dfSSM.RData')






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

load('dflittersGAMMA-UNIFprior.RData')
niter <- 400000
dflitters$timePer10k <- dflitters$timing *10000/niter
dflitters$essPer10k  <- dflitters$ess    *10000/niter * 2
dflitters$Efficiency <- dflitters$essPer10k / dflitters$timePer10k
print(max(abs(dflitters$Efficiency - dflitters$essPT*2)))
dflitters$mcmc <- gsub('-.+', '', dflitters$blocking)
dfLit <- printMinTimeABS(dflitters, round=FALSE)
save(dflitters, dfLit, file='dflittersGAMMA-UNIFprior.RData')





## spatial
spatialCode <- quote({
    control$niter <- 200000
    control$saveSamples <- TRUE   ##### changed #####
    runList <- list('all', 'default', 'auto')
    abspatial <- autoBlock(code=code_spatial, constants=constants_spatial, data=data_spatial, inits=inits_spatial, control=control)
    abspatial$run(runList)
    abList <- list(spatial=abspatial)
    dfspatial <- createDFfromABlist(abList)
    filename <- file.path(path, 'dfspatialWithSamples.RData')   ##### changed #####
    abspatial$abModel <- NULL
    abspatial$Cmcmcs <- list()
    save(dfspatial, abspatial, file = filename)   ##### changed #####
    if(abspatial$makePlots) plotABS(dfspatial, xlimToMin=FALSE)
    if(abspatial$makePlots) plotABS(dfspatial, xlimToMin=TRUE)
    ##printMinTimeABS(dfspatial)
})
filename <- file.path(path, 'runspatial.R')
cat(codeToText(preCode), file=filename)
cat(codeToText(spatialCode), file=filename, append=TRUE)

load('dfspatialWithSamples.RData')
niter <- 200000
dfspatial$timePer10k <- dfspatial$timing *10000/niter
dfspatial$essPer10k  <- dfspatial$ess    *10000/niter * 2
dfspatial$Efficiency <- dfspatial$essPer10k / dfspatial$timePer10k
print(max(abs(dfspatial$Efficiency - dfspatial$essPT*2)))
dfspatial$mcmc <- gsub('-.+', '', dfspatial$blocking)
dfSpat <- printMinTimeABS(dfspatial, round=FALSE)
save(dfspatial, dfSpat, abspatial, file='dfspatialWithSamples.RData')
## now, make a dataframe of samples, sorted better
keepInd <- (dim(abspatial$samples[[1]])[1]/2+1):dim(abspatial$samples[[1]])[1]
nKeep <- length(keepInd)
params <- c('rho', 'sigma', 'mu', 'g[66]')
mcmcs <- c('all', 'default', 'auto0', 'autoMax')
mcmcInds <- c(1, 2, 3, 5)  ### CAREFUL ### need to set these manually to get them right
mcmcCol <- paramCol <- character();   sampCol <- numeric()
for(i in seq_along(mcmcs)) {
    mcmc <- mcmcs[i]
    mcmcInd <- mcmcInds[i]
    for(param in params) {
        mcmcCol  <- c(mcmcCol,  rep(mcmc,  nKeep))
        paramCol <- c(paramCol, rep(param, nKeep))
        sampCol  <- c(sampCol,  abspatial$samples[[mcmcInd]][keepInd, param])
    }
}
dfSpatSamples <- data.frame(mcmc=mcmcCol, param=paramCol, samp=sampCol)
save(dfspatial, dfSpat, abspatial, dfSpatSamples, file='dfspatialWithSamples.RData')






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



