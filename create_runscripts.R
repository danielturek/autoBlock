

codeToText <- function(code) {
    a <- deparse(code, width.cutoff=500L)
    a <- a[c(-1, -length(a))]
    a <- sub('^    ', '', a)
    a[length(a) + 1] <- ''
    a[length(a) + 1] <- ''
    a <- paste0(a, collapse='\n')
    return(a)
}

makeRunScript <- function(modelName, niter = 200000) {
    scriptCode <- substitute(
        {
            source('~/GitHub/autoBlock/autoBlock_utils.R')
            load(MODELFILE)
            control <- list(niter = NITER)
            ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
            ab$run(runList)
            abList <- list(ab)
            names(abList) <- MODELNAME
            DF <- createDFfromABlist(abList, NITER)
            DFSUMMARY <- printMinTimeABS(DF, round=FALSE)
            save(DF, DFSUMMARY, file = RESULTSFILE)            
        },
        list(MODELNAME = modelName,
             NITER     = niter,
             DF        = as.name(paste0('df', modelName)),
             DFSUMMARY = as.name(paste0('df', modelName, '_summary')),
             MODELFILE   = paste0('~/GitHub/autoBlock/modelfiles/model_', modelName, '.RData'),
             RESULTSFILE = paste0('~/GitHub/autoBlock/results/results_',  modelName, '.RData')
             )
    )
    filename <- paste0('~/GitHub/autoBlock/runscripts/run_', modelName, '.R')
    cat(codeToText(scriptCode), file = filename)
}


makeRunScript('litters')
makeRunScript('ice')
makeRunScript('SSMindependent')
makeRunScript('SSMcorrelated')
makeRunScript('spatial')
makeRunScript('mhp')
#makeRunScript('redblue')



## shows how bad sampling efficiency can be, as correlation increases
## output used in Figure 'dfsampEff'
samplingEfficiencyCode <- quote({
    library(nimble)
    library(coda)
    source('~/GitHub/autoBlock/autoBlock_utils.R')
    kValues <- 0:3
    Nvalues <- c(2, 4, 8)
    niter <- 200000
    keepInd <- (niter/2+1):niter
    dfsamplingEfficiency <- data.frame()
    for(expDecay in c(FALSE, TRUE)) {
        for(k in kValues) {
            for(N in Nvalues) {
                rho <- 1 - (1-0.8)^k
                candc <- createCodeAndConstants(N, list(1:N), rho, expDecay=expDecay)
                code <- candc$code
                constants <- candc$constants
                data <- list()
                inits <- list(x = rep(0,N))
                Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
                nodeNames <- Rmodel$expandNodeNames('x', returnScalarComponents = TRUE)
                spec <- configureMCMC(Rmodel, nodes = NULL)
                for(node in nodeNames) spec$addSampler('RW', list(targetNode=node), print=FALSE)
                Rmcmc <- buildMCMC(spec)
                compiledList <- compileNimble(list(Rmodel, Rmcmc))
                Cmodel <- compiledList[[1]]; Cmcmc <- compiledList[[2]]
                Cmodel$setInits(inits)
                set.seed(0)
                timing <- as.numeric(system.time(Cmcmc$run(niter))[1])
                timePer10kN <- timing / (niter/10000)
                samples <- as.matrix(Cmcmc$mvSamples)
                samples <- samples[keepInd, , drop = FALSE]
                ess <- apply(samples, 2, effectiveSize)
                meanESS <- mean(ess)
                essPerN <- meanESS / length(keepInd)
                samples <- NULL; Cmcmcs <- NA; gc()
                thisDF <- data.frame(expDecay=expDecay, k=k, rho=rho, N=N, timePer10kN=timePer10kN, essPerN=essPerN)
                dfsamplingEfficiency <- rbind(dfsamplingEfficiency, thisDF)
                save(dfsamplingEfficiency, file = '~/GitHub/autoBlock/results/results_samplingEfficiency.RData')
            }
        }
    }
})
filename <- '~/GitHub/autoBlock/runscripts/run_samplingEfficiency.R'
cat(codeToText(samplingEfficiencyCode), file = filename)






## assesses the adapted scale, acceptance rates, ESS, and timing
## achieved by scalar/block samplers of various sizes, and underlying
## univariate or multivariate distributions
## used in Figure: 'computationalRequirement'
computationalRequirementCode <- quote({
    library(nimble)
    source('~/GitHub/autoBlock/autoBlock_utils.R')
    niter <- 50000
    keepInd <- (niter/2+1):niter
    dfcomputationalRequirement <- data.frame()
    Nvalues <- c(2, 3, 5)
    for(dist in c('uni', 'multi', 'gamma')) {
        for(N in Nvalues) {
            if(dist == 'uni')   candc <- createCodeAndConstants(N)
            if(dist == 'multi') candc <- createCodeAndConstants(N, list(1:N), 0)
            if(dist == 'gamma') candc <- createCodeAndConstants(N, gammaScalars = TRUE)
            code <- candc$code
            constants <- candc$constants
            data <- list()
            inits <- list(x = rep(1,N))
            Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
            nodeNames <- Rmodel$expandNodeNames('x', returnScalarComponents = TRUE)
            specList <- list()  # ordering: scalar, blockNoAdapt, blockAdapt
            for(i in 1:3) specList[[i]] <- configureMCMC(Rmodel, nodes = NULL)
            for(node in nodeNames) specList[[1]]$addSampler('RW', list(targetNode=node), print=FALSE)
            specList[[2]]$addSampler('RW_block', list(targetNodes=nodeNames, adaptScaleOnly=TRUE), print=FALSE)
            specList[[3]]$addSampler('RW_block', list(targetNodes=nodeNames), print=FALSE)
            toCompileList <- list(Rmodel)
            for(i in 1:3) toCompileList[[i+1]] <- buildMCMC(specList[[i]])
            compiledList <- compileNimble(toCompileList)
            Cmodel <- compiledList[[1]]
            Cmcmcs <- compiledList[2:4]  # ordering: scalar, blockNoAdapt, blockAdapt
            timePer10kN <- numeric(0)
            for(i in 1:3) {
                Cmodel$setInits(inits)
                set.seed(0)
                timing <- as.numeric(system.time(Cmcmcs[[i]]$run(niter))[1])
                timePer10kN[i] <- timing / (niter/10000)
            }
            thisDF <- data.frame(
                N = rep(N, 3),
                dist = rep(dist, 3),
                blocking = c('scalar', 'blockNoAdapt', 'blockAdapt'),
                timePer10kN = timePer10kN
            )
            dfcomputationalRequirement <- rbind(dfcomputationalRequirement, thisDF)
            save(dfcomputationalRequirement, file = '~/GitHub/autoBlock/results/results_computationalRequirement.RData')
            cat('\n'); print(dfcomputationalRequirement)
        }
    }
})
filename <- '~/GitHub/autoBlock/runscripts/run_computationalRequirement.R'
cat(codeToText(computationalRequirementCode), file = filename)

## LEGACY
## (from running 'tag's of blockTesting code, then combining them)
## filename <- 'runblockTesting.sh'
## cat('#!/bin/bash\n\n', file=filename)
## for(tag in tagValues) {
##     cat(paste0('R CMD BATCH --vanilla runblockTesting', tag, '.R\n'), file=filename, append=TRUE)
##     cat('git add --all\n', file=filename, append=TRUE)
##     cat(paste0('git commit -a -m\'ran runBlockTesting', tag, '.R\'\n'), file=filename, append=TRUE)
##     cat('git push\n\n', file=filename, append=TRUE)
## }
## system(paste0('chmod 777 ', filename))
## ## combining the A, B, C, ...  dataframes from blockTesting
## rm(list=ls())
## dfCombined <- data.frame()
## tagValues <- LETTERS[1:13]
## for(tag in tagValues) {
##     load(paste0('dfblockTesting', tag, '.RData'))
##     dfCombined <- rbind(dfCombined, dfblockTesting)
## }
## load('dfblockTesting.RData')
## dfblockTestingGamma <- dfCombined
## save(dfblockTesting, dfblockTestingUni, dfblockTestingMulti, dfblockTestingGamma,
##      file = 'dfblockTesting.RData')







## 'varyingBlksFixedCorr' Simulated Data example
## N = 2^k, constant values of rho
## Used in a Figure of Simulated Data results, and probably a table
varyingBlksFixedCorrCode <- quote({
    library(nimble)
    source('~/GitHub/autoBlock/autoBlock_utils.R')
    k <- 6
    N <- 2^k
    rhoVector <- c(0.2, 0.5, 0.8)
    niter <- 200000
    control <- list(niter = niter)
    runList <- list('all', 'auto')
    abList <- list()
    for(rho in rhoVector) {
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
        abList[[paste0('varyingBlksFixedCorr', rho)]] <- ab
    }
    dfVaryingBlksFixedCorr <- createDFfromABlist(abList, niter)
    dfVaryingBlksFixedCorr$rho <- as.numeric(gsub('.*Corr(.+)', '\\1', dfVaryingBlksFixedCorr$model))
    dfVaryingBlksFixedCorr_summary <- printMinTimeABS(dfVaryingBlksFixedCorr, round=FALSE)
    save(dfVaryingBlksFixedCorr, dfVaryingBlksFixedCorr_summary, file = '~/GitHub/autoBlock/results/results_varyingBlksFixedCorr.RData')
})
filename <- '~/GitHub/autoBlock/runscripts/run_varyingBlksFixedCorr.R'
cat(codeToText(varyingBlksFixedCorrCode), file = filename)





## 'fixedBlksVaryingCorr' Simualted Data example
## mixed, overlapping, rhos
## used in a figure of simulated results, and also a table
fixedBlksVaryingCorrCode <- quote({
    library(nimble)
    source('~/GitHub/autoBlock/autoBlock_utils.R')
    Nvalues <- c(20, 50, 100)
    niter <- 200000
    control <- list(niter = niter)
    runList <- list('all', 'auto')
    abList <- list()
    for(N in Nvalues) {
        blockSize <- N/10
        numberOfBlocks <- 9
        indList <- lapply(((1:numberOfBlocks)-1)*blockSize, function(x) x+(1:blockSize))
        rhoVector <- seq(from=0.9, to=0.1, by=-0.1)
        codeAndConstants <- createCodeAndConstants(N, indList, rhoVector)
        code <- codeAndConstants$code
        constants <- codeAndConstants$constants
        data <- list()
        inits <- list(x=rep(0,N))
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[[paste0('fixedBlksVaryingCorrN', N)]] <- ab
    }
    dfFixedBlksVaryingCorr <- createDFfromABlist(abList, niter)
    dfFixedBlksVaryingCorr$N <- as.numeric(gsub('.*N(.+)', '\\1', dfFixedBlksVaryingCorr$model))
    dffixedBlksVaryingCorr$model <- gsub('fixedBlksVaryingCorr', '', dffixedBlksVaryingCorr$model)
    dffixedBlksVaryingCorr$model <- gsub('N([23456789]0)$', 'N0\\1', dffixedBlksVaryingCorr$model)
    dfFixedBlksVaryingCorr_summary <- printMinTimeABS(dfFixedBlksVaryingCorr, round=FALSE)
    save(dfFixedBlksVaryingCorr, dfFixedBlksVaryingCorr_summary, file = '~/GitHub/autoBlock/results/results_fixedBlksVaryingCorr.RData')
})
filename <- '~/GitHub/autoBlock/runscripts/run_fixedBlksVaryingCorr.R'
cat(codeToText(fixedBlksVaryingCorrCode), file = filename)



## LEGACY
## filename <- 'runfixedBlksVaryingCorr.sh'
## cat('#!/bin/bash\n\n', file=filename)
## for(tag in tagValues) {
##     cat(paste0('R CMD BATCH --vanilla runfixedBlksVaryingCorr', tag, '.R\n'), file=filename, append=TRUE)
##     cat('git add --all\n', file=filename, append=TRUE)
##     cat(paste0('git commit -a -m\'ran runfixedBlksVaryingCorr', tag, '.R\'\n'), file=filename, append=TRUE)
##     cat('git push\n\n', file=filename, append=TRUE)
## }
## system(paste0('chmod 777 ', filename))
## ## combining the A, B, C, ...  dataframes from fixedBlksVaryingCorr
## rm(list=ls())
## dfCombined <- data.frame()
## tagValues <- LETTERS[1:3]
## for(tag in tagValues) {
##     load(paste0('dffixedBlksVaryingCorr', tag, '.RData'))
##     dfCombined <- rbind(dfCombined, dffixedBlksVaryingCorr)
## }
## dffixedBlksVaryingCorr <- dfCombined
## save(dffixedBlksVaryingCorr, file = 'dffixedBlksVaryingCorr.RData')





















