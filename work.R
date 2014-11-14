preCode <- Quote({
    path <- '~/GitHub/autoBlock'
    control <- list(setSeed0 = TRUE)
    source(file.path(path, 'autoBlock_utils.R'))
})
eval(preCode)
preCode[[length(preCode)+1]] <- quote(control$makePlots <- FALSE)


## assesses the adapted scale, acceptance rates, ESS, and timing
## achieved by scalar/block samplers of various sizes, and underlying
## univariate or multivariate distributions
tagValues <- c(LETTERS[13:21])
for(tag in tagValues) {
    blockTestingCode <- substitute({
        tag <- TAG
        switch(tag,
               A = { dist <- 'uni';   Nvalues <- c(2, 3, 4, 5, 10, 20, 30, 40, 50) }, # 8 mins
               B = { dist <- 'uni';   Nvalues <- c(100, 150, 200, 250, 300) }, # 49 mins
               C = { dist <- 'uni';   Nvalues <- c(350, 400, 450, 500) }, # 105 mins
               D = { dist <- 'uni';   Nvalues <- c(600) }, # 34 mins
               E = { dist <- 'uni';   Nvalues <- c(700) }, # 47 mins
               F = { dist <- 'uni';   Nvalues <- c(800) }, # 62 mins
               G = { dist <- 'uni';   Nvalues <- c(900) }, # 72 mins
               H = { dist <- 'uni';   Nvalues <- c(1000) }, # 86 mins
               I = { dist <- 'multi'; Nvalues <- c(2, 3, 4, 5, 10, 20, 30, 40, 50) }, # 10 mins
               J = { dist <- 'multi'; Nvalues <- c(100, 150) }, # 51 mins
               K = { dist <- 'multi'; Nvalues <- c(200, 250) }, # 4 hours
               L = { dist <- 'multi'; Nvalues <- c(300) }, # ??
               ## different after here!  start here
               M = { dist <- 'multi'; Nvalues <- c(350) }, #
               N = { dist <- 'multi'; Nvalues <- c(400) }, #
               O = { dist <- 'multi'; Nvalues <- c(450) }, #
               P = { dist <- 'multi'; Nvalues <- c(500) }, #
               Q = { dist <- 'multi'; Nvalues <- c(600) }, #
               R = { dist <- 'multi'; Nvalues <- c(700) }, #
               S = { dist <- 'multi'; Nvalues <- c(800) }, #
               T = { dist <- 'multi'; Nvalues <- c(900) }, #
               U = { dist <- 'multi'; Nvalues <- c(1000) } # 
               )
        niter <- 400000
        keepInd <- (niter/2+1):niter
        ##optimalRates <- c(0.44, 0.35, 0.32, 0.25, 0.234)
        dfblockTesting <- data.frame()
        for(N in Nvalues) {
            cat(paste0('\nN = ', N, '\n'))
            cat(paste0('\ndist = ', dist, '\n'))
            candc <- if(dist == 'uni') createCodeAndConstants(N) else createCodeAndConstants(N, list(1:N), 0)
            code <- candc$code
            constants <- candc$constants
            data <- list()
            inits <- list(x = rep(0,N))
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
                ##aRateHistory <- sampler1$acceptanceRateHistory
                ##acceptRate[i] <- aRateHistory[length(aRateHistory)]
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
                ##adaptedScale = adaptedScale,
                ##adaptedPropSD = adaptedPropSD,
                derivedScale = c(adaptedScale[1], adaptedScale[2:3] * adaptedPropSD[2:3]),
                ##acceptRate = acceptRate,
                ##optRate = c(0.44, rep(optimalRates[if(N>5) 5 else N], 2)),
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
tagValues <- LETTERS[1:12]
for(tag in tagValues) {
    load(paste0('dfblockTesting', tag, '.RData'))
    dfCombined <- rbind(dfCombined, dfblockTesting)
}
dfblockTesting <- dfCombined
save(dfblockTesting, file = 'dfblockTesting.RData')

## make a plot of timing from blockTesting
rm(list=ls())
load('dfblockTesting.RData')
qplot(data=dfblockTesting, x=N, y=timePer10kN, color=blocking, geom='line', facet=dist)

## interesting relationships I noticed from blockTesting
y <- sqrt(df$N) * df$adaptedScale  ###### Interesting ! ! !
y <- df$meanESS * df$N



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



