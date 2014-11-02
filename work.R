preCode <- quote({
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
Nvalues <- c(30, 100)   ## multiples of 10
for(Nval in Nvalues) {
    mixedRhosCode <- substitute({
        N <- NNN
        tag <- paste0('mixedRhosN', N)
        blockSize <- N/10
        numberOfBlocks <- 9
        indList <- lapply(((1:numberOfBlocks)-1)*blockSize, function(x) x+(1:blockSize))
        rhoVector <- seq(from=0.9, to=0.1, by=-0.1)
        control$niter <- 400000
        runList <- list('all', 'auto')
        codeAndConstants <- createCodeAndConstants(N, indList, rhoVector)
        code <- codeAndConstants$code
        constants <- codeAndConstants$constants
        data <- list()
        inits <- list(x=rep(0,N))
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList <- list(mixedRhos=ab)
        dfText <- paste0('df', tag)
        eval(substitute(DF <- createDFfromABlist(abList), list(DF=as.name(dfText))))
        filename <- file.path(path, paste0(dfText, '.RData'))
        eval(substitute(save(DF, file = filename), list(DF=as.name(dfText))))
        if(control$makePlots) eval(substitute(plotABS(DF), list(DF=as.name(dfText))))
        eval(substitute(printMinTimeABS(DF), list(DF=as.name(dfText))))
    }, list(NNN=Nval))
    filename <- file.path(path, paste0('runMixedRhosN', Nval, '.R'))
    cat(codeToText(preCode), file=filename)
    cat(codeToText(mixedRhosCode), file=filename, append=TRUE)
}






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








