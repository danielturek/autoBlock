
resultsDirectoryName <- 'results_hclust_single'
fast <- TRUE

if(fast) {
    niter_examples <- 5000
    niter_compReq  <- 1000
} else {
    niter_examples <- 200000
    niter_compReq  <- 50000
}

codeToText <- function(code) {
    a <- deparse(code, width.cutoff=500L)
    a <- a[c(-1, -length(a))]
    a <- sub('^    ', '', a)
    a[length(a) + 1] <- ''
    a[length(a) + 1] <- ''
    a <- paste0(a, collapse='\n')
    return(a)
}

makeRunScript <- function(modelName) {
    scriptCode <- substitute(
        {
            source('autoBlock.R')
            load(file.path('data', MODELFILE))
            OUT <- autoBlock(code, constants, data, inits, NITER, runList)$summary
            save(OUT, file = file.path(RESULTSDIRECTORY, RESULTSFILE))
        },
        list(OUT              = as.name(paste0('df', modelName)),
             MODELFILE        = paste0('model_',   modelName, '.RData'),
             NITER            = niter_examples,
             RESULTSDIRECTORY = resultsDirectoryName,
             RESULTSFILE      = paste0('results_', modelName, '.RData')
             )
    )
    filename <- paste0('~/GitHub/legacy/autoBlock/code/run_', modelName, '.R')
    cat(codeToText(scriptCode), file = filename)
}


makeRunScript('litters')
makeRunScript('ice')
makeRunScript('SSMindependent')
makeRunScript('SSMcorrelated')
makeRunScript('spatial')
makeRunScript('mhp')
##makeRunScript('redblue')
makeRunScript('test')





## shows how bad sampling efficiency can be, as correlation increases
## output used in Figure 'dfsampEff'
samplingEfficiencyCode <- substitute(
    {
        library(nimble)
        library(coda)
        source('autoBlock.R')
        kValues <- 0:3
        Nvalues <- c(2, 4, 8, 16)
        niter <- NITER
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
                    Cmodel <- compiledList[[1]]
                    Cmcmc <- compiledList[[2]]
                    Cmodel$setInits(inits)
                    set.seed(0)
                    timing <- as.numeric(system.time(Cmcmc$run(niter))[1])
                    timePer10kN <- timing / (niter/10000)
                    samples <- as.matrix(Cmcmc$mvSamples)
                    samples <- samples[keepInd, , drop = FALSE]
                    ess <- apply(samples, 2, effectiveSize)
                    meanESS <- mean(ess)
                    essPerN <- meanESS / length(keepInd)
                    thisDF <- data.frame(expDecay=expDecay, k=k, rho=rho, N=N, timePer10kN=timePer10kN, essPerN=essPerN)
                    dfsamplingEfficiency <- rbind(dfsamplingEfficiency, thisDF)
                    save(dfsamplingEfficiency, file = file.path(RESULTSDIRECTORY, 'results_samplingEfficiency.RData'))
                }
            }
        }
    },
    list(NITER            = niter_examples,
         RESULTSDIRECTORY = resultsDirectoryName
         )
)
filename <- '~/GitHub/legacy/autoBlock/code/run_samplingEfficiency.R'
cat(codeToText(samplingEfficiencyCode), file = filename)






## assesses the adapted scale, acceptance rates, ESS, and timing
## achieved by scalar/block samplers of various sizes, and underlying
## univariate or multivariate distributions
## used in Figure: 'computationalRequirement'
computationalRequirementCode <- substitute(
    {
        library(nimble)
        source('autoBlock.R')
        niter <- NITER
        keepInd <- (niter/2+1):niter
        dfcomputationalRequirement <- data.frame()
        Nvalues <- c(2, 3)
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
                save(dfcomputationalRequirement, file = file.path(RESULTSDIRECTORY, 'results_computationalRequirement.RData'))
                cat('\n'); print(dfcomputationalRequirement)
            }
        }
    },
    list(
        NITER            = niter_compReq,
        RESULTSDIRECTORY = resultsDirectoryName
    )
)
filename <- '~/GitHub/legacy/autoBlock/code/run_computationalRequirement.R'
cat(codeToText(computationalRequirementCode), file = filename)




## 'varyingBlksFixedCorr' Simulated Data example
## N = 2^k, constant values of rho
## Used in a Figure of Simulated Data results, and probably a table
varyingBlksFixedCorrCode <- substitute(
    {
        library(nimble)
        source('autoBlock.R')
        k <- 6
        N <- 2^k
        rhoVector <- c(0.2, 0.5, 0.8)
        niter <- NITER
        runList <- list('all', 'auto')
        dfVaryingBlksFixedCorr <- NULL
        for(rho in rhoVector) {
            blockLengths <- c(1, 2^(0:(k-1)))
            indList <- list(); cur <- 1
            for(len in blockLengths) { indList <- c(indList, list(cur:(cur+len-1))); cur <- cur+len }
            data <- list()
            inits <- list(x=rep(0,N))
            codeAndConstants <- createCodeAndConstants(N, indList, rep(rho,length(indList)))
            code <- codeAndConstants$code
            constants <- codeAndConstants$constants
            dfTEMP <- autoBlock(code=code, constants=constants, data=data, inits=inits, niter=niter, run=runList)$summary
            dfTEMP <- cbind(data.frame(rho=rho), dfTEMP)
            if(is.null(dfVaryingBlksFixedCorr)) dfVaryingBlksFixedCorr <- dfTEMP
            else dfVaryingBlksFixedCorr <- rbind(dfVaryingBlksFixedCorr, dfTEMP)
        }
        save(dfVaryingBlksFixedCorr, file = file.path(RESULTSDIRECTORY, 'results_varyingBlksFixedCorr.RData'))
    },
    list(
        NITER            = niter_compReq,
        RESULTSDIRECTORY = resultsDirectoryName
    )
)
filename <- '~/GitHub/legacy/autoBlock/code/run_varyingBlksFixedCorr.R'
cat(codeToText(varyingBlksFixedCorrCode), file = filename)





## 'fixedBlksVaryingCorr' Simualted Data example
## mixed, overlapping, rhos
## used in a figure of simulated results, and also a table
fixedBlksVaryingCorrCode <- substitute(
    {
        library(nimble)
        source('autoBlock.R')
        Nvalues <- c(20, 30, 50)
        niter <- NITER
        runList <- list('all', 'auto')
        dfFixedBlksVaryingCorr <- NULL
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
            dfTEMP <- autoBlock(code=code, constants=constants, data=data, inits=inits, niter=niter, run=runList)$summary
            dfTEMP <- cbind(data.frame(N=N), data.frame(model=paste0('N',N)), dfTEMP)
            if(is.null(dfFixedBlksVaryingCorr)) dfFixedBlksVaryingCorr <- dfTEMP
            else dfFixedBlksVaryingCorr <- rbind(dfFixedBlksVaryingCorr, dfTEMP)
        }
        save(dfFixedBlksVaryingCorr, file = file.path(RESULTSDIRECTORY, 'results_fixedBlksVaryingCorr.RData'))
    },
    list(
        NITER            = niter_compReq,
        RESULTSDIRECTORY = resultsDirectoryName
    )
)
filename <- '~/GitHub/legacy/autoBlock/code/run_fixedBlksVaryingCorr.R'
cat(codeToText(fixedBlksVaryingCorrCode), file = filename)



