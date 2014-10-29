path <- '~/GitHub/autoBlock'
#control <- list(cutree_heights = seq(0, 1, by=0.1))
source(file.path(path, 'autoBlock_utils.R'))





## state space models
if(FALSE) {
    runListMUB <- list('all', 'auto', blockMUB = list(c('mu','b')), 'default')
    runListAB  <- list('all', 'auto', blockAB  = list(c('a', 'b')), 'default')
    abSSMmub <- autoBlock(code=code_SSMmub, constants=constants_SSMmub, data=data_SSMmub, inits=inits_SSMmub, control=control)
    abSSMab <- autoBlock(code=code_SSMab, constants=constants_SSMab, data=data_SSMab, inits=inits_SSMab, control=control)
    abSSMmub$run(runListMUB)
    abSSMab$run(runListAB)
    abList <- list(independent=abSSMmub, correlated=abSSMab)
    dfSSM <- createDFfromABlist(abList)
    save('dfSSM', file='dfSSM.RData')
    plotABS(dfSSM, xlimToMin=TRUE)
}


## litters
if(FALSE) {
    runList <- list('all',
                    'auto',
                    blockAB = list(c('a[1]','b[1]'), c('a[2]','b[2]')),
                    crossLevel = quote({
                        spec <- MCMCspec(Rmodel, nodes=NULL)
                        spec$addSampler('crossLevel', list(topNodes = c('a[1]', 'b[1]')), print=FALSE)
                        spec$addSampler('crossLevel', list(topNodes = c('a[2]', 'b[2]')), print=FALSE)
                        spec
                    }),
                    'default')
    ablitters <- autoBlock(code=code_litters, constants=constants_litters, data=data_litters, inits=inits_litters, control=control)
    ablitters$run(runList)
    abList <- list(litters=ablitters)
    dflitters <- createDFfromABlist(abList)
    save('dflitters', file='dflitters.RData')
    plotABS(dflitters, xlimToMin=TRUE)
}



## partitions of N = 2^k
if(FALSE) {
    rho <- 0.7
    k <- 3
    N <- 2^k
    tag <- paste0('N', N, 'rho', rho)
    runList <- list(
        givenCov = quote({
            spec <- MCMCspec(Rmodel, nodes = NULL)
            spec$addSampler('RW_block', control=list(targetNodes='x', adaptive=TRUE, adaptScaleOnly=TRUE, propCov = Sigma), print=FALSE)
            spec }),
        'all', 'auto')
    blockSizes <- 2^(0:k)
    code <- quote({ x[1:N] ~ dmnorm(mu[1:N], cov = Sigma[1:N,1:N]) })
    data <- list()
    inits <- list(x=rep(0,N))
    abList <- list()
    for(blockSize in blockSizes) {
        numberOfBlocks <- N / blockSize
        listOfBlockIndexs <- lapply(((1:numberOfBlocks)-1)*blockSize, function(x) x+(1:blockSize))
        Sigma <- createCov(N, indList=listOfBlockIndexs, rho=rho)
        constants <- list(N=N, mu=rep(0,N), Sigma=Sigma)
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[[paste0('blockSz', blockSize)]] <- ab
    }
    if(k > 1) {
        blockLengths <- c(1, 2^(0:(k-1)))
        indList <- list(); cur <- 1
        for(len in blockLengths) { indList <- c(indList, list(cur:(cur+len-1))); cur <- cur+len }
        Sigma <- createCov(N, indList=indList, rho=rho)
        constants <- list(N=N, mu=rep(0,N), Sigma=Sigma)
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[['blockSzMixed']] <- ab
    }
    dfText <- paste0('df', tag)
    eval(substitute(DF <- createDFfromABlist(abList, rho=rho), list(DF=as.name(dfText))))
    eval(substitute(save(DF, file=FILE), list(DF=as.name(dfText), FILE=paste0(dfText,'.RData'))))
    eval(substitute(plotABS(DF), list(DF=as.name(dfText))))
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
    littersTraceplots('a[1]', samples)   # 50,000
    littersTraceplots('b[1]', samples)   # 5,000
    littersTraceplots('a[2]', samples)   # 100  (200 too high)
    littersTraceplots('b[2]', samples)   # 50 (100 to high)

}










## testing the new adaptation options; and a nice N=3 example
## N <- 3
## Sigma <- createCov(N, rho=0.8)
## code <- quote({ x[1:N] ~ dmnorm(mu[1:N], cov = Sigma[1:N,1:N]) })
## constants <- list(N=N, mu=rep(0,N), Sigma=Sigma)
## data <- list()
## inits <- list(x=rep(0,N))
## runList <- list(adaptTT = quote({ spec <- MCMCspec(Rmodel, nodes = NULL)
##                                   spec$addSampler('RW_block', control=list(targetNodes='x', adaptive=TRUE, adaptScaleOnly=TRUE))
##                                   spec}),
##                 adaptTF = quote({ spec <- MCMCspec(Rmodel, nodes = NULL)
##                                   spec$addSampler('RW_block', control=list(targetNodes='x', adaptive=TRUE, adaptScaleOnly=FALSE))
##                                   spec}),
##                 adaptFT = quote({ spec <- MCMCspec(Rmodel, nodes = NULL)
##                                   spec$addSampler('RW_block', control=list(targetNodes='x', adaptive=FALSE, adaptScaleOnly=TRUE))
##                                   spec}),
##                 adaptFF = quote({ spec <- MCMCspec(Rmodel, nodes = NULL)
##                                   spec$addSampler('RW_block', control=list(targetNodes='x', adaptive=FALSE, adaptScaleOnly=FALSE))
##                                   spec})
##                 )
## ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
## ab$run(runList)
## abList <- list(N3 = ab)
## df <- createDFfromABlist(abList)
## plotABS(df)




## some general code for testing run times
## library(nimble)
## N <- 200
## points <- c(10,50,100,150,200)
## niter <- 10000
## code <- quote({
##     x[1:N] ~ dmnorm(mu[1:N], cov = Sigma[1:N,1:N])
## })
## constants <- list(N=N, mu=rep(0,N), Sigma=diag(N))
## data <- list()
## inits <- list(x=rep(0,N))
## md <- nimbleModel(code=code, constants=constants, data=data, inits=inits, returnDef=TRUE)
## timing <- numeric(0)
## for(i in seq_along(points)) {
##     pt <- points[i]
##     cat(paste0(i, '\n'))
##     Rmodel <- md$newModel()
##     spec <- MCMCspec(Rmodel, nodes=NULL)
##     targetNodes <- paste0('x[1:', pt, ']')
##     spec$addSampler('RW_block', control=list(targetNodes=targetNodes))
##     spec$addMonitors('x', print=FALSE)
##     Rmcmc <- buildMCMC(spec)
##     Cmodel <- compileNimble(Rmodel)
##     Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
##     timing[i] <- system.time(Cmcmc(niter))[3]
## }
## plot(points, timing, pch=19, ylim=c(0, max(timing)))






