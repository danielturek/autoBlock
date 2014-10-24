path <- '/GitHub/autoBlock'
source(file.path(path, 'autoBlock_utils.R'))
source(file.path(path, 'autoBlock_models.R'))
control <- list(
    sparsifyCov = TRUE,
    sparseCovThreshold = 0.001,
    cutree_method = 'custom',
    cutree_h = 0.99,
    cutree_maxGroupSize = 5,
    cutree_maxHeightRelativeFromBase = 0.2,
    cutree_maxGroupSizeRelHeightOveride = 0.1,
    verbose = TRUE
    )


## partitions of N = 2^k
if(FALSE) {
    k <- 6  ### next to do: k=6
    rho <- 0.8
    runList <- list(givenCov = quote({ spec <- MCMCspec(Rmodel, nodes = NULL)
                                       spec$addSampler('RW_block', control=list(targetNodes='x', adaptive=TRUE, adaptScaleOnly=TRUE, propCov = Sigma), print=FALSE); spec }),
                    'all',
                    'auto')
    N <- 2^k
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
        blockLengths <- c(1, 2^(0:(k-1)));    if(k==0) stop('k=0 doesnt work here')
        indList <- list(); cur <- 1
        for(len in blockLengths) { indList <- c(indList, list(cur:(cur+len-1))); cur <- cur+len }
        Sigma <- createCov(N, indList=indList, rho=rho)
        constants <- list(N=N, mu=rep(0,N), Sigma=Sigma)
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[['blockSzMixed']] <- ab
    }
    df <- createDFfromABlist(abList, rho=rho)
    plotABS(df)
}
if(FALSE) {
    load('dfs.RData')
    dfs[[k]] <- df
    save('dfs', file='dfs.RData')
}
#plotABS(dfs[[5]])


## testing of SSMs
if(FALSE) {
    runListMUB <- list('all', 'auto', blockMUB = list(c('mu','b')))
    runListAB  <- list('all', 'auto', blockAB   = list(c('a', 'b')))
    abSSMmub$run(runListMUB)
    abSSMab$run(runListAB)
    abList <- list(independent=abSSMmub, correlated=abSSMab)
    dfSSM <- createDFfromABlist(abList)
    plotABS(dfSSM, xlimToMin=TRUE)
    save('dfSSM', file='dfSSM.RData')
}
    


## testing of litters
if(FALSE) {
    runList <- list('all', 'auto', blockAB = list(c('a[1]','b[1]'), c('a[2]','b[2]')))
    ablitters$run(runList)
    abList <- list(litters=ablitters)
    dflitters <- createDFfromABlist(abList)
    plotABS(dflitters, xlimToMin=TRUE)
    save('dflitters', file='dflitters.RData')
}


load('dfs.RData')
dfs
plotABS(dfs[[3]])
dfs[[2]]$naming
df <- dfs[[2]]
df[df$model=='blockSz2' , 'essPT']


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








