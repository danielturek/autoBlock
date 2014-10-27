path <- '~/GitHub/autoBlock'
control <- list(
    cutree_method = 'custom',
    cutree_h = 0.9,
    cutree_maxGroupSize = 8,
    cutree_maxGroupSizeDynamic = TRUE,
    cutree_maxHeightRelativeFromBase = 0.3,
    verbose = TRUE
    )
save(control, file=file.path(path, '.control.RData'))
source(file.path(path, 'autoBlock_utils.R'))
load('results.RData')





## partitions of N = 2^k
rho <- 0.8
kValues <- 1:7
RscriptName <- 'gen_partitionsN.R'
shellScriptName <- 'gen_partitionsNwrapper.sh'
partitionsNcode <- substitute({
    library(R.utils)
    source('autoBlock_utils.R')
    load('.control.RData')
    control$verbose <- FALSE
    args <- commandArgs(trailingOnly=TRUE, asValue=TRUE)
    k <- as.numeric(args[1])
    rho <- as.numeric(args[2])
    runList <- list(
        givenCov = quote({
            spec <- MCMCspec(Rmodel, nodes = NULL)
            spec$addSampler('RW_block', control=list(targetNodes='x', adaptive=TRUE, adaptScaleOnly=TRUE, propCov = Sigma), print=FALSE)
            spec }),
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
        blockLengths <- c(1, 2^(0:(k-1)))
        indList <- list(); cur <- 1
        for(len in blockLengths) { indList <- c(indList, list(cur:(cur+len-1))); cur <- cur+len }
        Sigma <- createCov(N, indList=indList, rho=rho)
        constants <- list(N=N, mu=rep(0,N), Sigma=Sigma)
        ab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)
        ab$run(runList)
        abList[['blockSzMixed']] <- ab
    }
    df <- createDFfromABlist(abList, rho=rho)
    for(ab in abList) { ab$abModel <- NULL; ab$Cmcmcs <- list() }
    DF <- list()
    ABLIST <- list()
    if(file.exists('results.RData')) load('results.RData')
    DF[[k]] <- df
    ABLIST[[k]] <- abList
    save(DF, ABLIST, file='results.RData')
},
                              list(DF = as.name(paste0('dfNrho', rho)),
                                   ABLIST = as.name(paste0('abListNrho', rho))))
cat(formatForFile(partitionsNcode), file=RscriptName)
f <- file(shellScriptName)
writeLines(c(
    paste0('for k in ', paste(kValues, collapse=' ')),
    'do',
        paste0('    R CMD BATCH --no-save --no-restore --no-timing --args -k=$k -rho=', rho, ' ', RscriptName, ' out'),
    'done'
##    'rm out'
), f)
close(f)
system(paste0('chmod 755 ', shellScriptName))



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
    runList <- list(
        'all',
        'auto',
        blockAB = list(c('a[1]','b[1]'), c('a[2]','b[2]')),
        blockABP = list(c('a[1]','b[1]'), c('a[2]','b[2]'), 'p[1,1:16]', 'p[2,1:16]'),
        blockABPhalf = list(c('a[1]','b[1]'), c('a[2]','b[2]'), 'p[1,1:8]', 'p[1,9:16]', 'p[2,1:8]', 'p[2,9:16]')
        )
    ablitters$run(runList)
    abList <- list(litters=ablitters)
    dflitters <- createDFfromABlist(abList)
    plotABS(dflitters, xlimToMin=TRUE)
    save('dflitters', file='dflitters.RData')
}












Cmcmcs <- ablitters$Cmcmcs
samples <- ablitters$samples
empCov <- ablitters$empCov
empCovSparse <- ablitters$empCovSparse
empCovSparseThresh <- ablitters$empCovSparseThresh
grouping <- ablitters$grouping
length(ablitters$empCov)
length(samples)
names(samples)
dim(samples[[2]])
cov(samples[[2]]) - empCovSparseThresh[[3]]
length(empCov)
length(empCovSparse)
empCovSparseThresh[[3]]
ablitters$sparseCovThreshold
cov2cor(empCov[[3]])[1:4,1:4]
nfVar(Cmcmcs[[3]], 'samplerFunctions')$contentsList[[1]]$propCov
cov2cor(nfVar(Cmcmcs[[3]], 'samplerFunctions')$contentsList[[3]]$propCov)
cov(samples[[3]])[1:4,1:4]
cov2cor(cov(samples[[3]]))[1:4,1:4]
length(Cmcmcs)
grouping[[3]]


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




library(cluster)
Rmodel <- ablitters$abModel$newModel()
spec <- MCMCspec(Rmodel, onlyRW=TRUE)                      
spec$addSampler('RW_block', list(targetNodes=c('a[1]','b[1]')))
spec$addSampler('RW_block', list(targetNodes=c('a[2]','b[2]')))
spec$setSamplers(c(37,38,5:36))
spec <- MCMCspec(Rmodel, nodes=NULL)
spec$addSampler('RW_block', list(targetNodes=c('a[1]','b[1]')))
spec$addSampler('RW_block', list(targetNodes=c('a[2]','b[2]')))
spec$addSampler('RW_block', list(targetNodes=c('p[1, 1:16]')))
spec$addSampler('RW_block', list(targetNodes=c('p[2, 1:16]')))
spec$getSamplers()
spec$addMonitors(c('a','b','p'))
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
set.seed(0)
Cmcmc(100000)
samples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))
Cor <- cor(samples)
cov(samples)[1:4,1:4]
cor(samples)[1:4,1:4]
empCor <- cor(samples)
dd <- as.dist(1-abs(empCor))
hTree <- hclust(dd)
plot(hTree)
groups <- cutree_custom(hTree, maxHeight=0.9, maxGroupSize=10, maxHeightRelativeFromBase=0.3, maxGroupSizeRelHeightOveride=0.1)
groups
i <- 4
nfVar(Cmcmc, 'samplerFunctions')$contentsList[[i]]$d
nfVar(Cmcmc, 'samplerFunctions')$contentsList[[i]]$acceptanceRateHistory
nfVar(Cmcmc, 'samplerFunctions')$contentsList[[i]]$propCov
cov2cor(nfVar(Cmcmc, 'samplerFunctions')$contentsList[[i]]$propCov)

ag <- agnes(dd)
plot(ag)



runList <- list(
    'auto', 'all'

    )
ablitters$run(runList)
abList <- list(litters=ablitters)
df <- createDFfromABlist(abList)
plotABS(df, xlimToMin=TRUE)
for(i in 4:7) {
    print(ablitters$naming[[i]])
    print(ablitters$essPT[[i]][c('a[1]','a[2]','b[1]','b[2]')])
    print(cor(ablitters$samples[[i]])[1:4,1:4])
}



