path <- '~/GitHub/autoBlock'
control <- list(
    cutree_heights = seq(0, 1, by=0.1),
    verbose = TRUE
    )
save(control, file=file.path(path, '.control.RData'))
source(file.path(path, 'autoBlock_utils.R'))
#load('results.RData')





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

    dfT <- df
    abListT <- abList
    load('results.RData')
    dim(df)
    dim(dfT)
    length(abList)
    length(abListT)
    df <- rbind(df, dfT)
    abList <- c(abList, abListT)
    save(df, abList, file='results.RData')
    plotABS(df)
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
    abSSMmub <- autoBlock(code=code_SSMmub, constants=constants_SSMmub, data=data_SSMmub, inits=inits_SSMmub, control=control)
    abSSMab <- autoBlock(code=code_SSMab, constants=constants_SSMab, data=data_SSMab, inits=inits_SSMab, control=control)
    abSSMmub$run(runListMUB)
    abSSMab$run(runListAB)
    abList <- list(independent=abSSMmub, correlated=abSSMab)
    dfSSM <- createDFfromABlist(abList)
    plotABS(dfSSM, xlimToMin=FALSE)
    save('dfSSM', file='dfSSM.RData')
}
    


## testing of litters
if(FALSE) {
    runList <- list('all', 'auto', blockAB = list(c('a[1]','b[1]'), c('a[2]','b[2]')))
    ablitters <- autoBlock(code=code_litters, constants=constants_litters, data=data_litters, inits=inits_litters, control=control)
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

bsamp <- list()
for(i in 1:8) bsamp[[i]] <- samples[[i]][50001:100000,]
names(bsamp) <- names(samples)


a1b1 <- c('a[1]', 'b[1]')
a2b2 <- c('a[2]', 'b[2]')
ab <- c('a[1]', 'a[2]', 'b[1]', 'b[2]')

for(i in 1:8) {
    print(names(bsamp)[i])
    print('mean')
    print(apply(bsamp[[i]][, ab], 2, mean))
    print('var')
    print(apply(bsamp[[i]][, ab], 2, var))
}



1


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






