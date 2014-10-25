

autoBlockModel <- setRefClass(
    Class = 'autoBlockModel',
    fields = list(
        code = 'ANY',
        constants = 'list',
        data = 'list',
        inits = 'list',
        md = 'ANY',
        Rmodel = 'ANY',
        scalarNodeVector = 'character',
        nodeGroupScalars = 'list',
        nodeGroupAllBlocked = 'list',
        nodeGroupStochNodes = 'list'
        ),
    methods = list(
        initialize = function(code, constants, data, inits) {
            library(nimble)
            code <<- code
            constants <<- if(missing(constants)) list() else constants
            data <<- if(missing(data)) list() else data
            inits <<- if(missing(inits)) list() else inits
            md <<- nimbleModel(code=code, constants=constants, returnDef=TRUE)
            Rmodel <<- md$newModel(data=data, inits=inits)
            scalarNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
            nodeGroupScalars <<- lapply(scalarNodeVector, function(x) x)
            nodeGroupAllBlocked <<- list(scalarNodeVector)
            stochNodeVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=FALSE)
            nodeGroupStochNodes <<- lapply(stochNodeVector, function(x) x)
        },
        createGroups = function(listOfBlocks = list()) {
            listOfBlocks <- lapply(listOfBlocks, function(blk) Rmodel$expandNodeNames(blk, returnScalarComponents=TRUE))
            nodes <- scalarNodeVector
            nodes <- setdiff(nodes, unlist(listOfBlocks))
            nodeList <- lapply(nodes, function(x) x)
            for(ng in listOfBlocks) nodeList[[length(nodeList)+1]] <- ng
            return(nodeList)
        },
        newModel = function() {
            newRmodel <- md$newModel(data=data, inits=inits)
            return(newRmodel)
        }
        )
    )

createCov <- function(N, indList=list(1:N), rho=0.8, indList2=list(), rho2=0.3) {
    Sigma <- diag(N)
    for(gp in indList)  { for(i1 in gp) for(i2 in gp) Sigma[i1,i2] <- Sigma[i2,i1] <- rho  }
    for(gp in indList2) { for(i1 in gp) for(i2 in gp) Sigma[i1,i2] <- Sigma[i2,i1] <- rho2 }
    diag(Sigma) <- 1
    Sigma
}


autoBlockParamDefaults <- function() {
    list(
        adaptInterval = 200,
        adaptIntervalBlock = 200,
        cutree_h = 0.9,
        cutree_maxGroupSize = 10,
        cutree_maxGroupSizeDynamic = FALSE,
##        cutree_maxGroupSizeRelHeightOveride = 0.2,
        cutree_maxHeightRelativeFromBase = 0.3,
        cutree_method = 'custom',
        debug = FALSE,
        niter = 100000,
        setSeed0 = TRUE,
        sparsifyCov = FALSE,
        sparseCovThreshold = -1,
        spcov_lambda = 0.06,
        spcov_step = 100,
        useConjugacy = FALSE,
        verbose = TRUE
        )
}


autoBlock <- setRefClass(

    Class = 'autoBlock',

    fields = list(
        
        abModel = 'ANY',

        ## overall control
        adaptInterval = 'numeric',
        adaptIntervalBlock = 'numeric',
        cutree_h = 'numeric',
        cutree_maxGroupSize = 'numeric',
        cutree_maxGroupSizeDynamic = 'logical',
##        cutree_maxGroupSizeRelHeightOveride = 'numeric',
        cutree_maxHeightRelativeFromBase = 'numeric',
        cutree_method = 'character',
        debug = 'logical',
        niter = 'numeric',
        setSeed0 = 'logical',
        sparsifyCov = 'logical',
        sparseSuccess = 'logical',   ## special one
        sparseCovThreshold = 'numeric',
        spcov_lambda = 'numeric',
        spcov_step = 'numeric',
        useConjugacy = 'logical',
        verbose = 'logical',

        it = 'numeric',
        
        ## persistant LISTS of historical data
        empCov = 'list',
        empCovSparse = 'list',
        empCovSparseThresh = 'list',
        empCor = 'list',
        distMatrix = 'list',
        hTree = 'list',
        naming = 'list',
        grouping = 'list',
        groupSizes = 'list',
        Cmcmcs = 'list',
        timing = 'list',
        samples = 'list',
        ess = 'list',
        essPT = 'list'
        ),

    methods = list(

        initialize = function(code, constants=list(), data=list(), inits=list(), control=list()) {
            library(lattice)
            library(coda)
            library(dynamicTreeCut)
            library(nimble)
            library(spcov)
            abModel <<- autoBlockModel(code=code, constants=constants, data=data, inits=inits)
            defaultsList <- autoBlockParamDefaults()
            for(i in seq_along(defaultsList)) if(is.null(control[[names(defaultsList)[i]]])) control[[names(defaultsList)[i]]] <- defaultsList[[i]]
            for(i in seq_along(control)) eval(substitute(verbose <<- VALUE, list(verbose=as.name(names(control)[i]), VALUE=control[[i]])))
            it <<- 0
        },

        run = function(runList) {
            if(!is.list(runList)) stop('runList argument should be a list')
            if(is.null(names(runList))) names(runList) <- rep('', length(runList))
            
            for(i in seq_along(runList)) {
                el <- runList[[i]]
                if(is.character(el)) {
                    type <- el
                } else if(is.list(el)) {
                    type <- 'blocks'
                    lst <- el
                    name <- if(names(runList)[i] == '') 'customBlocks' else names(runList)[i]
                } else if(class(el) == '{') {
                    type <- 'spec'
                    lst <- el
                    name <- if(names(runList)[i] == '') 'customSpec' else names(runList)[i]
                } else stop('dont understand element in run list')
                
                switch(type,
                       none    = runGroups(abModel$nodeGroupScalars,    'none'),
                       blocks  = runGroups(abModel$createGroups(lst),    name),
                       all     = runGroups(abModel$nodeGroupAllBlocked, 'all'),
                       default = runGroups(abModel$nodeGroupStochNodes, 'default', conjOveride=TRUE),
                       spec    = runGroups(lst,                          name, customSpec=TRUE),
                       auto    = runAutoBlock(),
                       stop('invalid argument'))
            }
        },

        runAutoBlock = function() {
            runGroups(abModel$nodeGroupScalars, 'auto0')
            groupsTagList <- determineGroupsFromPreviousSample()
            runGroups(groupsTagList$groups, paste0('auto1',groupsTagList$tag), printTree=TRUE)
            autoIt <- 1
            while(min(essPT[[it]]) > min(essPT[[it-1]])) {
                autoIt <- autoIt + 1
                groupsTagList <- determineGroupsFromPreviousSample()
                runGroups(groupsTagList$groups, paste0('auto',autoIt,groupsTagList$tag), printTree=TRUE)
            }
        },

        determineGroupsFromPreviousSample = function() {
            empCov[[it]] <<- cov(samples[[it]])
            empCovSparse[[it]] <<- empCov[[it]]
            if(sparsifyCov) {
                sparseSuccess <<- TRUE
                sparse_out <- try(spcov(Sigma=diag(diag(empCov[[it]])), S=empCov[[it]], lambda=spcov_lambda, step.size=spcov_step)$Sigma, silent=TRUE)
                if(inherits(sparse_out, 'try-error')) { sparseSuccess <<- FALSE
                                                    } else { empCovSparse[[it]] <<- sparse_out
                                                             dimnames(empCovSparse[[it]]) <<- dimnames(empCov[[it]]) }
            }
            empCovSparseThresh[[it]] <<- empCovSparse[[it]]
            if(sparseCovThreshold > 0) {
                ind <- abs(empCovSparseThresh[[it]]) > sparseCovThreshold
                diag(ind) <- true
                empCovSparseThresh[[it]] <<- empCovSparseThresh[[it]] * ind
            }
            empCor[[it]] <<- cov2cor(empCovSparseThresh[[it]])
            distMatrix[[it]] <<- as.dist(1 - abs(empCor[[it]]))
            hTree[[it]] <<- hclust(distMatrix[[it]])
            tag <- ''
            switch(cutree_method,
                   height = {
                       ct <- cutree(hTree[[it]], h = cutree_h)
                       groups <- lapply(unique(ct), function(x) names(ct)[ct==x]) },
                   dynamic = {
                       hybrid <- cutreeHybrid(dendro=hTree[[it]], distM=as.matrix(distMatrix[[it]]), minClusterSize=1, cutHeight=cutree_h, deepSplit=4, verbose=0)
                       ct <- hybrid$labels
                       names(ct) <- hTree[[it]]$labels
                       groups <- lapply(unique(ct), function(x) names(ct)[ct==x]) },
                   custom = {
                       if(!cutree_maxGroupSizeDynamic) {
                           groups <- cutree_custom(hTree[[it]], maxHeight=cutree_h, maxGroupSize=cutree_maxGroupSize, maxHeightRelativeFromBase=cutree_maxHeightRelativeFromBase)
                       } else {
                           nNodes <- length(abModel$nodeGroupScalars)
                           if(nNodes == 1) stop()
                           maxSizes <- unique(c(2^(1:floor(log(nNodes)/log(2))), nNodes))
                           Rmodel <- abModel$newModel()
                           groupsTmp <- specsTmp <- RmcmcsTmp <- CmcmcsTmp <- timingTmp <- samplesTmp <- essTmp <- essPTTmp <- essPTminTmp <-list()
                           for(iSize in seq_along(maxSizes)) {
                               groupsTmp[[iSize]] <- cutree_custom(hTree[[it]], maxHeight=cutree_h, maxGroupSize=maxSizes[iSize], maxHeightRelativeFromBase=cutree_maxHeightRelativeFromBase)
                               specsTmp[[iSize]] <- MCMCspec(Rmodel, nodes=NULL, monitors=character(0))
                               for(nodeGroup in groupsTmp[[iSize]]) addSamplerToSpec(Rmodel, specsTmp[[iSize]], nodeGroup, conjOveride=FALSE)
                               specsTmp[[iSize]]$addMonitors(Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE), print=FALSE)
                               RmcmcsTmp[[iSize]] <- buildMCMC(specsTmp[[iSize]])
                           }
                           Cmodel <- compileNimble(Rmodel)
                           CmcmcsTmp_temp <- compileNimble(RmcmcsTmp, project = Rmodel)
                           if(length(RmcmcsTmp) == 1) { CmcmcsTmp[[1]] <- CmcmcsTmp_temp
                           } else                     { CmcmcsTmp      <- CmcmcsTmp_temp }
                           if(setSeed0) set.seed(0)
                           for(iSize in seq_along(maxSizes)) {
                               timingTmp[[iSize]] <- as.numeric(system.time(CmcmcsTmp[[iSize]](niter))[3])
                               samplesTmp[[iSize]] <- as.matrix(nfVar(CmcmcsTmp[[iSize]], 'mvSamples'))
                               essTmp[[iSize]] <- apply(samplesTmp[[iSize]], 2, effectiveSize)
                               essPTTmp[[iSize]] <- essTmp[[iSize]]/timingTmp[[iSize]]
                               essPTminTmp[[iSize]] <- sort(essPTTmp[[iSize]])[1]
                           }
                           bestInd <- which(unlist(essPTminTmp) == max(unlist(essPTminTmp)))
                           groups <- groupsTmp[[bestInd]]
                           tag <- paste0('-max', maxSizes[bestInd])
                       }
                   },
                   stop('cutree method invalid'))
            return(list(groups=groups, tag=tag))
        },

        runGroups = function(groups, name, printTree=FALSE, conjOveride=FALSE, customSpec=FALSE) {
            if(setSeed0) set.seed(0)
            it <<- it + 1
            naming[[it]] <<- name
            Rmodel <- abModel$newModel()
            if(customSpec) {
                spec <- eval(groups, envir=environment())
                groups <- list()
                for(ss in spec$samplerSpecs) {
                    nodes <- if(!is.null(ss$control$targetNode)) ss$control$targetNode else ss$control$targetNodes
                    nodes <- Rmodel$expandNodeNames(nodes, returnScalarComponents=TRUE)
                    groups[[length(groups)+1]] <- nodes
                }
            } else {
                spec <- MCMCspec(Rmodel, nodes=NULL, monitors=character(0))
                for(nodeGroup in groups) addSamplerToSpec(Rmodel, spec, nodeGroup, conjOveride)
            }
            grouping[[it]] <<- lapply(groups, function(gp) Rmodel$expandNodeNames(gp, returnScalarComponents=TRUE))
            groupSizes[[it]] <<- numeric(0)
            for(gp in grouping[[it]]) for(node in gp) groupSizes[[it]][[node]] <<- if(node %in% names(groupSizes[[it]])) max(length(gp), groupSizes[[it]][[node]]) else length(gp)
            spec$addMonitors(Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE), print=FALSE)
            Rmcmc <- buildMCMC(spec)
            Cmodel <- compileNimble(Rmodel)
            Cmcmcs[[it]] <<- compileNimble(Rmcmc, project = Rmodel)
            timing[[it]] <<- as.numeric(system.time(Cmcmcs[[it]](niter))[3])
            samples[[it]] <<- as.matrix(nfVar(Cmcmcs[[it]], 'mvSamples'))
            ess[[it]] <<- round(apply(samples[[it]], 2, effectiveSize), 0)
            essPT[[it]] <<- round(ess[[it]]/timing[[it]], 1)
            essPT[[it]] <<- sort(essPT[[it]])
            updateNames()
            if(verbose) printCurrent(name, spec, printTree)
        },

        addSamplerToSpec = function(Rmodel, spec, nodeGroup, conjOveride) {
            if(length(nodeGroup) > 1) {
                spec$addSampler(type = 'RW_block', control = list(targetNodes=nodeGroup, adaptInterval=adaptIntervalBlock), print = FALSE); return()
            }
            if(!(nodeGroup %in% Rmodel$getNodeNames())) {
                spec$addSampler(type = 'RW', control = list(targetNode=nodeGroup, adaptInterval=adaptInterval), print = FALSE); return()
            }
            if(nodeGroup %in% Rmodel$getMaps('nodeNamesEnd')) {
                spec$addSampler(type = 'end', control = list(targetNode=nodeGroup), print = FALSE); return()
            }
            conjugacyResult <- Rmodel$checkConjugacy(nodeGroup)
            if((!is.null(conjugacyResult)) && (conjOveride || useConjugacy)) {
                spec$addSampler(type = conjugacyResult$samplerType, control = conjugacyResult$control, print = FALSE); return()
            }
            if(Rmodel$getNodeInfo()[[nodeGroup]]$isDiscrete()) {
                spec$addSampler(type = 'slice', control = list(targetNode=nodeGroup), print = FALSE); return()
            }
            if(length(Rmodel$expandNodeNames(nodeGroup, returnScalarComponents = TRUE)) > 1) {
                spec$addSampler(type = 'RW_block', control = list(targetNodes=nodeGroup, adaptInterval=adaptIntervalBlock), print = FALSE); return()
            }
            spec$addSampler(type = 'RW', control = list(targetNode=nodeGroup, adaptInterval=adaptInterval), print = FALSE); return()
        },

        updateNames = function() {
            names(grouping) <<- naming
            names(groupSizes) <<- naming
            names(Cmcmcs) <<- naming
            names(timing) <<- naming
            names(samples) <<- naming
            names(ess) <<- naming
            names(essPT) <<- naming
        },

        printCurrent = function(name, spec, printTree) {
            cat(paste0('\n######################\nITERATION ', it, ': ', name, '\n######################\n'))
            if(printTree && sparsifyCov) {
                if(sparseSuccess) cat('\nsparsifying empirical covariance matrix.....\n') else cat('\nsparsification failed, covariance matrix nearly singular\n') }
            cat('\ngroups:\n'); g <- grouping[[it]]
            if(is.list(g)) { for(i in seq_along(g)) cat(paste0('[', i, '] ', paste0(g[[i]], collapse=', '), '\n'))
                         } else { cat('custom MCMC specification\n') }
            if(printTree) { dev.new(); if(inherits(try(plot(as.dendrogram(hTree[[it-1]]), ylim=c(0,1), main=name), silent=TRUE), 'try-error')) dev.off() }
            cat('\nsamplers:\n'); spec$getSamplers()
            cat(paste0('\nMCMC runtime: ', round(timing[[it]], 2), ' seconds\n'))
            cat('\nESS:\n'); print(ess[[it]])
            cat('\nESS/time:\n'); print(essPT[[it]])
            cat(paste0('\n######################\nITERATION ', it, ': ', name, '\n######################\n'))
            cat('\n')
        }
        )
    )



cutree_custom <- function(ht, maxHeight, maxGroupSize, maxHeightRelativeFromBase) {
    labels <- ht$labels;     height <- ht$height;     merge <- ht$merge
    nNodes <- length(labels)
    nMerges <- dim(merge)[1]
    if(nMerges+1 != nNodes) stop('something fishy')
    df <- data.frame(baseHeight=numeric(nNodes), recentHeight = numeric(nNodes), num=numeric(nNodes))
    for(i in 1:nMerges) { for(j in 1:2) {
        if(merge[i,j] < 0) df$baseHeight[abs(merge[i,j])] <- df$recentHeight[abs(merge[i,j])] <- height[i] } }
    df$num <- rep(1, nNodes)
    nodes <- lapply(1:nNodes, function(x) x)
    for(i in 1:nMerges) {
        ind1 <- if(merge[i,1] < 0) abs(merge[i,1]) else (merge[i,1]+nNodes)
        ind2 <- if(merge[i,2] < 0) abs(merge[i,2]) else (merge[i,2]+nNodes)
        indout <- i + nNodes
        df[indout,] <- NA
        if((df$num[ind1]==0) || (df$num[ind2]==0) || is.na(df$num[ind1]) || is.na(df$num[ind2])) stop('something wrong')
        canMerge <- TRUE
        if(height[i] > maxHeight) canMerge <- FALSE
        if(df$num[ind1] + df$num[ind2] > maxGroupSize) {
##            rh <- max(df$recentHeight[c(ind1,ind2)])
##            if(height[i] > rh + maxGroupSizeRelHeightOveride*(1-rh)) canMerge <- FALSE
            canMerge <- FALSE
        }
        bh <- min(df$baseHeight[c(ind1,ind2)])
        if(height[i] > bh + maxHeightRelativeFromBase*(1-bh)) canMerge <- FALSE
        if(canMerge) {
            ## merge
            df$baseHeight[indout] <- min(df$baseHeight[c(ind1,ind2)])
            df$recentHeight[indout] <- height[i]
            df$num[indout] <- df$num[ind1] + df$num[ind2]
            nodes[[indout]] <- c(nodes[[ind2]], nodes[[ind1]])
            if(df$num[indout] != length(nodes[[indout]])) stop('something wrong')
            df[ind1,] <- NA;     df[ind2,] <- NA
            df$num[c(ind1,ind2)] <- c(0,0)
            nodes[[ind1]] <- NA;     nodes[[ind2]] <- NA
        } else {
            ## don't merge; lower branch becomes a group, higher branch propagates up
            indhigher <- if(df$recentHeight[ind1] > df$recentHeight[ind2]) ind1 else ind2
            df[indout,] <- df[indhigher,]
            df[indhigher,] <- NA
            df$num[indhigher] <- 0
            nodes[[indout]] <- nodes[[indhigher]]
            nodes[[indhigher]] <- NA
        }
    }
    groupInd <- which(df$num > 0)
    groupNodeNumbers <- nodes[groupInd]
    groups <- lapply(groupNodeNumbers, function(gnn) labels[gnn])
    return(groups)
}



createDFfromABlist <- function(lst, rho=NA, rho2=NA) {
    df <- data.frame(model=character(), rho=numeric(), rho2=numeric(), blocking=character(), timing=numeric(), node=character(), groupSize=numeric(), ess=numeric(), essPT=numeric(), stringsAsFactors=FALSE)
    for(iAB in seq_along(lst)) {
        ab <- lst[[iAB]]
        abName <- names(lst)[iAB]
        for(iBlock in seq_along(ab$naming)) {
            blocking <- ab$naming[[iBlock]]
            timing <- ab$timing[[iBlock]]
            ess <- ab$ess[[iBlock]]
            nodes <- names(ess)
            nNodes <- 
            essPT <- ab$essPT[[iBlock]][nodes]            ## sort
            groupSizes <- ab$groupSizes[[iBlock]][nodes]  ## sort
            newIndDF <- (1:length(nodes)) + dim(df)[1]
            df[newIndDF,] <- NA
            df[newIndDF,]$model <- abName
            df[newIndDF,]$blocking <- blocking
            df[newIndDF,]$timing <- timing
            df[newIndDF,]$node <- nodes
            df[newIndDF,]$groupSize <- groupSizes
            df[newIndDF,]$ess <- ess
            df[newIndDF,]$essPT <- essPT
        }
    }
    df$rho  <- rho
    df$rho2 <- rho2
    return(df)
}



plotABS <- function(df, plotGroupSizes=TRUE, xlimToMin=FALSE, together) {
    models <- unique(df$model)
    nModels <- length(models)
    if(missing(together)) together <- if(nModels <= 5) TRUE else FALSE
    nVertPlots <- if(together) nModels*2 else nModels
    xVarNames <- c('ess', 'essPT')
    parCmd <- quote(par(mfrow=c(nVertPlots,1),mar=c(1,0,1,0),tcl=-.1,mgp=c(3,0,0),cex.axis=.7))
    if(together) { quartz(); eval(parCmd) }
    for(xVarName in xVarNames) {
        if(!together) { quartz(); eval(parCmd) }
        maxMinXVar<-0; for(mod in models) {dfMod<-df[df$model==mod,]; blks<-unique(dfMod$blocking); for(blk in blks) {maxMinXVar<-max(maxMinXVar,min(dfMod[dfMod$blocking==blk,xVarName]))}}
        maxXVar <- if(xlimToMin) maxMinXVar else max(df[, xVarName])
        xlim <- c(maxXVar*-0.05, maxXVar)
        maxTiming <- max(df[, 'timing'])
        for(mod in models) {
            dfMod <- df[df$model==mod,]
            blockings <- unique(dfMod$blocking)
            nBlockings <- length(blockings)
            bestBlk<-''; bestEssPT<-0; for(blk in blockings) { if(min(dfMod[dfMod$blocking==blk,'essPT'])>bestEssPT && ((blk=='all')||grepl('^auto',blk))) {bestEssPT<-min(dfMod[dfMod$blocking==blk,'essPT']); bestBlk<-blk} }
            plot(-100,-100,xlim=xlim,ylim=c(0,nBlockings+1),xlab='',ylab='',main=paste0(xVarName, ' for model ', mod))
            for(iBlocking in 1:nBlockings) {
                blocking <- blockings[iBlocking]
                dfModBlock <- dfMod[dfMod$blocking==blocking,]
                xVarValues <- dfModBlock[,xVarName]
                groupSizes <- dfModBlock[,'groupSize']
                timing <- dfModBlock[,'timing'][1]   # only first element
                timingOnXaxis <- timing/maxTiming * xlim[2]
                yCoord <- nBlockings+1-iBlocking
                lines(x=c(0,timingOnXaxis), y=rep(yCoord,2), lty=1, lwd=2, col='lightgrey')
                if(plotGroupSizes) { text(x=xVarValues, y=yCoord, labels=groupSizes, cex=0.7)
                } else { points(x=xVarValues, y=rep(yCoord,length(xVarValues)), pch=20) }
                col <- if(blocking == bestBlk) 'green' else 'blue'
                text(x=xlim[1], y=yCoord, labels=blocking, col=col)
                if(timing==maxTiming) text(xlim[2], yCoord+1, paste0('t = ',round(timing,1)))
            }
        }
    }
}


formatForFile <- function(code) {
    a <- deparse(code, width.cutoff=500L)
    a <- a[c(-1, -length(a))]
    a <- sub('^    ', '', a)
    a[length(a) + 1] <- ''
    a <- paste0(a, collapse='\n')
    a
}




########################################################
#######################################################
### models after here
########################################################

library(nimble)


################
### tester
################

## code <- modelCode({
##     x[1:5] ~ dmnorm(mu[1:5], cov=Sigma[1:5,1:5])
##     y[1:5] ~ dmnorm(x[1:5], cov=Ident[1:5,1:5])
## })

## mu <- rep(0,5)

## Chol <- matrix(0, nrow=5, ncol=5)
## Chol[1,1:5] <- c(1,   .99, .9, -.8, 0)
## Chol[2,1:5] <- c(0, 1,   .8, 0,   0)
## Chol[3,1:5] <- c(0,  0,  1,  -.9, 0)
## Chol[4,1:5] <- c(0, 0,  0,  1,   0)
## Chol[5,1:5] <- c(0,   0,   0,  0,   1)
## Sigma <- t(Chol) %*% Chol
## Prec <- solve(Sigma)
## cov2cor(Sigma)
## Ident <- diag(5)
## y <- c(3, 3, 2, 0, -4)
## postPrec <- Prec + Ident
## postSigma <- solve(postPrec)
## postCor <- cov2cor(postSigma)
## postMean <- postSigma %*% t(t(y))

## constants <- list(mu=mu, Sigma=Sigma, Ident=Ident)
## data <- list(y = y)
## inits <- list()

## abtester <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)


################
### litters
################

G <- 2
N <- 16
n <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13, 7, 5, 10, 7, 6, 10, 10, 10, 7), dim = c(2, 16))
r <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9, 12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0), dim = c(2, 16))
p <- array(0.5, dim = c(2, 16))

constants_litters<- list(G=G, N=N, n=n)
data_litters      <- list(r=r)
inits_litters     <- list(p=p)

code_litters <- modelCode({
     for (i in 1:G) {
         a[i] ~ dgamma(1, 0.001)
         b[i] ~ dgamma(1, 0.001)
         # mu[i] <- a[i] / (a[i] + b[i])
         # theta[i] <- 1 / (a[i] + b[i])
         for (j in 1:N) {
             r[i,j] ~ dbin(p[i,j], n[i,j])
             p[i,j] ~ dbeta(a[i], b[i])
         }
     }
})

rm(list = c('G','N','n','r','p'))

################
### SSMmub
################

## better parameterization: mean and autocorrelation
code_SSMmub<- modelCode({
    mu ~ dnorm(0, sd = 1000)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(mu, sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    a <- 1-(b/mu)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})

t <- 30
constants_SSMmub <- list(t = t)
Rmodel <- nimbleModel(code_SSMmub, constants = constants_SSMmub)
Rmodel$mu <- 10/(1-.5)
Rmodel$b <- 10
Rmodel$sigPN <- .1
Rmodel$sigOE <- .1
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data_SSMmub <- list(y = Rmodel$y)
inits_SSMmub <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)

rm(list = c('t','Rmodel'))

################
### SSMab
################

## parameterization in terms of slope (autocorelation) and intercept,
## which are highly correlated in mixing
code_SSMab <- modelCode({
    a ~ dunif(-0.9999, 0.9999)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(b/(1-a), sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})

t <- 30
constants_SSMab <- list(t = t)
Rmodel <- nimbleModel(code_SSMab, constants = constants_SSMab)
Rmodel$a <- .5
Rmodel$b <- 10
Rmodel$sigPN <- .1
Rmodel$sigOE <- .1
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('a','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data_SSMab <- list(y = Rmodel$y)
inits_SSMab <- list(a = Rmodel$a, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)

rm(list = c('t','Rmodel'))
